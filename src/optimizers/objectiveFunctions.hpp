#ifndef ULOCATOR_OPTIMIZERS_OBJECTIVE_FUNCTIONS_HPP
#define ULOCATOR_OPTIMIZERS_OBJECTIVE_FUNCTIONS_HPP
#include <iostream>
#include <cmath>
#include <array>
#include <string>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/station.hpp"

namespace
{

enum class Norm
{
    L1,
    LeastSquares,
    Lp
};

enum class Measurement
{
    Standard,
    DoubleDifference
};

template<typename T>
T sgn(const T x)
{
    constexpr T zero{0};
    if (x == zero){return 0;} 
    if (x > zero)
    {   
        return 1;
    }   
    else //if (x < zero)
    {   
        return -1; 
    }   
}

template<typename T>
[[nodiscard]] 
T leastSquares(const std::vector<T> &weights,
               const std::vector<T> &observations,
               const std::vector<T> &estimates,
               const Measurement measurement)
{
     const T *__restrict__ weightsPtr = weights.data();
     const T *__restrict__ observationsPtr = observations.data();
     const T *__restrict__ estimatesPtr = estimates.data();
     auto n = static_cast<int> (weights.size());
     double sumSquaredOfResiduals{0};
     if (measurement == Measurement::Standard)
     {
         for (int i = 0; i < n; ++i)
         {
             double weightedResidual
                = weightsPtr[i]*(observationsPtr[i] - estimatesPtr[i]);
             sumSquaredOfResiduals = sumSquaredOfResiduals
                                   + weightedResidual*weightedResidual;
         }
     }
     else if (measurement == Measurement::DoubleDifference)
     {   
         for (int i = 0; i < n; ++i)
         {
             for (int j = i + 1; j < n; ++j)
             {
                 double doubleDifferenceResidual
                    = (observationsPtr[i] - observationsPtr[j])
                    - (estimatesPtr[i]    - estimatesPtr[j]);
                 double weight = std::sqrt(weightsPtr[i]*weightsPtr[j]);
                 double weightedResidual = weight*doubleDifferenceResidual;
                 sumSquaredOfResiduals = sumSquaredOfResiduals
                                       + weightedResidual*weightedResidual;
             }
         }
     }
     return static_cast<T> (sumSquaredOfResiduals);
}

template<typename T>
[[nodiscard]]
T l1(const std::vector<T> &weights,
     const std::vector<T> &observations,
     const std::vector<T> &estimates,
     const Measurement measurement)
{
     const T *__restrict__ weightsPtr = weights.data();
     const T *__restrict__ observationsPtr = observations.data();
     const T *__restrict__ estimatesPtr = estimates.data();
     auto n = static_cast<int> (weights.size());
     double sumAbsoluteResiduals{0};
     if (measurement == Measurement::Standard)
     {
         for (int i = 0; i < n; ++i)
         {
             double weightedAbsoluteResidual
                = weightsPtr[i]*std::abs(observationsPtr[i] - estimatesPtr[i]);
             sumAbsoluteResiduals = sumAbsoluteResiduals
                                  + weightedAbsoluteResidual;
         }
     }
     else if (measurement == Measurement::DoubleDifference)
     {
         for (int i = 0; i < n; ++i)
         {
             for (int j = i + 1; j < n; ++j)
             {
                 double doubleDifference 
                    = (observationsPtr[i] - observationsPtr[j])
                    - (estimatesPtr[i]    - estimatesPtr[j]);
                 double weight = std::sqrt(weightsPtr[i]*weightsPtr[j]);
                 double weightedAbsoluteResidual
                    = weight*std::abs(doubleDifference);
                 sumAbsoluteResiduals = sumAbsoluteResiduals
                                      + weightedAbsoluteResidual;
             }
         }
     }
#ifndef NDEBUG
     else
     {
         assert(false);
     }
#endif
     return static_cast<T> (sumAbsoluteResiduals);
}

template<typename T>
[[nodiscard]] 
T lp(const std::vector<T> &weights,
     const std::vector<T> &observations,
     const std::vector<T> &estimates,
     const double p,
     const Measurement measurement)
{
#ifndef NDEBUG
     assert(p >= 1);
#endif
     const T *__restrict__ weightsPtr = weights.data();
     const T *__restrict__ observationsPtr = observations.data();
     const T *__restrict__ estimatesPtr = estimates.data();
     auto n = static_cast<int> (weights.size());
     double sumAbsoluteResidualsToP{0};
     if (measurement == Measurement::Standard)
     { 
         for (int i = 0; i < n; ++i)
         {
             double weightedAbsoluteResidual
                = weightsPtr[i]*std::abs(observationsPtr[i] - estimatesPtr[i]);
             sumAbsoluteResidualsToP = sumAbsoluteResidualsToP
                                     + std::pow(weightedAbsoluteResidual, p);
         }
     }
     else if (measurement == Measurement::DoubleDifference)
     {
         for (int i = 0; i < n; ++i)
         {
             for (int j = i + 1; j < n; ++j)
             {
                 double doubleDifference 
                     = (observationsPtr[i] - observationsPtr[j])
                     - (estimatesPtr[i]    - estimatesPtr[j]);
                 double weight = std::sqrt(weightsPtr[i]*weightsPtr[j]);
                 double weightedAbsoluteResidual
                     = weight*std::abs(doubleDifference);
                 sumAbsoluteResidualsToP
                     = sumAbsoluteResidualsToP
                     + std::pow(weightedAbsoluteResidual, p);
             }
         }
     }
#ifndef NDEBUG
     else
     {
         assert(false);
     }
#endif
     return std::pow(sumAbsoluteResidualsToP, 1./p);;
}

//----------------------------------------------------------------------------//
//                                 3D + Time                                  //
//----------------------------------------------------------------------------//
/// @brief Cost function for the general problem of an unknown hypocenter
///        and origin time.
struct Problem3DAndTime
{
    /// @result The log-likelihood function for maximizing something of
    ///         the form: LL = -C(observed, estimates) - penalty
    /// @param[in]
    double logLikelihood(const int n, const double x[]) const
    {
        return logLikelihoodAndGradient(n, x, nullptr);
    }
    /// @param[out] gradient  A finite difference approximation of the gradient.
    /// @result The log-likelihood function evaluated at x.
    /// @note This is for testing purposes only.
    double finiteDifference(const int n, const double x[],
                            double gradient[])
    {
#ifndef NDEBUG
        assert(n == 4);
        assert(gradient != nullptr);
#endif
        std::array<double, 4> xCopy;
        std::array<double, 4> perturbation{0.0001,  // This is about 1.e-14 for
                                                    // origin times.  Any
                                                    // smaller is dangerous.
                                           0.001,   // This is 1 mm
                                           0.001,   // This is 1 mm
                                           0.001    // This is 1 mm
                                           };
        auto f = logLikelihood(n, x);
        for (int i = 0; i < n; ++i)
        {
            std::copy(x, x + n, xCopy.data());
            xCopy[i] = xCopy[i] + perturbation[i];
            auto fPert = logLikelihood(xCopy.size(), xCopy.data());
            gradient[i] = (fPert - f)/perturbation[i];
        }
        return f;
    }
    /// @brief The actual workhorse that computes the log-likelihood function
    ///        and, optionally, the gradient.
    /// @param[in] n   The number of parameters to optimize.
    /// @param[in] x   The current position of the solution vector.  This
    ///                is an array with dimension [n].
    /// @param[in,out] gradient  If not a nullptr then this 
    /// @result The objective function and the inequality constraint that 
    ///         requires the source be below the topography. 
    double logLikelihoodAndGradient(
        const int n, const double x[], double gradient[]) const
    {
#ifndef NDEBUG
        assert(n == 4);
        assert(x != nullptr);
#endif
        bool doGradient{false};
        if (gradient != nullptr){doGradient = true;}
        double originTime = x[0];
        double xSource = x[1];
        double ySource = x[2];
        double zSource = x[3];
        std::vector<double> arrivalTimes(mStationPhases.size());
        double costFunction{0};
        if (doGradient)
        {
            std::fill(gradient, gradient + n, 0);
            auto nArrivals = static_cast<int> (mStationPhases.size());
            std::vector<double> dtdt0s(nArrivals);
            std::vector<double> dtdxs(nArrivals);
            std::vector<double> dtdys(nArrivals);
            std::vector<double> dtdzs(nArrivals);
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               &dtdt0s, &dtdxs, &dtdys, &dtdzs,
                                               mApplyCorrection);
            if (mNorm == ::Norm::LeastSquares)
            {
                constexpr double two{2};
                for (int i = 0; i < nArrivals; ++i)
                {
                    double residual = mObservations[i] - arrivalTimes[i];
                    double twoWeightSquared = two*(mWeights[i]*mWeights[i]);
                    double twoWeightSquaredResidual = twoWeightSquared*residual;
                    gradient[0] = gradient[0]
                                - twoWeightSquaredResidual*dtdt0s[i];
                    gradient[1] = gradient[1]
                                - twoWeightSquaredResidual*dtdxs[i];
                    gradient[2] = gradient[2]
                                - twoWeightSquaredResidual*dtdys[i];
                    gradient[3] = gradient[3]
                                - twoWeightSquaredResidual*dtdzs[i];
                }
            }
            else if (mNorm == ::Norm::L1)
            {
                for (int i = 0; i < nArrivals; ++i)
                {
                    double weightedResidual
                        = mWeights[i]*(mObservations[i] - arrivalTimes[i]);
                    auto sgnWeightedResidual = ::sgn(weightedResidual);
                    double weightedSgnWeightedResidual
                        = mWeights[i]*sgnWeightedResidual;
                    gradient[0] = gradient[0]
                                - weightedSgnWeightedResidual*dtdt0s[i];
                    gradient[1] = gradient[1]
                                - weightedSgnWeightedResidual*dtdxs[i];
                    gradient[2] = gradient[2]
                                - weightedSgnWeightedResidual*dtdys[i];
                    gradient[3] = gradient[3]
                                - weightedSgnWeightedResidual*dtdzs[i];
                }
            }
            else if (mNorm == ::Norm::Lp)
            {
                // p-norm: ||w*(T - T(x))||_p
                //       = ( sum_i |w_i*(T_i - T_i(x))|^p )^(1/p)
                // Chain rule:
                //       = 1/p (( sum_i |w_i*(T_i - T_i(x))|^p )^(-(p - 1)/p)
                //        *p*sum_j | w_j*(T_j - T_j(x))|^(p - 1) delta_{i,j}
                //        *sgn(w_i*(T_i - T_i(x)))
                //        *-w_i dT/dx
                //       =-(|w_i*(T_i - T_i(x))|/||w*(T - T(x))||_p))^(p - 1)
                //        w_i*sgn(w_i*(T_i - T_i(x)))*dT_i/dx
                // I need the cost function as a normalization term
                double costFunctionInverse
                    = ::lp(mWeights, mObservations, arrivalTimes,
                           mPNorm, Measurement::Standard);
                if (costFunctionInverse > 0)
                {
                    costFunctionInverse = 1./costFunctionInverse;
                } 
                for (int i = 0; i < nArrivals; ++i)
                {
                    double weightedResidual
                        = mWeights[i]*(mObservations[i] - arrivalTimes[i]);
                    double weightedAbsoluteResidual
                        = std::abs(weightedResidual);
                    double argument
                        = weightedAbsoluteResidual*costFunctionInverse; 
                    double scalar
                        =-(mWeights[i]*std::pow(argument, mPNorm - 1))
                         *::sgn(weightedResidual);
                    gradient[0] = gradient[0] + scalar*dtdt0s[i];
                    gradient[1] = gradient[1] + scalar*dtdxs[i];
                    gradient[2] = gradient[2] + scalar*dtdys[i];
                    gradient[3] = gradient[3] + scalar*dtdzs[i];
                }
            }
#ifndef NDEBUG
            else
            {
                assert(false);
            }
#endif
        }
        else
        {
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               mApplyCorrection);
        }
        if (mNorm == ::Norm::LeastSquares)
        {
            costFunction
                = ::leastSquares(mWeights, mObservations, arrivalTimes,
                                 Measurement::Standard);
        }
        else if (mNorm == ::Norm::L1)
        {
            costFunction
                = ::l1(mWeights, mObservations, arrivalTimes,
                       Measurement::Standard);
        }
        else if (mNorm == ::Norm::Lp)
        {
            costFunction
                = ::lp(mWeights, mObservations, arrivalTimes,
                       mPNorm, Measurement::Standard);
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        // Use a penalty formulation for topography.  So we have
        //    L(x,y,z) = E(x,y) - zSource when zSource > E(x,y) and 0 otherwise
        // so the partials are as follows:
        //    dL/dx = dE/dx
        //    dL/dy = dE/dy
        //    dL/dz =-1
        double elevationPenalty{0};
        if (doGradient)
        {
            double dEdx, dEdy;
            double elevation = mTopography->evaluate(xSource, ySource,
                                                     &dEdx, &dEdy);
            elevation =-elevation; // Elevation increases +up; make it +down
            dEdx =-dEdx;
            dEdy =-dEdy;
            if (zSource < elevation)
            {
                elevationPenalty = elevation - zSource;
                gradient[1] = gradient[1] + mPenaltyCoefficient*dEdx;
                gradient[2] = gradient[2] + mPenaltyCoefficient*dEdy;
                gradient[3] = gradient[3] - 1*mPenaltyCoefficient;
            }
        }
        else
        {
            double elevation = mTopography->evaluate(xSource, ySource);
            elevation =-elevation; // Elevation increases +up; make it +down
            if (zSource < elevation)
            {
                elevationPenalty = elevation - zSource;
            }
        }
#ifndef NDEBUG
        assert(elevationPenalty >= 0);
#endif
        // Maximizing the posterior probability amounts to minimizing something
        // of the form exp(-arg).  The log is simply -arg.  Additionally,
        // the arg = Cost + Penalty (i.e., the goal is to make small the
        // cost AND the penalty)  
        double logLikelihood
             =-(costFunction + mPenaltyCoefficient*elevationPenalty);
        if (doGradient)
        {
            std::transform(gradient, gradient + n, gradient,
                           [](const auto &gi)
                           {
                               return -gi;
                           });
        } 
        return logLikelihood;
    }
    /// @result The hard bounds on the search region.  Most optimization tools
    ///         can handle linear constraints.  Otherwise, I could add a penalty 
    ///         to the objective function that requires the origin time be
    ///         less than the first arrival time.
    std::pair<std::vector<double>, std::vector<double>> getLinearBounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    const ULocator::TravelTimeCalculatorMap *mTravelTimeCalculatorMap{nullptr};
    const ULocator::Topography::ITopography *mTopography{nullptr};
    std::vector<std::pair<ULocator::Station, std::string>> mStationPhases;
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    ::Norm mNorm{Norm::LeastSquares};
    double mPNorm{1.5};
    double mPenaltyCoefficient{1};
    bool mApplyCorrection{true};
};

//----------------------------------------------------------------------------//
//                            Fixed Depth + Time                              //
//----------------------------------------------------------------------------//
/// @brief Cost function for the problem of a an unknown epicenter and origin
///        time but specified depth.
struct Problem2DAndTimeAndFixedDepth
{
    /// @result The log-likelihood function for maximizing something of
    ///         the form: LL = -C(observed, estimates)
    /// @param[in]
    double logLikelihood(const int n, const double x[]) const
    {
        return logLikelihoodAndGradient(n, x, nullptr);
    }
    /// @param[out] gradient  A finite difference approximation of the gradient.
    /// @result The log-likelihood function evaluated at x.
    /// @note This is for testing purposes only.
    double finiteDifference(const int n, const double x[],
                            double gradient[])
    {
#ifndef NDEBUG
        assert(n == 3);
        assert(gradient != nullptr);
#endif
        std::array<double, 3> xCopy;
        std::array<double, 3> perturbation{0.0001,  // Origin time
                                           0.001,   // This is 1 mm
                                           0.001    // This is 1 mm
                                           };
        auto f = logLikelihood(n, x);
        for (int i = 0; i < n; ++i)
        {
            std::copy(x, x + n, xCopy.data());
            xCopy[i] = xCopy[i] + perturbation[i];
            auto fPert = logLikelihood(xCopy.size(), xCopy.data());
            gradient[i] = (fPert - f)/perturbation[i];
        }
        return f;
    }
    /// @brief The actual workhorse that computes the log-likelihood function
    ///        and, optionally, the gradient.
    /// @param[in] n   The number of parameters to optimize.
    /// @param[in] x   The current position of the solution vector.  This
    ///                is an array with dimension [n].
    /// @param[in,out] gradient  If not a nullptr then this 
    /// @result The objective function and the inequality constraint that 
    ///         requires the source be below the topography. 
    double logLikelihoodAndGradient(
        const int n, const double x[], double gradient[]) const
    {
#ifndef NDEBUG
        assert(n == 3);
        assert(x != nullptr);
#endif
        bool doGradient{false};
        if (gradient != nullptr){doGradient = true;}
        double originTime = x[0];
        double xSource = x[1];
        double ySource = x[2];
        double zSource = mDepth;
        std::vector<double> arrivalTimes(mStationPhases.size());
        double costFunction{0};
        if (doGradient)
        {
            std::fill(gradient, gradient + n, 0); 
            auto nArrivals = static_cast<int> (mStationPhases.size());
            std::vector<double> dtdt0s(nArrivals);
            std::vector<double> dtdxs(nArrivals);
            std::vector<double> dtdys(nArrivals);
            std::vector<double> dtdzs(nArrivals);
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               &dtdt0s, &dtdxs, &dtdys, &dtdzs,
                                               mApplyCorrection);
            if (mNorm == ::Norm::LeastSquares)
            {
                constexpr double two{2};
                for (int i = 0; i < nArrivals; ++i)
                {
                    double residual = mObservations[i] - arrivalTimes[i];
                    double twoWeightSquared = two*(mWeights[i]*mWeights[i]);
                    double twoWeightSquaredResidual = twoWeightSquared*residual;
                    gradient[0] = gradient[0]
                                - twoWeightSquaredResidual*dtdt0s[i];
                    gradient[1] = gradient[1]
                                - twoWeightSquaredResidual*dtdxs[i];
                    gradient[2] = gradient[2]
                                - twoWeightSquaredResidual*dtdys[i];
                }
            }
            else if (mNorm == ::Norm::L1)
            {
                for (int i = 0; i < nArrivals; ++i)
                {
                    double weightedResidual
                        = mWeights[i]*(mObservations[i] - arrivalTimes[i]);
                    auto sgnWeightedResidual = ::sgn(weightedResidual);
                    double weightedSgnWeightedResidual
                        = mWeights[i]*sgnWeightedResidual;
                    gradient[0] = gradient[0]
                                - weightedSgnWeightedResidual*dtdt0s[i];
                    gradient[1] = gradient[1]
                                - weightedSgnWeightedResidual*dtdxs[i];
                    gradient[2] = gradient[2]
                                - weightedSgnWeightedResidual*dtdys[i];
                }
            }
            else if (mNorm == ::Norm::Lp)
            {
                double costFunctionInverse
                    = ::lp(mWeights, mObservations, arrivalTimes,
                           mPNorm, Measurement::Standard);
                if (costFunctionInverse > 0)
                {
                    costFunctionInverse = 1./costFunctionInverse;
                }
                for (int i = 0; i < nArrivals; ++i)
                {
                    double weightedResidual
                        = mWeights[i]*(mObservations[i] - arrivalTimes[i]);
                    double weightedAbsoluteResidual
                        = std::abs(weightedResidual);
                    double argument
                        = weightedAbsoluteResidual*costFunctionInverse;
                    double scalar
                        =-(mWeights[i]*std::pow(argument, mPNorm - 1))
                         *::sgn(weightedResidual);
                    gradient[0] = gradient[0] + scalar*dtdt0s[i];
                    gradient[1] = gradient[1] + scalar*dtdxs[i];
                    gradient[2] = gradient[2] + scalar*dtdys[i];
                }
            }
#ifndef NDEBUG
            else
            {
                assert(false);
            }
#endif
        }
        else
        {
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               mApplyCorrection);
        }
        if (mNorm == ::Norm::LeastSquares)
        {
            costFunction
                = ::leastSquares(mWeights, mObservations, arrivalTimes,
                                 Measurement::Standard);
        }
        else if (mNorm == ::Norm::L1)
        {
            costFunction
                = ::l1(mWeights, mObservations, arrivalTimes,
                       Measurement::Standard);
        }
        else if (mNorm == ::Norm::Lp)
        {
            costFunction
                = ::lp(mWeights, mObservations, arrivalTimes,
                       mPNorm, Measurement::Standard);
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        double logLikelihood =-costFunction;
        if (doGradient)
        {
            std::transform(gradient, gradient + n, gradient,
                           [](const auto &gi)
                           {
                               return -gi;
                           });
        }
        return logLikelihood;
    }
    /// @result The hard bounds on the search region.  Most optimization tools
    ///         can handle linear constraints.
    std::pair<std::vector<double>, std::vector<double>> getLinearBounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    const ULocator::TravelTimeCalculatorMap *mTravelTimeCalculatorMap{nullptr};
    std::vector<std::pair<ULocator::Station, std::string>> mStationPhases;
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    ::Norm mNorm{Norm::LeastSquares};
    double mDepth{6000};
    double mPNorm{1.5};
    double mPenaltyCoefficient{1};
    bool mApplyCorrection{true};
};

//----------------------------------------------------------------------------//
//                   2D + Time + Depth at Free Surface                        //
//----------------------------------------------------------------------------//
/// @brief Cost function for the problem of unknown epicenter and origin time
///        but with the depth fixed to free surface.
struct Problem2DAndTimeAndDepthAtFreeSurface
{
    /// @result The log-likelihood function for maximizing something of
    ///         the form: LL = -C(observed, estimates) - penalty
    /// @param[in]
    double logLikelihood(const int n, const double x[]) const
    {
        return logLikelihoodAndGradient(n, x, nullptr);
    }
    /// @param[out] gradient  A finite difference approximation of the gradient.
    /// @result The log-likelihood function evaluated at x.
    /// @note This is for testing purposes only.
    double finiteDifference(const int n, const double x[],
                            double gradient[])
    {
#ifndef NDEBUG
        assert(n == 3);
        assert(gradient != nullptr);
#endif
        std::array<double, 3> xCopy;
        std::array<double, 3> perturbation{0.0001,  // This is 0.1 millisec
                                           0.001,   // This is 1 mm
                                           0.001    // This is 1 mm
                                           };
        auto f = logLikelihood(n, x);
        for (int i = 0; i < n; ++i)
        {
            std::copy(x, x + n, xCopy.data());
            xCopy[i] = xCopy[i] + perturbation[i];
            auto fPert = logLikelihood(xCopy.size(), xCopy.data());
            gradient[i] = (fPert - f)/perturbation[i];
        }
        return f;
    }
    /// @brief The actual workhorse that computes the log-likelihood function
    ///        and, optionally, the gradient.
    /// @param[in] n   The number of parameters to optimize.
    /// @param[in] x   The current position of the solution vector.  This
    ///                is an array with dimension [n].
    /// @param[in,out] gradient  If not a nullptr then this 
    /// @result The objective function and the inequality constraint that 
    ///         requires the source be below the topography. 
    double logLikelihoodAndGradient(
        const int n, const double x[], double gradient[]) const
    {
#ifndef NDEBUG
        assert(n == 3);
        assert(x != nullptr);
#endif
        bool doGradient{false};
        if (gradient != nullptr){doGradient = true;}
        double originTime = x[0];
        double xSource = x[1];
        double ySource = x[2];
        // Source is fixed to elevation.  Note, sign convention for
        // elevation is +up whereas solver is +down hence the minus signs.
        double zSource{0};
        double dEdx{0}, dEdy{0};
        if (doGradient)
        {
            zSource =-mTopography->evaluate(xSource, ySource, &dEdx, &dEdy);
            dEdx =-dEdx;
            dEdy =-dEdy;
        }
        else
        {
            zSource =-mTopography->evaluate(xSource, ySource);
        }
        std::vector<double> arrivalTimes(mStationPhases.size());
        double costFunction{0};
        if (doGradient)
        {
            std::fill(gradient, gradient + n, 0);
            auto nArrivals = static_cast<int> (mStationPhases.size());
            std::vector<double> dtdt0s(nArrivals);
            std::vector<double> dtdxs(nArrivals);
            std::vector<double> dtdys(nArrivals);
            std::vector<double> dtdzs(nArrivals);
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               &dtdt0s, &dtdxs, &dtdys, &dtdzs,
                                               mApplyCorrection);
            // In the following we leverage implicit partial derivatives:
            //  d f(x, y, z = E(x,y))/dx = df/dx + df/dE dE/dx
            //  d f(x, y, z = E(x,y))/dy = df/dy + df/dE dE/dy
            if (mNorm == ::Norm::LeastSquares)
            {
                constexpr double two{2};
                for (int i = 0; i < nArrivals; ++i)
                {
                    double residual = mObservations[i] - arrivalTimes[i];
                    double twoWeightSquared = two*(mWeights[i]*mWeights[i]);
                    double twoWeightSquaredResidual = twoWeightSquared*residual;
                    gradient[0] = gradient[0]
                                - twoWeightSquaredResidual*dtdt0s[i];
                    gradient[1] = gradient[1]
                                - twoWeightSquaredResidual
                                 *(dtdxs[i] + dtdzs[i]*dEdx);
                    gradient[2] = gradient[2]
                                - twoWeightSquaredResidual
                                 *(dtdys[i] + dtdzs[i]*dEdy);
                }
            }
            else if (mNorm == Norm::L1)
            {
                for (int i = 0; i < nArrivals; ++i)
                {
                    double weightedResidual
                        = mWeights[i]*(mObservations[i] - arrivalTimes[i]);
                    auto sgnWeightedResidual = ::sgn(weightedResidual);
                    double weightedSgnWeightedResidual
                        = mWeights[i]*sgnWeightedResidual;
                    gradient[0] = gradient[0]
                                - weightedSgnWeightedResidual*dtdt0s[i];
                    gradient[1] = gradient[1]
                                - weightedSgnWeightedResidual
                                 *(dtdxs[i] + dtdzs[i]*dEdx);
                    gradient[2] = gradient[2]
                                - weightedSgnWeightedResidual
                                 *(dtdys[i] + dtdzs[i]*dEdy);
                }
            }
            else if (mNorm == Norm::Lp)
            {
                double costFunctionInverse
                    = ::lp(mWeights, mObservations, arrivalTimes,
                           mPNorm, Measurement::Standard);
                if (costFunctionInverse > 0)
                {
                    costFunctionInverse = 1./costFunctionInverse;
                } 
                for (int i = 0; i < nArrivals; ++i)
                {
                    double weightedResidual
                        = mWeights[i]*(mObservations[i] - arrivalTimes[i]);
                    double weightedAbsoluteResidual
                        = std::abs(weightedResidual);
                    double argument
                        = weightedAbsoluteResidual*costFunctionInverse; 
                    double scalar
                        =-(mWeights[i]*std::pow(argument, mPNorm - 1))
                         *::sgn(weightedResidual);
                    gradient[0] = gradient[0] + scalar*dtdt0s[i];
                    gradient[1] = gradient[1]
                                + scalar*(dtdxs[i] + dtdzs[i]*dEdx);
                    gradient[2] = gradient[2]
                                + scalar*(dtdys[i] + dtdzs[i]*dEdy);
                }
            }
#ifndef NDEBUG
            else
            {
                assert(false);
            }
#endif
        }
        else
        {
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               mApplyCorrection);
        }
        if (mNorm == Norm::LeastSquares)
        {
            costFunction
                = ::leastSquares(mWeights, mObservations, arrivalTimes,
                                 Measurement::Standard);
        }
        else if (mNorm == Norm::L1)
        {
            costFunction
                = ::l1(mWeights, mObservations, arrivalTimes,
                       Measurement::Standard);
        }
        else if (mNorm == Norm::Lp)
        {
            costFunction
                = ::lp(mWeights, mObservations, arrivalTimes,
                       mPNorm, Measurement::Standard);
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        double logLikelihood =-costFunction;
        if (doGradient)
        {
            std::transform(gradient, gradient + n, gradient,
                           [](const auto &gi)
                           {
                               return -gi;
                           });
        } 
        return logLikelihood;
    }
    /// @result The hard bounds on the search region.  Most optimization tools
    ///         can handle linear constraints.  Otherwise, I could add a penalty 
    ///         to the objective function that requires the origin time be
    ///         less than the first arrival time.
    std::pair<std::vector<double>, std::vector<double>> getLinearBounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    const ULocator::TravelTimeCalculatorMap *mTravelTimeCalculatorMap{nullptr};
    const ULocator::Topography::ITopography *mTopography{nullptr};
    std::vector<std::pair<ULocator::Station, std::string>> mStationPhases;
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    Norm mNorm{::Norm::LeastSquares};
    double mPNorm{1.5};
    double mPenaltyCoefficient{1};
    bool mApplyCorrection{true};
};

}
#endif
