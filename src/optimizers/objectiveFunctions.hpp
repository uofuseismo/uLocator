#ifndef ULOCATOR_OPTIMIZERS_OBJECTIVE_FUNCTIONS_HPP
#define ULOCATOR_OPTIMIZERS_OBJECTIVE_FUNCTIONS_HPP
#include <iostream>
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
    L2
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
               const std::vector<T> &estimates)
{
     const T *__restrict__ weightsPtr = weights.data();
     const T *__restrict__ observationsPtr = observations.data();
     const T *__restrict__ estimatesPtr = estimates.data();
     auto n = static_cast<int> (weights.size());
     double sumSquaredOfResiduals{0};
     for (int i = 0; i < n; ++i)
     {
         double weightedResidual
            = weightsPtr[i]*(observationsPtr[i] - estimatesPtr[i]);
         sumSquaredOfResiduals = sumSquaredOfResiduals
                               + weightedResidual*weightedResidual;
     }
     return static_cast<T> (sumSquaredOfResiduals);
}

template<typename T>
[[nodiscard]]
T l1(const std::vector<T> &weights,
     const std::vector<T> &observations,
     const std::vector<T> &estimates)
{
     const T *__restrict__ weightsPtr = weights.data();
     const T *__restrict__ observationsPtr = observations.data();
     const T *__restrict__ estimatesPtr = estimates.data();
     double sumAbsoluteResiduals{0};
     for (size_t i = 0; i < observations.size(); ++i)
     {
         double weightedResidual
            = weightsPtr[i]*std::abs(observationsPtr[i] - estimatesPtr[i]);
         sumAbsoluteResiduals = sumAbsoluteResiduals + weightedResidual;
     }
     return static_cast<T> (sumAbsoluteResiduals);
}

//----------------------------------------------------------------------------//

struct LeastSquaresProblem3DAndTime
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
            std::vector<double> dtdt0s(mStationPhases.size());
            std::vector<double> dtdxs(mStationPhases.size());
            std::vector<double> dtdys(mStationPhases.size());
            std::vector<double> dtdzs(mStationPhases.size());
            mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                               originTime,
                                               xSource, ySource, zSource,
                                               &arrivalTimes,
                                               &dtdt0s, &dtdxs, &dtdys, &dtdzs,
                                               mApplyCorrection);
            if (mNorm == Norm::L2)
            {
                constexpr double two{2};
                for (int i = 0; i < static_cast<int> (arrivalTimes.size()); ++i)
                {
                    double residual = mObservations[i] - arrivalTimes[i];
                    double twoWeightSquared = two*(mWeights[i]*mWeights[i]);
                    double twoWeightSquaredResidual = twoWeightSquared*residual;
                    gradient[0] = gradient[0] - twoWeightSquaredResidual*dtdt0s[i];
                    gradient[1] = gradient[1] - twoWeightSquaredResidual*dtdxs[i];
                    gradient[2] = gradient[2] - twoWeightSquaredResidual*dtdys[i];
                    gradient[3] = gradient[3] - twoWeightSquaredResidual*dtdzs[i];
                }
            }
            else if (mNorm == Norm::L1)
            {
                for (int i = 0; i < static_cast<int> (arrivalTimes.size()); ++i)
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
        if (mNorm == Norm::L2)
        {
            costFunction
                 = ::leastSquares(mWeights, mObservations, arrivalTimes);
        }
        else if (mNorm == Norm::L1)
        {
            costFunction
                 = ::l1(mWeights, mObservations, arrivalTimes);
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
    Norm mNorm{Norm::L2};
    double mPenaltyCoefficient{1};
    bool mApplyCorrection{true};
};

//----------------------------------------------------------------------------//

/*
struct LeastSquaresProblem2DAndTimeAfixToFreeSurface
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {   
#ifndef NDEBUG
        assert(dv.size() == 3); 
#endif
        double originTime = dv[0];
        double xSource = dv[1];
        double ySource = dv[2];
        auto [latitude, longitude]
            = mRegion->localToGeographicCoordinates(xSource, ySource);
        auto elevation = mTopography->evaluate(latitude, longitude);
        double zSource =-elevation; // Elevation increases +up; make it +down
        std::vector<double> arrivalTimes(mStationPhases.size());
        mTravelTimeCalculatorMap->evaluate(mStationPhases,
                                           originTime,
                                           xSource, ySource, zSource,
                                           &arrivalTimes,
                                           mApplyCorrection);
        double fitnessScore{0};
        if (mNorm == Norm::L2)
        {
            fitnessScore
                 = ::leastSquares(mWeights, mObservations, arrivalTimes);
        }
        else if (mNorm == Norm::L1)
        {
            fitnessScore
                 = ::l1(mWeights, mObservations, arrivalTimes);
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        return pagmo::vector_double {fitnessScore};
    }
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    const ULocator::TravelTimeCalculatorMap *mTravelTimeCalculatorMap{nullptr};
    std::shared_ptr<ULocator::Topography::ITopography> mTopography{nullptr};
    std::unique_ptr<ULocator::Position::IGeographicRegion> mRegion{nullptr};
    std::vector<std::pair<ULocator::Station, std::string>> mStationPhases;
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
    Norm mNorm{Norm::L2};
    bool mApplyCorrection{true};
};

//----------------------------------------------------------------------------//

struct LeastSquaresProblem2DAndTimeFixedDepth
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {   
#ifndef NDEBUG
        assert(dv.size() == 3); 
#endif
        double originTime = dv[0];
        double xSource = dv[1];
        double ySource = dv[2];
        std::vector<double> arrivalTimes(mStationPhases.size());
        mTravelTimeCalculatorMap->evaluate(mStationPhases,
                       originTime, xSource, ySource, mSourceDepth,
                       &arrivalTimes,
                       mApplyCorrection);
        double fitnessScore{0};
        if (mNorm == Norm::L2)
        {
            fitnessScore
                 = ::leastSquares(mWeights, mObservations, arrivalTimes);
        }
        else if (mNorm == Norm::L1)
        {
            fitnessScore
                 = ::l1(mWeights, mObservations, arrivalTimes);
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        return pagmo::vector_double {fitnessScore};
    }   
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {   
        return std::pair {mLowerBounds, mUpperBounds};
    }   
    const ULocator::TravelTimeCalculatorMap *mTravelTimeCalculatorMap{nullptr};
    std::vector<std::pair<ULocator::Station, std::string>> mStationPhases;
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
    double mSourceDepth{0};
    Norm mNorm{Norm::L2};
    bool mApplyCorrection{true};
};
*/
}
#endif
