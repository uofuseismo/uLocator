#ifndef ULOCATOR_OPTIMIZERS_PRIMA_OBJECTIVE_FUNCTION_HPP
#define ULOCATOR_OPTIMIZERS_PRIMA_OBJECTIVE_FUNCTION_HPP
#include <iostream>
#include <functional>
#include <vector>
#include <array>
#include <prima/prima.h>
#ifndef NDEBUG
#include <cassert>
#endif
#include "optimizers/objectiveFunctions.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/station.hpp"
namespace
{

void getReturnCodeAndThrow(const int returnCode)
{
    if (returnCode == PRIMA_SMALL_TR_RADIUS){return;} 
    if (returnCode == PRIMA_FTARGET_ACHIEVED){return;}
    if (returnCode == PRIMA_INVALID_INPUT)
    {
        throw std::runtime_error("Internal error - Prima invalid input");
    }
    if (returnCode == PRIMA_NULL_OPTIONS)
    {
        throw std::runtime_error("Internal error - Prima NULL options");
    }
    if (returnCode == PRIMA_NULL_X0)
    {
        throw std::runtime_error("Internal error - Prima NULL problem");
    }
    if (returnCode == PRIMA_NULL_RESULT)
    {
        throw std::runtime_error("Internal error - Prima NULL result");
    }
    if (returnCode == PRIMA_NULL_FUNCTION)
    {
        throw std::runtime_error("Internal error - Prima null function");
    } 
std::cerr << "Unhandled prima error: " << static_cast<int> (returnCode) << std::endl;
/*
typedef enum
{
  PRIMA_SMALL_TR_RADIUS = 0,
  PRIMA_FTARGET_ACHIEVED = 1,
  PRIMA_TRSUBP_FAILED = 2,
  PRIMA_MAXFUN_REACHED = 3,
  PRIMA_MAXTR_REACHED = 20, 
  PRIMA_NAN_INF_X = -1, 
  PRIMA_NAN_INF_F = -2, 
  PRIMA_NAN_INF_MODEL = -3, 
  PRIMA_NO_SPACE_BETWEEN_BOUNDS = 6,
  PRIMA_DAMAGING_ROUNDING = 7,
  PRIMA_ZERO_LINEAR_CONSTRAINT = 8,
  PRIMA_CALLBACK_TERMINATE = 30, 
  PRIMA_INVALID_INPUT = 100,
  PRIMA_ASSERTION_FAILS = 101,
  PRIMA_VALIDATION_FAILS = 102,
  PRIMA_MEMORY_ALLOCATION_FAILS = 103,
  PRIMA_NULL_OPTIONS = 110,
  PRIMA_NULL_PROBLEM = 111,
  PRIMA_NULL_X0 = 112,
  PRIMA_NULL_RESULT = 113,
  PRIMA_NULL_FUNCTION = 114,
*/
}

static void primaProblem3DAndTimeCallback(const double x[],
                                          double *f, 
                                          const void *data);

struct PrimaProblem3DAndTime : public ::Problem3DAndTime
{
    typedef void(PrimaProblem3DAndTime::*Callback)(const double [], double *, const void *);
    explicit PrimaProblem3DAndTime(const prima_algorithm_t algorithm) :
        mAlgorithm(algorithm)
    {
        // Create the problem
        prima_init_problem(&mPrimaProblem, 4);
        mPrimaProblem.x0 = mX0.data();
        mPrimaProblem.calfun = &primaProblem3DAndTimeCallback;
        // And the base options
        prima_init_options(&mPrimaOptions);
        mPrimaOptions.rhobeg = 100;   // Reasonable starting search condition
        mPrimaOptions.rhoend = 1.e-5; // Reasonable terminating search condition
        mPrimaOptions.maxfun = 1500;
        mPrimaOptions.data = this;
    }
    ~PrimaProblem3DAndTime()
    {
        //prima_free_problem(&mPrimaProblem);
        //prima_free_options(&mPrimaOptions);
    }
    /// @brief Applies model for objective function and, optionally, gradient.
    double evaluate(const double x[]) const
    {
        double objectiveFunction{0};
        objectiveFunction
            =-::Problem3DAndTime::logLikelihood(mParameters, x);
        if (mSaveHistory)
        {
            for (int i = 0; i < mParameters; ++i)
            {
                mModelHistory.push_back(x[i]);
            }
            mObjectiveFunctionHistory.push_back(objectiveFunction);
        }
        return objectiveFunction;
    }
    /// @result The location vector converted to an origin.
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const std::vector<double> &xLocation,
                     const ULocator::Position::IGeographicRegion &region)
    {
        ULocator::Origin origin; 
        // Extract the origin information
        origin.setTime(mReductionTime + xLocation.at(0));
        double xSource{xLocation.at(1)};
        double ySource{xLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region.localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(xLocation.at(3));
        return origin;
    }
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const prima_result_t &primaResult,
                     const ULocator::Position::IGeographicRegion &region)
    {
        std::vector<double> x(mParameters);
        std::copy(primaResult.x, primaResult.x + x.size(), x.begin());
        return locationToOrigin(x, region);
    }
    /// @brief Sets the search boundaries.
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        { 
            throw std::invalid_argument("lowerBoundaries size != 4");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 4");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
        mPrimaProblem.xl = mLowerBounds.data();
        mPrimaProblem.xu = mUpperBounds.data(); 
    }
    /// @brief Sets the maximum number of objective function eavluations.
    void setMaximumNumberOfObjectiveFunctionEvaluations(const int nEvaluations)
    {
        if (nEvaluations < mParameters)
        {
            throw std::invalid_argument("nEvaluations must be at least 4");
        }
        mPrimaOptions.maxfun = nEvaluations;
    }
    void setLocationAndTimeTolerance(const double xTolerance,
                                     const double timeTolerance)
    {
        mPrimaOptions.rhoend = std::min(xTolerance, timeTolerance);
    }
    void setInitialGuess(
        const ULocator::Origin &initialGuess,
        const ULocator::Position::IGeographicRegion &region)
    {
        if (mLowerBounds.empty())
        {
            throw std::runtime_error("Search boundaries not set");
        }
        // Initial guess
        std::vector<double> xStartingLocation(mParameters);
        for (int i = 0; i < mParameters; ++i)
        {
            mX0.at(i) = 0.5*(mLowerBounds.at(i)
                           + mUpperBounds.at(i));
        }
        if (initialGuess.haveTime())
        {
            auto time = initialGuess.getTime() - mReductionTime;
            if (time >= mLowerBounds.at(0) &&
                time <= mUpperBounds.at(0))
            {
                mX0.at(0) = time;
            }
        }
        if (initialGuess.haveEpicenter())
        {
            auto epicenter = initialGuess.getEpicenter(); 
            auto [xGuess, yGuess] = region.geographicToLocalCoordinates(
                epicenter.getLatitude(), epicenter.getLongitude());
            if (xGuess >= mLowerBounds.at(1) &&
                xGuess <= mUpperBounds.at(1))
            {
                mX0.at(1) = xGuess;
            }
            if (yGuess >= mLowerBounds.at(2) &&
                yGuess <= mUpperBounds.at(2))
            {
                mX0.at(2) = yGuess;
            }
        }
        if (initialGuess.haveDepth())
        {
            auto depth = initialGuess.getDepth();
            if (depth >= mLowerBounds.at(3) &&
                depth <= mUpperBounds.at(3))
            {
                mX0.at(3) = depth;
            }
        }
    }
    prima_problem_t mPrimaProblem;
    prima_options_t mPrimaOptions;
    std::array<double, 4> mX0{0, 0, 0, 0};
    mutable std::vector<double> mModelHistory;
    mutable std::vector<double> mObjectiveFunctionHistory;
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    int mParameters{4};
    prima_algorithm_t mAlgorithm;
    bool mSaveHistory{false};
};

void primaProblem3DAndTimeCallback(const double x[],
                                   double *f,
                                   const void *data)
{
#ifndef NDEBUG
     assert(x != nullptr);
     assert(f != nullptr);
     assert(data != nullptr);
#endif
     auto objectiveFunction
         = reinterpret_cast<const PrimaProblem3DAndTime *> (data);
     *f = objectiveFunction->evaluate(x); 
}

/*
///--------------------------------------------------------------------------///

struct PrimaProblem2DAndTimeAndFixedDepth :
     public ::Problem2DAndTimeAndFixedDepth
{
    explicit PrimaProblem2DAndTimeAndFixedDepth(
        const nlopt::algorithm algorithm) :
        mOptimizer(algorithm, 3)
    {   
        mOptimizer.set_min_objective(mObjectiveFunction);
    }
    double operator()(unsigned int n, const double *x, double *gradient)
    {
#ifndef NDEBUG
        assert(n == mParameters);
        assert(x != nullptr);
#endif
        double objectiveFunction{0};
        if (gradient == nullptr)
        {
            mObjectiveFunctionEvaluations = mObjectiveFunctionEvaluations + 1;
            objectiveFunction
                =-::Problem2DAndTimeAndFixedDepth::logLikelihood(n, x);
        }
        else
        {
            mGradientEvaluations = mGradientEvaluations + 1;
            objectiveFunction
                =-::Problem2DAndTimeAndFixedDepth
                  ::logLikelihoodAndGradient(n, x, gradient);
            std::transform(gradient, gradient + n, gradient,
                           [](const double gi) 
                           {
                               return -gi;
                           });
        }
        if (mSaveHistory)
        {
            for (int i = 0; i < mParameters; ++i)
            {
                mModelHistory.push_back(x[i]);
            }
            mObjectiveFunctionHistory.push_back(objectiveFunction);
        }
        return objectiveFunction;
    }
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const std::vector<double> &xLocation,
                     const ULocator::Position::IGeographicRegion &region)
    {   
        ULocator::Origin origin; 
        // Extract the origin information
        origin.setTime(mReductionTime + xLocation.at(0));
        double xSource{xLocation.at(1)};
        double ySource{xLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region.localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(mDepth);
        return origin;
    }
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {   
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 3");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 3");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mOptimizer.set_lower_bounds(lowerBoundaries);
        mOptimizer.set_upper_bounds(upperBoundaries);
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }
    void setMaximumNumberOfObjectiveFunctionEvaluations(const int nEvaluations)
    {
        if (nEvaluations < 1)
        {
            throw std::invalid_argument("nEvaluations must be positive");
        }
        mOptimizer.set_maxeval(nEvaluations);
    }
    void setLocationAndTimeTolerance(const double xTolerance,
                                     const double timeTolerance)
    {
        std::vector<double> tolerance(mParameters, 0);
        tolerance.at(0) = timeTolerance;
        tolerance.at(1) = xTolerance;
        tolerance.at(2) = xTolerance;
        mOptimizer.set_xtol_abs(tolerance);
    }
    void resetCounters()
    {
        mObjectiveFunctionEvaluations = 0;
        mGradientEvaluations = 0;
    }
    [[nodiscard]] std::vector<double> createInitialGuess(
        const ULocator::Origin &initialGuess,
        const ULocator::Position::IGeographicRegion &region)
    {
        if (mLowerBounds.empty())
        {
            throw std::runtime_error("Search boundaries not set");
        }
        // Initial guess
        std::vector<double> xStartingLocation(mParameters, 0);
        for (int i = 0; i < mParameters; ++i)
        {
            xStartingLocation.at(i) = 0.5*(mLowerBounds.at(i)
                                         + mUpperBounds.at(i));
        }
        if (initialGuess.haveTime())
        {
            auto time = initialGuess.getTime() - mReductionTime;
            if (time >= mLowerBounds.at(0) &&
                time <= mUpperBounds.at(0))
            {
                xStartingLocation.at(0) = time;
            }
        }
        if (initialGuess.haveEpicenter())
        {
            auto epicenter = initialGuess.getEpicenter(); 
            auto [xGuess, yGuess] = region.geographicToLocalCoordinates(
                epicenter.getLatitude(), epicenter.getLongitude());
            if (xGuess >= mLowerBounds.at(1) &&
                xGuess <= mUpperBounds.at(1))
            {
                xStartingLocation.at(1) = xGuess;
            }
            if (yGuess >= mLowerBounds.at(2) &&
                yGuess <= mUpperBounds.at(2))
            {
                xStartingLocation.at(2) = yGuess;
            }
        }
        return xStartingLocation;
    }
    std::function<double (unsigned int, const double *, double *) >
        mObjectiveFunction
        {
            std::bind(&PrimaProblem2DAndTimeAndFixedDepth::operator(),
                      this,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3)
        };
    nlopt::opt mOptimizer;
    std::vector<double> mModelHistory;
    std::vector<double> mObjectiveFunctionHistory;
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    int mObjectiveFunctionEvaluations{0};
    int mGradientEvaluations{0};
    int mParameters{3};
    bool mSaveHistory{false};
};

///--------------------------------------------------------------------------///

struct PrimaProblem2DAndTimeAndDepthAtFreeSurface :
     public ::Problem2DAndTimeAndDepthAtFreeSurface
{
    explicit PrimaProblem2DAndTimeAndDepthAtFreeSurface(
        const nlopt::algorithm algorithm) :
        mOptimizer(algorithm, 3)
    {
        mOptimizer.set_min_objective(mObjectiveFunction);
    }
    double operator()(unsigned int n, const double *x, double *gradient)
    {   
#ifndef NDEBUG
        assert(n == mParameters);
        assert(x != nullptr);
#endif
        double objectiveFunction{0};
        if (gradient == nullptr)
        {
            mObjectiveFunctionEvaluations = mObjectiveFunctionEvaluations + 1;
            objectiveFunction
                =-::Problem2DAndTimeAndDepthAtFreeSurface::logLikelihood(n, x);
        }
        else
        {
            mGradientEvaluations = mGradientEvaluations + 1;
            objectiveFunction
                =-::Problem2DAndTimeAndDepthAtFreeSurface
                  ::logLikelihoodAndGradient(n, x, gradient);
            std::transform(gradient, gradient + n, gradient,
                           [](const double gi) 
                           {
                               return -gi;
                           });
        }
        if (mSaveHistory)
        {
            for (int i = 0; i < mParameters; ++i)
            {
                mModelHistory.push_back(x[i]);
            }
            mObjectiveFunctionHistory.push_back(objectiveFunction);
        }
        return objectiveFunction;
    }
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const std::vector<double> &xLocation,
                     const ULocator::Position::IGeographicRegion &region)
    {   
        ULocator::Origin origin; 
        // Extract the origin information
        origin.setTime(mReductionTime + xLocation.at(0));
        double xSource{xLocation.at(1)};
        double ySource{xLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region.localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        double sourceDepth =-mTopography->evaluate(xSource, ySource);
        origin.setDepth(sourceDepth);
        return origin;
    }
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 3");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 3");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mOptimizer.set_lower_bounds(lowerBoundaries);
        mOptimizer.set_upper_bounds(upperBoundaries);
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }   
    void setMaximumNumberOfObjectiveFunctionEvaluations(const int nEvaluations)
    {
        if (nEvaluations < 1)
        {
            throw std::invalid_argument("nEvaluations must be positive");
        }
        mOptimizer.set_maxeval(nEvaluations);
    }
    void setLocationAndTimeTolerance(const double xTolerance,
                                     const double timeTolerance)
    {
        std::vector<double> tolerance(mParameters, 0);
        tolerance.at(0) = timeTolerance;
        tolerance.at(1) = xTolerance;
        tolerance.at(2) = xTolerance;
        mOptimizer.set_xtol_abs(tolerance);
    }
    void resetCounters()
    {
        mObjectiveFunctionEvaluations = 0;
        mGradientEvaluations = 0;
    }
    [[nodiscard]] std::vector<double> createInitialGuess(
        const ULocator::Origin &initialGuess,
        const ULocator::Position::IGeographicRegion &region)
    {
        if (mLowerBounds.empty())
        {
            throw std::runtime_error("Search boundaries not set");
        }
        // Initial guess
        std::vector<double> xStartingLocation(mParameters, 0);
        for (int i = 0; i < mParameters; ++i)
        {
            xStartingLocation.at(i) = 0.5*(mLowerBounds.at(i)
                                         + mUpperBounds.at(i));
        }
        if (initialGuess.haveTime())
        {
            auto time = initialGuess.getTime() - mReductionTime;
            if (time >= mLowerBounds.at(0) &&
                time <= mUpperBounds.at(0))
            {
                xStartingLocation.at(0) = time;
            }
        }
        if (initialGuess.haveEpicenter())
        {
            auto epicenter = initialGuess.getEpicenter();
            auto [xGuess, yGuess] = region.geographicToLocalCoordinates(
                epicenter.getLatitude(), epicenter.getLongitude());
            if (xGuess >= mLowerBounds.at(1) &&
                xGuess <= mUpperBounds.at(1))
            {
                xStartingLocation.at(1) = xGuess;
            }
            if (yGuess >= mLowerBounds.at(2) &&
                yGuess <= mUpperBounds.at(2))
            {
                xStartingLocation.at(2) = yGuess;
            }
        }
        return xStartingLocation;
    }
    std::function<double (unsigned int, const double *, double *) >
        mObjectiveFunction
        {
            std::bind(&PrimaProblem2DAndTimeAndDepthAtFreeSurface::operator(),
                      this,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3)
        };
    nlopt::opt mOptimizer;
    std::vector<double> mModelHistory;
    std::vector<double> mObjectiveFunctionHistory;
    std::vector<double> mLowerBounds;
    std::vector<double> mUpperBounds;
    int mObjectiveFunctionEvaluations{0};
    int mGradientEvaluations{0};
    int mParameters{3};
    bool mSaveHistory{false};
};
*/

}
#endif
