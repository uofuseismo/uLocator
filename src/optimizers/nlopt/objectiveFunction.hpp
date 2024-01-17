#ifndef ULOCATOR_OPTIMIZERS_NLOPT_OBJECTIVE_FUNCTION_HPP
#define ULOCATOR_OPTIMIZERS_NLOPT_OBJECTIVE_FUNCTION_HPP
#include <iostream>
#include <functional>
#include <vector>
#include <nlopt.hpp>
#ifndef NDEBUG
#include <cassert>
#endif
#include "../objectiveFunctions.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/station.hpp"
namespace
{

struct NLOptProblem3DAndTime : public ::Problem3DAndTime
{
    explicit NLOptProblem3DAndTime(const nlopt::algorithm algorithm) :
        mOptimizer(algorithm, 4)
    {
        mOptimizer.set_min_objective(mObjectiveFunction);
    }
    /// @brief Applies model for objective function and, optionally, gradient.
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
                =-::Problem3DAndTime::logLikelihood(n, x);
        }
        else
        {
            mGradientEvaluations = mGradientEvaluations + 1;
            objectiveFunction
                =-::Problem3DAndTime::logLikelihoodAndGradient(n, x, gradient);
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
    /// @brief Resets objective function counters.
    void resetCounters()
    {
        mObjectiveFunctionEvaluations = 0;
        mGradientEvaluations = 0;
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
        mOptimizer.set_lower_bounds(lowerBoundaries);
        mOptimizer.set_upper_bounds(upperBoundaries);
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }
    /// @brief Sets the maximum number of objective function eavluations.
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
        tolerance.at(3) = xTolerance;
        mOptimizer.set_xtol_abs(tolerance);
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
        std::vector<double> xStartingLocation(mParameters);
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
        if (initialGuess.haveDepth())
        {
            auto depth = initialGuess.getDepth();
            if (depth >= mLowerBounds.at(3) &&
                depth <= mUpperBounds.at(3))
            {
                xStartingLocation.at(3) = depth;
            }
        }
        return xStartingLocation;
    }
    std::function<double (unsigned int, const double *, double *) >
        mObjectiveFunction
        {
            std::bind(&NLOptProblem3DAndTime::operator(),
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
    int mParameters{4};
    bool mSaveHistory{false};
};

///--------------------------------------------------------------------------///

struct NLOptProblem2DAndTimeAndFixedDepth :
     public ::Problem2DAndTimeAndFixedDepth
{
    explicit NLOptProblem2DAndTimeAndFixedDepth(
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
            std::bind(&NLOptProblem2DAndTimeAndFixedDepth::operator(),
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

struct NLOptProblem2DAndTimeAndDepthAtFreeSurface :
     public ::Problem2DAndTimeAndDepthAtFreeSurface
{
    explicit NLOptProblem2DAndTimeAndDepthAtFreeSurface(
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
            std::bind(&NLOptProblem2DAndTimeAndDepthAtFreeSurface::operator(),
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

}
#endif
