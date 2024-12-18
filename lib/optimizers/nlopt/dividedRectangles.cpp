#include <string>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include <nlopt.hpp>
#include "uLocator/optimizers/nlopt/dividedRectangles.hpp" 
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "objectiveFunction.hpp"

using namespace ULocator::Optimizers::NLOpt;

class DividedRectangles::DividedRectanglesImpl
{
public:
    DividedRectanglesImpl(
        std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    { 
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    double mMaximumDepth{65000}; // Maximum event depth
    double mInitialDepth{6000}; // Initial depth guess
    double mTimeWindow{120}; // The time window to search
    double mLocationTolerance{1}; // 1 meter is close enough
    double mOriginTimeTolerance{0.01}; // This is fine 
    double mOptimalValue{std::numeric_limits<double>::max()};
    int mMaximumObjectiveFunctionEvaluations{5000};
    int mNumberOfObjectiveFunctionEvaluations{0};
    int mNumberOfGradientEvaluations{0};
    bool mLocallyBias{true};
    bool mHaveOrigin{false};
    bool mNormalize{true};
};

/// Constructor
DividedRectangles::DividedRectangles() :
    ULocator::Optimizers::IOptimizer(),
    pImpl(std::make_unique<DividedRectanglesImpl> ()) 
{
}

/// Constructor
DividedRectangles::DividedRectangles(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    ULocator::Optimizers::IOptimizer(logger),
    pImpl(std::make_unique<DividedRectanglesImpl> (logger))
{
}

/// Maximum function evaluations
void DividedRectangles::setMaximumNumberOfObjectiveFunctionEvaluations(
    const int nEvaluations)
{
    if (nEvaluations < 1)
    {
        throw std::invalid_argument("Number of evaluations must be positive");
    }
    pImpl->mMaximumObjectiveFunctionEvaluations = nEvaluations;
}

int DividedRectangles::getMaximumNumberOfObjectiveFunctionEvaluations() const noexcept
{
    return pImpl->mMaximumObjectiveFunctionEvaluations;
}

/// Model tolerance
void DividedRectangles::setOriginTimeTolerance(const double tolerance)
{
    if (tolerance < 0)
    {
        throw std::invalid_argument("Origin time tolerance must be positive");
    }
    pImpl->mOriginTimeTolerance = tolerance;
}

double DividedRectangles::getOriginTimeTolerance() const noexcept
{
    return pImpl->mOriginTimeTolerance;
}

void DividedRectangles::setLocationTolerance(const double tolerance)
{
    if (tolerance < 0)
    {
        throw std::invalid_argument("Location tolerance must be positive");
    }
    pImpl->mLocationTolerance = tolerance;
}

double DividedRectangles::getLocationTolerance() const noexcept
{
    return pImpl->mLocationTolerance;
}

/// Locally bias?
void DividedRectangles::enableLocallyBias() noexcept
{
    pImpl->mLocallyBias = true;
}

void DividedRectangles::disableLocallyBias() noexcept
{
    pImpl->mLocallyBias = false;
}

bool DividedRectangles::locallyBias() const noexcept
{
    return pImpl->mLocallyBias;
}

/// Normalize?
void DividedRectangles::enableNormalization() noexcept
{
    pImpl->mNormalize = true;
}

void DividedRectangles::disableNormalization() noexcept
{
    pImpl->mNormalize = false;
}

bool DividedRectangles::normalize() const noexcept
{
    return pImpl->mNormalize;
}

void DividedRectangles::setOriginTimeSearchWindowDuration(
    const double duration)
{
    if (duration <= 0)
    {
        throw std::invalid_argument("Duration must be positive");
    }
    pImpl->mTimeWindow = duration;
}

double DividedRectangles::getOriginTimeSearchWindowDuration() const noexcept
{
    return pImpl->mTimeWindow;
}

/// Locates 
void DividedRectangles::locate(
    const ULocator::Origin &initialGuess,
    const IOptimizer::LocationProblem locationProblem,
    const IOptimizer::Norm norm)
{
    constexpr bool reduceTimes{true};
    pImpl->mHaveOrigin = false;
    pImpl->mOptimalValue = std::numeric_limits<double>::max();
    pImpl->mNumberOfObjectiveFunctionEvaluations = 0;
    pImpl->mNumberOfGradientEvaluations = 0;
    // Throws
    auto region = ULocator::Optimizers::IOptimizer::getGeographicRegion();
    const auto &arrivals
        = ULocator::Optimizers::IOptimizer::getArrivalsReference();
    if (arrivals.empty())
    {
        throw std::runtime_error("No arrivals");
    }
   
    // Figure out the bounds
    double z0 =-ULocator::Optimizers::IOptimizer::getTopography()->
                getMinimumAndMaximumElevation().first;
    double z1 = pImpl->mMaximumDepth;
    double t0 =-getOriginTimeSearchWindowDuration(); //pImpl->mTimeWindow;
    double t1 = 0; // Reduced arrival time - first arrival is t = 0
    auto [x0, x1] = region->getExtentInX();
    auto [y0, y1] = region->getExtentInY();

    // Figure out the directAlgorithm
    std::vector<double> estimateArrivalTimes;
    auto directAlgorithm = nlopt::GN_DIRECT_L;
    if (normalize())
    {
        if (locallyBias())
        {
            pImpl->mLogger->debug("Using DIRECT locally biased algorithm");
            directAlgorithm = nlopt::GN_DIRECT_L;
        }
        else
        {
            pImpl->mLogger->debug("Using DIRECT classic algorithm");
            directAlgorithm = nlopt::GN_DIRECT; 
        }
    }
    else
    {
        if (locallyBias())
        {
            pImpl->mLogger->debug(
                "Using DIRECT locally biased unscaled algorithm");
            directAlgorithm = nlopt::GN_DIRECT_L_NOSCAL;
        }
        else
        {
            pImpl->mLogger->debug("Using DIRECT unscaled algorithm");
            directAlgorithm = nlopt::GN_DIRECT_NOSCAL;
        }
    }

    // Will use origin information after convergence
    Origin origin;
    if (locationProblem == IOptimizer::LocationProblem::ThreeDimensionsAndTime)
    {
        pImpl->mLogger->debug("Initializing DIRECT 3D and time problem...");
        ::NLOptProblem3DAndTime optimizer(directAlgorithm);
        optimizer.mNorm = norm;
        optimizer.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        optimizer.mTopography
            = ULocator::Optimizers::IOptimizer::getTopography();
        // Boundaries
        optimizer.setSearchBoundaries(std::vector<double> {t0, x0, y0, z0},
                                      std::vector<double> {t1, x1, y1, z1} );
        // Convergence criteria
        optimizer.setMaximumNumberOfObjectiveFunctionEvaluations(
            getMaximumNumberOfObjectiveFunctionEvaluations());
        optimizer.setLocationAndTimeTolerance(getLocationTolerance(),
                                              getOriginTimeTolerance());
        // Copy observations/weights/etc.
        optimizer.fillObservationsWeightsStationPhasesArraysFromArrivals(
             arrivals, reduceTimes);
        // Initial guess
        auto xLocation = optimizer.createInitialGuess(initialGuess, *region);
        // Optimize
        double optimalValue{std::numeric_limits<double>::max()};
        try
        {
            pImpl->mLogger->debug("Optimizing...");
            optimizer.mOptimizer.optimize(xLocation, optimalValue);
            pImpl->mOptimalValue = optimalValue;
            pImpl->mNumberOfObjectiveFunctionEvaluations
                = optimizer.mObjectiveFunctionEvaluations;
            pImpl->mNumberOfGradientEvaluations
                = optimizer.mGradientEvaluations;
        }
        catch (const std::exception &e)
        {
            auto errorMessage = "3D and time DIRECT failed with: "
                              + std::string {e.what()};
            pImpl->mLogger->error(errorMessage);
            throw std::runtime_error(errorMessage);
        }
        // Extract the origin information
        origin = optimizer.locationToOrigin(xLocation, *region);
        // Estimate the travel times
        pImpl->mLogger->debug("Computing estimate travel times...");
        optimizer.mTravelTimeCalculatorMap->evaluate(
            optimizer.mStationPhases,
            origin.getTime(),
            xLocation.at(1),
            xLocation.at(2),
            origin.getDepth(),
            &estimateArrivalTimes,
            optimizer.mApplyCorrection);
    }
    else if (locationProblem ==
             IOptimizer::LocationProblem::FixedToFreeSurfaceAndTime)
    {
        pImpl->mLogger->debug("Initializing DIRECT free surface problem...");
        ::NLOptProblem2DAndTimeAndDepthAtFreeSurface optimizer(directAlgorithm);
        optimizer.mNorm = norm;
        optimizer.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        optimizer.mTopography
            = ULocator::Optimizers::IOptimizer::getTopography();
        // Boundaries
        optimizer.setSearchBoundaries(std::vector<double> {t0, x0, y0},
                                      std::vector<double> {t1, x1, y1} );
        // Convergence criteria
        optimizer.setMaximumNumberOfObjectiveFunctionEvaluations(
            getMaximumNumberOfObjectiveFunctionEvaluations());
        optimizer.setLocationAndTimeTolerance(getLocationTolerance(),
                                              getOriginTimeTolerance());
        // Copy observations/weights/etc.
        optimizer.fillObservationsWeightsStationPhasesArraysFromArrivals(
             arrivals, reduceTimes);
        // Initial guess
        auto xLocation = optimizer.createInitialGuess(initialGuess, *region);
        pImpl->mLogger->debug("Initial location: (t,x,y) = ("
                             + std::to_string(optimizer.mReductionTime
                                            + xLocation.at(0)) + "," 
                             + std::to_string(xLocation.at(1)) + "," 
                             + std::to_string(xLocation.at(2)) + ")");
        // Optimize
        double optimalValue{std::numeric_limits<double>::max()};
        try
        {
            pImpl->mLogger->debug("Optimizing...");
            optimizer.mOptimizer.optimize(xLocation, optimalValue);
            pImpl->mOptimalValue = optimalValue;
            pImpl->mNumberOfObjectiveFunctionEvaluations
                = optimizer.mObjectiveFunctionEvaluations;
            pImpl->mNumberOfGradientEvaluations
                = optimizer.mGradientEvaluations;
        }
        catch (const std::exception &e) 
        {
            auto errorMessage = "Fixed to surface and time DIRECT failed with: "
                              + std::string {e.what()};
            pImpl->mLogger->error(errorMessage);
            throw std::runtime_error(errorMessage);
        }
        // Extract the origin information
        origin = optimizer.locationToOrigin(xLocation, *region);
        // Estimate the travel times
        pImpl->mLogger->debug("Computing estimate travel times...");
        optimizer.mTravelTimeCalculatorMap->evaluate(
            optimizer.mStationPhases,
            origin.getTime(),
            xLocation.at(1),
            xLocation.at(2),
            origin.getDepth(),
            &estimateArrivalTimes,
            optimizer.mApplyCorrection);
    }
    else if (locationProblem ==
             IOptimizer::LocationProblem::FixedDepthAndTime)
    {
        pImpl->mLogger->debug("Initializing DIRECT fixed-depth problem...");
        ::NLOptProblem2DAndTimeAndFixedDepth optimizer(directAlgorithm);
        optimizer.mNorm = norm;
        optimizer.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        // Boundaries
        optimizer.setSearchBoundaries(std::vector<double> {t0, x0, y0},
                                      std::vector<double> {t1, x1, y1} );
        // Convergence criteria
        optimizer.setMaximumNumberOfObjectiveFunctionEvaluations(
            getMaximumNumberOfObjectiveFunctionEvaluations());
        optimizer.setLocationAndTimeTolerance(getLocationTolerance(),
                                              getOriginTimeTolerance());
        // Copy observations/weights/etc.
        optimizer.fillObservationsWeightsStationPhasesArraysFromArrivals(
             arrivals, reduceTimes);
        // Initial guess
        auto xLocation = optimizer.createInitialGuess(initialGuess, *region);
        double sourceDepth = std::min(std::max(z0, pImpl->mInitialDepth), z1);
        if (initialGuess.haveDepth())
        {
            sourceDepth = initialGuess.getDepth();
        }
        optimizer.mDepth = sourceDepth;
        pImpl->mLogger->debug("Initial location: (t,x,y,z) = ("
                             + std::to_string(optimizer.mReductionTime
                                            + xLocation.at(0)) + ","
                             + std::to_string(xLocation.at(1)) + ","
                             + std::to_string(xLocation.at(2)) + ","
                             + std::to_string(sourceDepth) + ")");
        // Copy observations/weights/etc.
        //optimizer.fillObservationsWeightsStationPhasesArraysFromArrivals(
        //     arrivals, reduceTimes);
        // Optimize
        double optimalValue{std::numeric_limits<double>::max()};
        try
        {
            pImpl->mLogger->debug("Optimizing...");
            optimizer.mOptimizer.optimize(xLocation, optimalValue);
            pImpl->mOptimalValue = optimalValue;
            pImpl->mNumberOfObjectiveFunctionEvaluations
                = optimizer.mObjectiveFunctionEvaluations;
            pImpl->mNumberOfGradientEvaluations
                = optimizer.mGradientEvaluations;
        }
        catch (const std::exception &e) 
        {
            auto errorMessage = "Fixed depth and time DIRECT failed with: "
                              + std::string {e.what()};
            pImpl->mLogger->error(errorMessage);
            throw std::runtime_error(errorMessage);
        }
        // Extract the origin information
        origin = optimizer.locationToOrigin(xLocation, *region);
        // Estimate the travel times
        pImpl->mLogger->debug("Computing estimate travel times...");
        optimizer.mTravelTimeCalculatorMap->evaluate(
            optimizer.mStationPhases,
            origin.getTime(),
            xLocation.at(1),
            xLocation.at(2),
            origin.getDepth(),
            &estimateArrivalTimes,
            optimizer.mApplyCorrection);
    }
    else
    {
#ifndef NDEBUG
        assert(false);
#else
        throw std::runtime_error("Unhandled location problem");
#endif
    }
    // Insert the arrival times
    pImpl->mLogger->debug("Computing residuals...");
    auto newArrivals = arrivals;
    for (int i = 0; i < static_cast<int> (newArrivals.size()); ++i)
    {
        newArrivals[i].setResidual(arrivals[i].getTime()
                                 - estimateArrivalTimes.at(i));
    }   
    origin.setArrivals(std::move(newArrivals));
    // Sets the origin
    pImpl->mLogger->debug("Setting the origin...");
    IOptimizer::setOrigin(origin);
    pImpl->mHaveOrigin = true;
}

/// Destructor
DividedRectangles::~DividedRectangles() = default;

/// Have origin?
bool DividedRectangles::haveOrigin() const noexcept
{
    return pImpl->mHaveOrigin;
}

/// Objective function evaluations
int DividedRectangles::getNumberOfObjectiveFunctionEvaluations() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mNumberOfObjectiveFunctionEvaluations;
}

/// Gradient evaluations
int DividedRectangles::getNumberOfGradientEvaluations() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mNumberOfGradientEvaluations;
}

/// Optimal objective function value
double DividedRectangles::getOptimalObjectiveFunction() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mOptimalValue;
}
