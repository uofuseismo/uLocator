#include <vector>
#include <string>
#include <limits>
#include <prima/prima.h>
#include <umps/logging/standardOut.hpp>
#include "uLocator/optimizers/prima/boundOptimizationByQuadraticApproximation.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "objectiveFunction.hpp"

using namespace ULocator::Optimizers::Prima;

class BoundOptimizationByQuadraticApproximation::
      BoundOptimizationByQuadraticApproximationImpl
{
public:
    BoundOptimizationByQuadraticApproximationImpl(
        std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {   
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> (); 
        }
    }   
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::pair<double, double> mExtentInX;
    std::pair<double, double> mExtentInY;
    double mOptimalValue{std::numeric_limits<double>::max()};
    double mMaximumDepth{65000};
    double mTimeWindow{250};// // Default time window to search
    double mLocationTolerance{1};
    double mOriginTimeTolerance{1.e-3};
    int mNumberOfObjectiveFunctionEvaluations{0};
    int mMaximumObjectiveFunctionEvaluations{1500};
    bool mHaveOrigin{false};
    bool mHaveExtentInX{false};
    bool mHaveExtentInY{false};
};

/// Constructor
BoundOptimizationByQuadraticApproximation::
    BoundOptimizationByQuadraticApproximation() :
    ULocator::Optimizers::IOptimizer(),
    pImpl(std::make_unique<BoundOptimizationByQuadraticApproximationImpl> ()) 
{
}

/// Constructor
BoundOptimizationByQuadraticApproximation::BoundOptimizationByQuadraticApproximation(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    ULocator::Optimizers::IOptimizer(logger),
    pImpl(std::make_unique<BoundOptimizationByQuadraticApproximationImpl> (logger))
{
}

/// Destructor
BoundOptimizationByQuadraticApproximation::
    ~BoundOptimizationByQuadraticApproximation() = default;

/// Maximum function evaluations
void BoundOptimizationByQuadraticApproximation::
    setMaximumNumberOfObjectiveFunctionEvaluations(const int nEvaluations)
{
    if (nEvaluations < 1)
    {   
        throw std::invalid_argument("Number of evaluations must be positive");
    }   
    pImpl->mMaximumObjectiveFunctionEvaluations = nEvaluations;
}

int BoundOptimizationByQuadraticApproximation::
    getMaximumNumberOfObjectiveFunctionEvaluations() const noexcept
{
    return pImpl->mMaximumObjectiveFunctionEvaluations;
}


void BoundOptimizationByQuadraticApproximation::setExtentInX(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {   
        throw std::invalid_argument("extent.second <= extent.first in x");
    }   
    pImpl->mExtentInX = extent;
    pImpl->mHaveExtentInX = true;
}

/// Model tolerance
void BoundOptimizationByQuadraticApproximation::setOriginTimeTolerance(
    const double tolerance)
{
    if (tolerance < 0)
    {   
        throw std::invalid_argument("Origin time tolerance must be positive");
    }   
    pImpl->mOriginTimeTolerance = tolerance;
}

double BoundOptimizationByQuadraticApproximation::getOriginTimeTolerance() 
    const noexcept
{
    return pImpl->mOriginTimeTolerance;
}

void  BoundOptimizationByQuadraticApproximation::setLocationTolerance(
    const double tolerance)
{
    if (tolerance < 0)
    {   
        throw std::invalid_argument("Location tolerance must be positive");
    }   
    pImpl->mLocationTolerance = tolerance;
}

double BoundOptimizationByQuadraticApproximation::
    getLocationTolerance() const noexcept
{
    return pImpl->mLocationTolerance;
}

void BoundOptimizationByQuadraticApproximation::
    setOriginTimeSearchWindowDuration(const double duration)
{
    if (duration <= 0)
    {
        throw std::invalid_argument("Duration must be positive");
    }
    pImpl->mTimeWindow = duration;
}

double BoundOptimizationByQuadraticApproximation::
    getOriginTimeSearchWindowDuration() const noexcept
{
    return pImpl->mTimeWindow;
}

/// Extent
void BoundOptimizationByQuadraticApproximation::setExtentInY(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {   
        throw std::invalid_argument("extent.second <= extent.first in y");
    }   
    pImpl->mExtentInY = extent;
    pImpl->mHaveExtentInY = true;
}

std::pair<double, double> BoundOptimizationByQuadraticApproximation::getExtentInX() const
{
    if (!pImpl->mHaveExtentInX)
    {   
        throw std::runtime_error("Extent in x not set");
    }   
    return pImpl->mExtentInX;
}

std::pair<double, double> BoundOptimizationByQuadraticApproximation::getExtentInY() const
{
    if (!pImpl->mHaveExtentInY)
    {   
        throw std::runtime_error("Extent in y not set");
    }   
    return pImpl->mExtentInY;
}

/// Have origin?
bool BoundOptimizationByQuadraticApproximation::haveOrigin() const noexcept
{
    return pImpl->mHaveOrigin;
}

/// Objective function evaluations
int
BoundOptimizationByQuadraticApproximation::
    getNumberOfObjectiveFunctionEvaluations() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mNumberOfObjectiveFunctionEvaluations;
}

/// Gradient evaluations
int BoundOptimizationByQuadraticApproximation::
    getNumberOfGradientEvaluations() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return 0; // Derivative free
}

/// Optimal objective function value
double BoundOptimizationByQuadraticApproximation::
    getOptimalObjectiveFunction() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mOptimalValue;
}


/// Locates 
void BoundOptimizationByQuadraticApproximation::locate(
    const ULocator::Origin &initialGuess,
    const IOptimizer::LocationProblem locationProblem,
    const IOptimizer::Norm norm)
{
    constexpr bool reduceTimes{true};
    pImpl->mHaveOrigin = false;
    pImpl->mOptimalValue = std::numeric_limits<double>::max();
    pImpl->mNumberOfObjectiveFunctionEvaluations = 0;
    // Throws
    auto region = ULocator::Optimizers::IOptimizer::getGeographicRegion();
    const auto &arrivals
        = ULocator::Optimizers::IOptimizer::getArrivalsReference();
    if (arrivals.empty())
    {
        throw std::runtime_error("No arrivals");
    }
    // Ensure extent is set
    if (!pImpl->mHaveExtentInX)
    {
        setExtentInX(region->getExtentInX());
    }
    if (!pImpl->mHaveExtentInY)
    {
        setExtentInY(region->getExtentInY());
    }
    // Extract and reduce observations
    std::vector<std::pair<ULocator::Station, std::string>>
        stationPhases(arrivals.size());
    std::vector<double> weights(arrivals.size());
    std::vector<double> observations(arrivals.size());
    std::vector<double> estimateArrivalTimes(arrivals.size(), 0);
    for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
    {
        stationPhases[i]
            = std::pair{arrivals[i].getStation(), arrivals[i].getPhase()};
        weights[i] = 1./arrivals[i].getStandardError();
        observations[i] = arrivals[i].getTime();
    }
    // Figure out the bounds
    double z0 =-ULocator::Optimizers::IOptimizer::getTopography()->
                getMinimumAndMaximumElevation().first;
    double z1 = pImpl->mMaximumDepth;
    double t0 =-getOriginTimeSearchWindowDuration();
    double t1 = 0; // Reduced arrival time - first arrival is t = 0
    auto [x0, x1] = getExtentInX();
    auto [y0, y1] = getExtentInY();

    // Will use origin information after convergence
    Origin origin;
    if (locationProblem == IOptimizer::LocationProblem::ThreeDimensionsAndTime)
    {   
        pImpl->mLogger->debug("Initializing BOBYQA 3D and time problem...");
        ::PrimaProblem3DAndTime optimizer(PRIMA_BOBYQA);
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
        optimizer.setInitialGuess(initialGuess, *region);
        // Optimize
        double optimalValue{std::numeric_limits<double>::max()};
        prima_result_t primaResult;
        try
        {
            pImpl->mLogger->debug("Optimizing...");
            prima_minimize(optimizer.mAlgorithm, 
                           &optimizer.mPrimaProblem,
                           &optimizer.mPrimaOptions,
                           &primaResult);
            ::getReturnCodeAndThrow(primaResult.status);
pImpl->mLogger->info(std::string {primaResult.message});
            //pImpl->mOptimalValue = optimalValue;
            pImpl->mNumberOfObjectiveFunctionEvaluations = primaResult.nf;
            origin = optimizer.locationToOrigin(primaResult, *region);
            pImpl->mOptimalValue = primaResult.f;
        }
        catch (const std::exception &e) 
        {
            auto errorMessage = "3D and time BOBYQA failed with: "
                              + std::string {e.what()};
            //pImpl->mLogger->error(errorMessage);
            prima_free_result(&primaResult);
            throw std::runtime_error(errorMessage);
        }
        prima_free_result(&primaResult);
        // Estimate the travel times
        pImpl->mLogger->debug("Computing estimate travel times...");
        auto [xLocation, yLocation]
            = region->geographicToLocalCoordinates(
                 origin.getEpicenter().getLatitude(),
                 origin.getEpicenter().getLongitude() );
        optimizer.mTravelTimeCalculatorMap->evaluate(
            optimizer.mStationPhases,
            origin.getTime(),
            xLocation,
            yLocation,
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
                                 - estimateArrivalTimes[i]);
    }   
    origin.setArrivals(std::move(newArrivals));
    // Sets the origin
    pImpl->mLogger->debug("Setting the origin...");
    IOptimizer::setOrigin(origin);
    pImpl->mHaveOrigin = true;

}


