#include <iostream>
#include <iomanip>
#include <string>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include <nlopt.hpp>
#include "uLocator/nlopt.hpp"
#include "uLocator/nloptOptions.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/objectiveFunctions/leastSquares.hpp"
#include "uLocator/objectiveFunctions/l1.hpp"
#include "uLocator/objectiveFunctions/lp.hpp"
#include "uLocator/objectiveFunctions/doubleDifferenceL1.hpp"
#include "uLocator/quarry.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/station.hpp"
#include "uLocator/topography.hpp"
#include "uLocator/position/wgs84.hpp"

using namespace ULocator;

namespace
{
struct TopographyConstraintData
{
    std::shared_ptr<UMPS::Logging::ILog> logger{nullptr};
    std::function<double (const double latitude, const double longitude)
                 > topographyCallbackFunction;
    double defaultElevation{0}; 
    int utmZone{12};
    bool north{true};
};

double topographyConstraintFunction(
    const std::vector<double> &x, std::vector<double> &gradient, void *inputData)
{   
    auto data = reinterpret_cast<::TopographyConstraintData *> (inputData); 
    Position::WGS84 position{data->utmZone, data->north, x.at(0), x.at(1)};
    auto latitude  = position.getLatitude();
    auto longitude = position.getLongitude();
    auto depth = x.at(2);
    double topography = data->defaultElevation;
    try
    {
        topography = data->topographyCallbackFunction(latitude, longitude);
    }
    catch (const std::exception &e) 
    {
        data->logger->warn(std::string{e.what()} + "; using default elevation");
    }
    // NLOpt example appears to be in form of: x_3 >= f(x_2, x_1) where,
    // for our purposes, f is a function of the topography.
    auto topographicDepth =-topography; // Negate for consistent coord system
    // Convert it to 0 >= f(x_1, x_2) - x_3 hence the feasable region is now
    // below the topography
    return topographicDepth - depth;
}

}

class NLOpt::NLOptImpl
{
public:
    NLOptImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    std::pair<double, double> toUTM(const double latitude, const double longitude)
    {
        Position::WGS84 position{latitude, longitude, mUTMZone}; 
        return std::pair {position.getEasting(), position.getNorthing()};
    }
    std::pair<double, double> toLatLon(const double utmX, const double utmY)
    {
        Position::WGS84 position{mUTMZone, mNorth, utmX, utmY};
        return std::pair {position.getLatitude(), position.getLongitude()};
    }
    void switchData()
    {
    }
    void setDirectOptions( )
    {
    }
    [[nodiscard]] double fixedDepthCallback() const
    {
        return mFixedDepth;
    }
    [[nodiscard]] double defaultTopographyCallback(const double, const double) const
    {
        return mDefaultElevation;
    }
    [[nodiscard]] double topographyCallback(const double latitude,
                                            const double longitude) const
    {
        if (mTopography != nullptr)
        {
            if (mTopography->haveTopography())
            {
                return mTopography->evaluate(latitude, longitude);
            }
        }
        return defaultTopographyCallback(latitude, longitude);
    }
    std::function<double () > mFixedDepthCallbackFunction
    {
        std::bind(&NLOptImpl::fixedDepthCallback, this)
    };
    std::function<double (const double latitude, const double longitude)
                 > mTopographyCallbackFunction
    {
        std::bind(&NLOptImpl::topographyCallback,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2)
    }; 
    std::vector<double> guessInitialEpicenter()
    {
        std::vector<double> x(2);
        auto arrivals = mObjectiveFunction->getArrivals();
        auto index = std::min_element(arrivals.begin(), arrivals.end(),
                                      [=](const Arrival &lhs, const Arrival &rhs)
                                      {
                                           return lhs.getTime() < rhs.getTime();
                                      });
#ifndef NDEBUG
        assert(index != arrivals.end());
#endif
        auto stationReference = index->getStationReference();
        auto stationPosition = stationReference.getGeographicPosition();
        auto latitudeBoundaries  = mOptions.getLatitudeBoundaries();
        auto longitudeBoundaries = mOptions.getLongitudeBoundaries();
        x[0] = std::min(std::max(stationPosition.getLatitude(),
                                     latitudeBoundaries[0] + 0.01),
                                     latitudeBoundaries[1] - 0.01);
        x[1] = std::min(std::max(stationPosition.getLongitude(),
                                 longitudeBoundaries[0] + 0.01),
                                 longitudeBoundaries[1] - 0.01);
         
        return x;
    }
    std::vector<double> searchQuarries()
    {    
        auto xInitial = guessInitialEpicenter();
        auto depth = fixedDepthCallback(); 
        Origin stationOrigin;
        stationOrigin.setEpicenter(
            Position::WGS84 {xInitial[0], xInitial[1], mUTMZone} ); 
        stationOrigin.setDepth(fixedDepthCallback());
        double maximumObjectiveFunction
             = mObjectiveFunction->evaluateLoss(stationOrigin);
        auto latitudeBoundaries  = mOptions.getLatitudeBoundaries();
        auto longitudeBoundaries = mOptions.getLongitudeBoundaries();
        for (const auto &quarry : mQuarries)
        {
            std::vector<double> x(2);
            auto candidateEpicenter = quarry.getGeographicPosition();
            auto latitude  = candidateEpicenter.getLatitude();
            auto longitude = candidateEpicenter.getLongitude();
            if (latitude < latitudeBoundaries[0] ||
                latitude > latitudeBoundaries[1])
            {
                continue;
            }
            if (longitude < longitudeBoundaries[0] ||
                longitude > longitudeBoundaries[1])
            {
                continue;
            }
            Origin quarryOrigin;
            quarryOrigin.setEpicenter(candidateEpicenter);
            double quarryDepth
                =-mTopographyCallbackFunction(latitude, longitude);
            quarryOrigin.setDepth(quarryDepth);
            auto objectiveFunction
                = mObjectiveFunction->evaluateLoss(quarryOrigin);
            if (objectiveFunction < maximumObjectiveFunction)
            {
                //std::cout << quarry.getName() << std::endl;
                maximumObjectiveFunction = objectiveFunction;
                xInitial[0] = latitude;
                xInitial[1] = longitude;
            }
        }
        return xInitial;
    }
    void locate(const SourceDepthConstraint sourceConstraint,
                const Origin &initialGuess,
                const bool isBlast)
    {
        mObjectiveFunctionCounter = 0;
        mGradientCounter = 0;
        if (!mHaveObjectiveFunction)
        {
            throw std::runtime_error("Objective function (options) not set");
        }   
        if (!mObjectiveFunction->haveTravelTimeCalculatorMap())
        {
            throw std::runtime_error("Travel time calculators not set");
        }
        mHaveOrigin = false;
#ifndef NDEBUG
        assert(mObjectiveFunction != nullptr);
#endif
/*
        ObjectiveFunctions::IObjectiveFunction *objectiveFunction{nullptr};
        if (mOptions.getObjectiveFunction() ==
            NLOptOptions::ObjectiveFunction::LeastSquares)
        {
            objectiveFunction = mLeastSquaresObjectiveFunction.get();
        }
        else if (mOptions.getObjectiveFunction() ==
                 NLOptOptions::ObjectiveFunction::L1)
        {
            objectiveFunction = mL1ObjectiveFunction.get();
        }
        else if (mOptions.getObjectiveFunction() ==
                 NLOptOptions::ObjectiveFunction::LeastSquaresDoubleDifference)
        {
            //objectiveFunction = mLeastSquaresDoubleDifferenceObjectiveFunction.get();
        }
        else
        {
            throw std::runtime_error("Unhandled objective function");
        }
*/
        if (mObjectiveFunction->getNumberOfArrivals() < 1)
        {
            throw std::runtime_error("Arrivals not set");
        }
        // Set the depth (won't matter for fixed to free surface)
        mFixedDepth = mOptions.getInitialEarthquakeSearchDepth();
        if (isBlast){mFixedDepth = mOptions.getInitialQuarryBlastSearchDepth();}
        if (initialGuess.haveDepth())
        {
            mFixedDepth = initialGuess.getDepth();
        }
        // Set the inversion strategy
        if (sourceConstraint == SourceDepthConstraint::FixedToFreeSurface ||
            sourceConstraint == SourceDepthConstraint::Fixed)
        {
            if (mObjectiveFunction->getNumberOfArrivals() < 3)
            {
                mLogger->warn(
                "Attempting to fixed depth location with less than 3 arrivals");
            }
            if (sourceConstraint == SourceDepthConstraint::FixedToFreeSurface)
            {
                mObjectiveFunction->setInversionStrategy(
                    InversionStrategy::LatitudeLongitudeDepthAtFreeSurface);
                mObjectiveFunction->setTopographyCallbackFunction(
                    mTopographyCallbackFunction);
            }
            else
            {
                mObjectiveFunction->setInversionStrategy(
                   InversionStrategy::LatitudeLongitudeFixedDepth);
            }
        }
        else
        {
            if (mObjectiveFunction->getNumberOfArrivals() < 4)
            {
                mLogger->warn("Attempting to locate with less than 4 arrivals");
            }
            mObjectiveFunction->setInversionStrategy(
                InversionStrategy::LatitudeLongitudeFixedDepth);
        }
        // Set up the initial optimization scheme
        mNorth = mOptions.isNorthernHemisphere();
        // Set up the optimization scheme
        mObjectiveFunction->toggleApplyCorrection(true); 
        mObjectiveFunction->resetCounters();
        auto nParameters2D = mObjectiveFunction->getNumberOfModelParameters();
#ifndef NDEBUG
        assert(nParameters2D == 2);
#endif
        auto nlOptObjectiveFunction
             = std::bind(&ObjectiveFunctions::IObjectiveFunction::operator(),
                         mObjectiveFunction,
                         std::placeholders::_1,
                         std::placeholders::_2,
                         std::placeholders::_3);
        // Get the boundaries
        auto latitudeBoundaries  = mOptions.getLatitudeBoundaries();
        auto longitudeBoundaries = mOptions.getLongitudeBoundaries();
        auto depthBoundaries     = mOptions.getDepthBoundaries();
        auto [minimumUTMX, minimumUTMY] = toUTM(latitudeBoundaries[0],
                                                longitudeBoundaries[0]);
        auto [maximumUTMX, maximumUTMY] = toUTM(latitudeBoundaries[1],
                                                longitudeBoundaries[1]);
        // Initial optimization
        nlopt::opt initialOptimizer(nlopt::GN_DIRECT_L, nParameters2D);
        initialOptimizer.set_min_objective(nlOptObjectiveFunction);
        initialOptimizer.set_maxeval(mOptions.getInitialMaximumNumberOfFunctionEvaluations());
        initialOptimizer.set_lower_bounds(
            std::vector<double> {minimumUTMX, minimumUTMY});
        initialOptimizer.set_upper_bounds(
            std::vector<double> {maximumUTMX, maximumUTMY});
        initialOptimizer.set_xtol_abs(
            mOptions.getInitialAbsoluteModelTolerance());
        std::vector<double> xInitialLatitudeAndLongitude{
             0.5*(latitudeBoundaries[0] + latitudeBoundaries[1]),
             0.5*(longitudeBoundaries[0] + longitudeBoundaries[1])
        };
        bool useHeuristic{false};
        if (initialGuess.haveEpicenter())
        {
            auto latitude  = initialGuess.getEpicenter().getLatitude();
            auto longitude = initialGuess.getEpicenter().getLongitude();
            if (latitude < latitudeBoundaries[0] || 
                latitude > latitudeBoundaries[1])
            {
                mLogger->warn("Initial latitude not in search boundaries.  Using heuristics...");
                useHeuristic = true;
            }
            if (longitude < longitudeBoundaries[0] || 
                longitude > longitudeBoundaries[1])
            {
                mLogger->warn("Initial longitude not in search boundaries.  Using heuristics...");
                useHeuristic = true;
            }
            if (!useHeuristic)
            {
                xInitialLatitudeAndLongitude[0] = latitude;
                xInitialLatitudeAndLongitude[1] = longitude;
            }
        }
        else
        {
            useHeuristic = true;
        }
        if (useHeuristic)
        {
            if (!isBlast)
            {
                xInitialLatitudeAndLongitude = guessInitialEpicenter();
            }
            else
            {
                xInitialLatitudeAndLongitude = searchQuarries();
            }
        }
        mLogger->debug("Using initial (lat,lon) = ("
                     + std::to_string( xInitialLatitudeAndLongitude[0] )
                     + ","
                     + std::to_string( xInitialLatitudeAndLongitude[1] )
                     + ")");
        auto [initialSourceUTMX, initialSourceUTMY]
            = toUTM(xInitialLatitudeAndLongitude[0],
                    xInitialLatitudeAndLongitude[1]);
        std::vector<double> xInitialLocation(nParameters2D);
        xInitialLocation.at(0) = initialSourceUTMX;
        xInitialLocation.at(1) = initialSourceUTMY;
        // Finally perform initial location
        double initialOptimalValue{0};
        try
        {
            mLogger->debug("Performing initial location...");
            initialOptimizer.optimize(xInitialLocation, initialOptimalValue);
        }
        catch (const std::exception &e)
        {
            mLogger->error(e.what());
            throw std::runtime_error("Initial optimization failed");
        }
        auto [initialLatitude, initialLongitude]
            = toLatLon(xInitialLocation[0], xInitialLocation[1]);
        Origin initialOrigin;
        initialOrigin.setEpicenter(
            Position::WGS84 {initialLatitude, initialLongitude, mUTMZone});
        initialOrigin.setDepth(fixedDepthCallback());
        initialOrigin = mObjectiveFunction->predict(initialOrigin);
        double initialDepth = fixedDepthCallback();
        if (sourceConstraint == SourceDepthConstraint::FixedToFreeSurface)
        {
            initialDepth =-topographyCallback(initialLatitude,
                                              initialLongitude);
        } 
        mLogger->debug("Initial DIRECT location (lat,lon,depth,time,value)=("
             +       std::to_string(initialOrigin.getEpicenter().getLatitude())
             + "," + std::to_string(initialOrigin.getEpicenter().getLongitude())
             + "," + std::to_string(initialOrigin.getDepth())
             + "," + std::to_string(initialOrigin.getTime())
             + "," + std::to_string(initialOptimalValue)
             + ")");
        //------------------------ Now locate for real -----------------------//
        mObjectiveFunctionCounter
             = mObjectiveFunction->getNumberOfObjectiveFunctionEvaluations();
        mGradientCounter
             = mObjectiveFunction->getNumberOfGradientEvaluations();
        mObjectiveFunction->resetCounters();
        // Output product
        double finalLatitude{initialOrigin.getEpicenter().getLatitude()};
        double finalLongitude{initialOrigin.getEpicenter().getLongitude()};
        double finalDepth{initialOrigin.getDepth()};
        double finalOptimalValue{initialOptimalValue};
        // Refine the search region
        double latitudeRefinement  = mOptions.getLatitudeRefinement();
        double longitudeRefinement = mOptions.getLongitudeRefinement();
        double refinedMinimumUTMX = xInitialLocation.at(0) - latitudeRefinement;
        double refinedMaximumUTMX = xInitialLocation.at(0) + latitudeRefinement;
        double refinedMinimumUTMY = xInitialLocation.at(1) - longitudeRefinement;
        double refinedMaximumUTMY = xInitialLocation.at(1) + longitudeRefinement;
        if (sourceConstraint == SourceDepthConstraint::FixedToFreeSurface ||
            sourceConstraint == SourceDepthConstraint::Fixed)
        {
            nlopt::opt optimizer(nlopt::LD_MMA, nParameters2D);
            optimizer.set_min_objective(nlOptObjectiveFunction);
            optimizer.set_lower_bounds(
                std::vector<double> {refinedMinimumUTMX, refinedMinimumUTMY});
            optimizer.set_upper_bounds(
                std::vector<double> {refinedMaximumUTMX, refinedMaximumUTMY});
            optimizer.set_maxeval(mOptions.getMaximumNumberOfFunctionEvaluations());
            optimizer.set_xtol_abs(mOptions.getAbsoluteModelTolerance());
            std::vector<double> xLocation(nParameters2D);
            xLocation.at(0) = xInitialLocation.at(0);
            xLocation.at(1) = xInitialLocation.at(1);
            double optimalValue{0};
            optimizer.optimize(xLocation, optimalValue);
            if (optimalValue <= initialOptimalValue)
            {
                auto [latitude, longitude]
                   = toLatLon(xLocation[0], xLocation[1]);
                double depth = fixedDepthCallback();
                if (sourceConstraint == SourceDepthConstraint::FixedToFreeSurface)
                {
                    depth =-topographyCallback(latitude, longitude);
                }
                mLogger->debug("Using MMA 2D location with objective function:"
                             + std::to_string(optimalValue));
                finalLatitude   = latitude;
                finalLongitude  = longitude;
                finalDepth      = depth; 
                finalOptimalValue = optimalValue;
            }
            else
            {
                mLogger->warn(
                   "MMA did not reduce obj function; (initial,refined)=("
                 + std::to_string(initialOptimalValue) 
                 + ","
                 + std::to_string(optimalValue)
                 + "); using initial location");
            }
        }
        else if (sourceConstraint == SourceDepthConstraint::Free)
        {
            mObjectiveFunction->setInversionStrategy(
                InversionStrategy::LatitudeLongitudeDepth);
            auto nParameters = mObjectiveFunction->getNumberOfModelParameters();
            nlopt::opt optimizer(nlopt::LD_MMA, nParameters);
            optimizer.set_min_objective(nlOptObjectiveFunction);
            optimizer.set_maxeval(mOptions.getMaximumNumberOfFunctionEvaluations());
            optimizer.set_xtol_abs(mOptions.getAbsoluteModelTolerance());
            // Topographic constraint
            ::TopographyConstraintData
                topographyConstraintData{mLogger,
                                         mTopographyCallbackFunction,
                                         mDefaultElevation, mUTMZone, mNorth};
            optimizer.set_lower_bounds(
                std::vector<double> {refinedMinimumUTMX,
                                     refinedMinimumUTMY,
                                     depthBoundaries[0]});
            optimizer.set_upper_bounds(
                std::vector<double> {refinedMaximumUTMX,
                                     refinedMaximumUTMY,
                                     depthBoundaries[1]});
            optimizer.add_inequality_constraint(topographyConstraintFunction,
                                                &topographyConstraintData,
                                                0.1);
            // Loop through trial depths
            double bestOptimalValue = std::numeric_limits<double>::max();
            std::vector<double> finalLocation;
            std::vector<double> candidateDepths 
                = mOptions.getRefinementSearchDepths();
            if (candidateDepths.empty())
            {
                candidateDepths.push_back(fixedDepthCallback());
            }
            for (const auto &candidateDepth : candidateDepths)
            {
                std::vector<double> xLocation(nParameters);
                xLocation.at(0) = xInitialLocation.at(0);
                xLocation.at(1) = xInitialLocation.at(1);
                xLocation.at(2) = candidateDepth;
                double optimalValue{0};
                optimizer.optimize(xLocation, optimalValue);
                if (optimalValue < bestOptimalValue)
                {
                    finalLocation = xLocation;
                    bestOptimalValue = optimalValue;
                } 
            }
            if (bestOptimalValue < initialOptimalValue)
            {
                mLogger->debug("Using MMA best objective function: "
                             + std::to_string(bestOptimalValue));
            }
            else
            {
                mLogger->warn(
                   "MMA did not reduce obj function; (initial,refined)=("
                 + std::to_string(initialOptimalValue) 
                 + "," 
                 + std::to_string(bestOptimalValue)
                 + "); using initial location");
                bestOptimalValue = initialOptimalValue;
                finalLocation.resize(3);
                finalLocation.at(0) = xInitialLocation.at(0);
                finalLocation.at(1) = xInitialLocation.at(1);
                finalLocation.at(2) = fixedDepthCallback();
            }
            auto [latitude, longitude]
                = toLatLon(finalLocation.at(0), finalLocation.at(1));
            finalLatitude   = latitude;
            finalLongitude  = longitude;
            finalDepth      = finalLocation.at(2);
            finalOptimalValue = bestOptimalValue;
        }
        else
        {
            throw std::runtime_error("Unhandled source depth constraint");
        }
        // Make final predictions
        Origin origin;
        origin.setEpicenter(
            Position::WGS84 {finalLatitude, finalLongitude, mUTMZone});
        origin.setDepth(finalDepth);
        // origin.setTime(finalOriginTime); DO NOT SET THIS!
        mOrigin = mObjectiveFunction->predict(origin); // Reoptimize for origin
        mHaveOrigin = true;
        // Write something for the users
        mObjectiveFunctionCounter
             = mObjectiveFunction->getNumberOfObjectiveFunctionEvaluations();
        mGradientCounter
             = mObjectiveFunction->getNumberOfGradientEvaluations();
        mLogger->debug("Location finished with (lat,lon,depth,time,value)=("
                   +       std::to_string(mOrigin.getEpicenter().getLatitude())
                   + "," + std::to_string(mOrigin.getEpicenter().getLongitude())
                   + "," + std::to_string(mOrigin.getDepth())
                   + "," + std::to_string(mOrigin.getTime())
                   + "," + std::to_string(finalOptimalValue)
                   + ")"); 
    }

    NLOptOptions mOptions;
    std::unique_ptr<const Topography> mTopography{nullptr};
    std::unique_ptr<const TravelTimeCalculatorMap> mTravelTimeCalculatorMap{nullptr};
    std::unique_ptr<ObjectiveFunctions::LeastSquares> mLeastSquaresObjectiveFunction{nullptr};
    std::unique_ptr<ObjectiveFunctions::L1> mL1ObjectiveFunction{nullptr};
    std::unique_ptr<ObjectiveFunctions::LP> mLPObjectiveFunction{nullptr};
    std::unique_ptr<ObjectiveFunctions::DoubleDifferenceL1> mDoubleDifferenceL1ObjectiveFunction{nullptr};
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::vector<Quarry> mQuarries;
    ObjectiveFunctions::IObjectiveFunction *mObjectiveFunction{nullptr};
    Origin mOrigin;
    double mDefaultElevation{1500};
    double mFixedDepth{0};
    int mUTMZone{12};
    int mGradientCounter{0};
    int mObjectiveFunctionCounter{0};
    bool mNorth{true};
    bool mHaveObjectiveFunction{false};
    bool mInitialized{false};
    bool mHaveOrigin{false};
};

/// Constructor
NLOpt::NLOpt() : 
    pImpl(std::make_unique<NLOptImpl> ())
{
}

/// Constructor
NLOpt::NLOpt(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<NLOptImpl> (logger))
{
}

/// Destructor
NLOpt::~NLOpt() = default;

/// Reset class
void NLOpt::clear() noexcept
{
    auto logger = pImpl->mLogger;
    pImpl = std::make_unique<NLOptImpl> (logger);
}

/// Initialized?
bool NLOpt::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Set options
void NLOpt::setOptions(const NLOptOptions &options)
{
    clear();
    // Unset the travel time calculator
    if (pImpl->mObjectiveFunction)
    {
        if (pImpl->mObjectiveFunction
                 ->haveTravelTimeCalculatorMap())
        {
            pImpl->mTravelTimeCalculatorMap
                = pImpl->mObjectiveFunction->releaseTravelTimeCalculatorMap();
        }
    }
    pImpl->mObjectiveFunction = nullptr;
    pImpl->mHaveObjectiveFunction = false;
    if (options.getObjectiveFunction() ==
        NLOptOptions::ObjectiveFunction::LeastSquares)
    {
        pImpl->mLogger->debug("Will use least-squares objective function...");
        pImpl->mLeastSquaresObjectiveFunction
            = std::make_unique<ULocator::ObjectiveFunctions::LeastSquares> (pImpl->mLogger);
        pImpl->mLeastSquaresObjectiveFunction->setUTMZone(
            options.getUTMZone(),
            options.isNorthernHemisphere());
        pImpl->mLeastSquaresObjectiveFunction->setFixedDepthCallbackFunction(
            pImpl->mFixedDepthCallbackFunction);
        pImpl->mLeastSquaresObjectiveFunction->setTopographyCallbackFunction(
            pImpl->mTopographyCallbackFunction);
        if (pImpl->mTravelTimeCalculatorMap)
        {
            pImpl->mLeastSquaresObjectiveFunction->setTravelTimeCalculatorMap(
                std::move(pImpl->mTravelTimeCalculatorMap));
            pImpl->mTravelTimeCalculatorMap = nullptr;
        }
        pImpl->mObjectiveFunction = pImpl->mLeastSquaresObjectiveFunction.get();
        pImpl->mHaveObjectiveFunction = true;
    }
    else if (options.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::L1)
    {
        pImpl->mLogger->debug("Will use L1 objective function...");
        pImpl->mL1ObjectiveFunction
            = std::make_unique<ULocator::ObjectiveFunctions::L1> (pImpl->mLogger);
        pImpl->mL1ObjectiveFunction->setUTMZone(
            options.getUTMZone(),
            options.isNorthernHemisphere());
        pImpl->mL1ObjectiveFunction->setFixedDepthCallbackFunction(
            pImpl->mFixedDepthCallbackFunction);
        pImpl->mL1ObjectiveFunction->setTopographyCallbackFunction(
            pImpl->mTopographyCallbackFunction);
        if (pImpl->mTravelTimeCalculatorMap)
        {
            pImpl->mL1ObjectiveFunction->setTravelTimeCalculatorMap(
                std::move(pImpl->mTravelTimeCalculatorMap));
            pImpl->mTravelTimeCalculatorMap = nullptr;
        }
        pImpl->mObjectiveFunction = pImpl->mL1ObjectiveFunction.get();
        pImpl->mHaveObjectiveFunction = true;
    }
    else if (options.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::LP)
    {
double pNorm = 1.5;
        pImpl->mLogger->debug("Will use LP objective function with p = "
                            + std::to_string(pNorm) + "...");
        pImpl->mLPObjectiveFunction
            = std::make_unique<ULocator::ObjectiveFunctions::LP> (pImpl->mLogger, pNorm);
        pImpl->mLPObjectiveFunction->setUTMZone(
            options.getUTMZone(),
            options.isNorthernHemisphere());
        pImpl->mLPObjectiveFunction->setFixedDepthCallbackFunction(
            pImpl->mFixedDepthCallbackFunction);
        pImpl->mLPObjectiveFunction->setTopographyCallbackFunction(
            pImpl->mTopographyCallbackFunction);
        if (pImpl->mTravelTimeCalculatorMap)
        {
            pImpl->mLPObjectiveFunction->setTravelTimeCalculatorMap(
                std::move(pImpl->mTravelTimeCalculatorMap));
            pImpl->mTravelTimeCalculatorMap = nullptr;
        }
        pImpl->mObjectiveFunction = pImpl->mLPObjectiveFunction.get();
        pImpl->mHaveObjectiveFunction = true;
    }
    else if (options.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::DoubleDifferenceL1)
    {
        pImpl->mLogger->debug(
            "Will use double difference L1 objective function...");
        pImpl->mDoubleDifferenceL1ObjectiveFunction
            = std::make_unique<ULocator::ObjectiveFunctions::DoubleDifferenceL1>
              (pImpl->mLogger);
        pImpl->mDoubleDifferenceL1ObjectiveFunction->setUTMZone(
            options.getUTMZone(),
            options.isNorthernHemisphere());
        pImpl->mDoubleDifferenceL1ObjectiveFunction
            ->setFixedDepthCallbackFunction(
                 pImpl->mFixedDepthCallbackFunction);
        pImpl->mDoubleDifferenceL1ObjectiveFunction
            ->setTopographyCallbackFunction(
                 pImpl->mTopographyCallbackFunction);
        if (pImpl->mTravelTimeCalculatorMap)
        {
            pImpl->mDoubleDifferenceL1ObjectiveFunction
               ->setTravelTimeCalculatorMap(
                   std::move(pImpl->mTravelTimeCalculatorMap));
            pImpl->mTravelTimeCalculatorMap = nullptr;
        }
        pImpl->mObjectiveFunction
            = pImpl->mDoubleDifferenceL1ObjectiveFunction.get();
        pImpl->mHaveObjectiveFunction = true;
    }    
    else
    {
        throw std::runtime_error("Unhandled objective function");
    }
    pImpl->mOptions = options;
    pImpl->mQuarries = options.getQuarries();
    pImpl->mUTMZone = options.getUTMZone();
    pImpl->mNorth = options.isNorthernHemisphere();
    pImpl->mDefaultElevation = pImpl->mOptions.getDefaultElevation();
}

void NLOpt::clearArrivals() noexcept
{
    pImpl->mHaveOrigin = false;
    if (pImpl->mOptions.getObjectiveFunction() ==
        NLOptOptions::ObjectiveFunction::LeastSquares)
    {
        //pImpl->mLeastSquaresObjectiveFunction()->clearArrivals();
    } 
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::L1)
    {

    }
}

void NLOpt::setArrivals(const std::vector<Arrival> &arrivals)
{
    pImpl->mHaveOrigin = false;
    if (pImpl->mOptions.getObjectiveFunction() ==
        NLOptOptions::ObjectiveFunction::LeastSquares)
    {
        pImpl->mLogger->debug(
           "Setting arrivals in least-squares objective function...");
    }
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::L1)
    {
        pImpl->mLogger->debug(
           "Setting arrivals in L1 objective function...");
    }
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::LP)
    {
        pImpl->mLogger->debug(
           "Setting arrivals in LP objective function...");
    }
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::DoubleDifferenceL1)
    {
        pImpl->mLogger->debug(
           "Setting arrivals in double difference L1 objective function...");
    }
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::DoubleDifferenceLeastSquares)
    {
        pImpl->mLogger->debug(
           "Setting arrivals in least-squares double difference objective function...");
    }   
    else
    {
        throw std::runtime_error("Objective function (options) not set");
    }
    pImpl->mObjectiveFunction->setArrivals(arrivals);
}

/// Sets the travel time calculator map
void NLOpt::setTravelTimeCalculatorMap(
    std::unique_ptr<const TravelTimeCalculatorMap> &&calculator)
{
    pImpl->mHaveOrigin = false;
    if (!pImpl->mHaveObjectiveFunction)
    {
        pImpl->mLogger->debug("Temporarily storing calculators");
        pImpl->mTravelTimeCalculatorMap = std::move(calculator);
    }
    else
    {
        if (pImpl->mOptions.getObjectiveFunction() ==
            NLOptOptions::ObjectiveFunction::LeastSquares)
        {
            pImpl->mLogger->debug(
                "Setting travel time calculators on least squares");
        }
        else if (pImpl->mOptions.getObjectiveFunction() ==
                 NLOptOptions::ObjectiveFunction::L1)
        {
            pImpl->mLogger->debug(
                "Setting travel time calculators on l1");
        }
        else if (pImpl->mOptions.getObjectiveFunction() ==
                 NLOptOptions::ObjectiveFunction::LP)
        {
            pImpl->mLogger->debug(
                "Setting travel time calculators on lp");
        }
        else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::DoubleDifferenceL1)
        {
            pImpl->mLogger->debug(
               "Setting travel time calculators on double difference L1");
        }
        else if (pImpl->mOptions.getObjectiveFunction() ==
                 NLOptOptions::ObjectiveFunction::DoubleDifferenceLeastSquares)
        {
            pImpl->mLogger->debug(
               "Setting travel time calculators on double difference L2");
        }
        else
        {
            pImpl->mLogger->error("Unhandled objective function");
        }
#ifndef NDEBUG
        assert(pImpl->mObjectiveFunction != nullptr);
#endif 
        pImpl->mObjectiveFunction->setTravelTimeCalculatorMap(
            std::move(calculator));

    }
}

std::unique_ptr<const TravelTimeCalculatorMap> 
NLOpt::releaseTravelTimeCalculatorMap()
{
    if (!pImpl->mHaveObjectiveFunction)
    {
        auto result = std::move(pImpl->mTravelTimeCalculatorMap);
        pImpl->mTravelTimeCalculatorMap = nullptr;
        return result;
    }
    return pImpl->mObjectiveFunction->releaseTravelTimeCalculatorMap();
}

/// Sets the topography
void NLOpt::setTopography(
    std::unique_ptr<const Topography> &&topography)
{
    pImpl->mHaveOrigin = false;
    pImpl->mLogger->debug("Setting topography...");
    pImpl->mTopography = std::move(topography);
}

std::unique_ptr<const Topography> NLOpt::releaseTopography()
{
    if (pImpl->mTopography == nullptr)
    {
        throw std::runtime_error("Topography not set");
    }
    auto result = std::move(pImpl->mTopography);
    pImpl->mTopography = nullptr;
    return result;
}

Origin NLOpt::getOrigin() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mOrigin;
}

bool NLOpt::haveOrigin() const noexcept
{
    return pImpl->mHaveOrigin;
}

void NLOpt::locateQuarryBlast()
{
    Origin emptyOrigin;
    pImpl->locate(SourceDepthConstraint::FixedToFreeSurface, emptyOrigin, true);
}

/*
void NLOpt::locateQuarryBlast(const Origin &origin)
{
    auto depth = pImpl->mOptions.getInitialQuarryBlastSearchDepth();
    //pImpl->locate(SourceDepthConstraint::FixedToFreeSurface, depth, true);
}
*/

Origin NLOpt::predict(const Origin &origin, const bool applyCorrection)
{
    auto oldApplyCorrection = pImpl->mObjectiveFunction->applyCorrection();
    pImpl->mObjectiveFunction->toggleApplyCorrection(applyCorrection);
    Origin newOrigin;
    try
    {
        newOrigin = pImpl->mObjectiveFunction->predict(origin);
        pImpl->mObjectiveFunction->toggleApplyCorrection(oldApplyCorrection);
        return newOrigin;
    }
    catch (const std::exception &e)
    {
        pImpl->mObjectiveFunction->toggleApplyCorrection(oldApplyCorrection); 
        throw;
    } 
}

double NLOpt::evaluateLoss(const Origin &origin) const
{
    if (!pImpl->mHaveObjectiveFunction)
    {
        throw std::runtime_error("Objective function not set");
    }
    if (!origin.haveEpicenter())
    {
        throw std::invalid_argument("Epicenter not set");
    }
    if (!origin.haveDepth())
    {
        throw std::invalid_argument("Depth not set");
    }
    if (!origin.haveTime())
    {
        throw std::invalid_argument("Origin time not set");
    }
    return pImpl->mObjectiveFunction->evaluateLoss(origin);
}


void NLOpt::locateEarthquake()
{
    Origin emptyOrigin;
    pImpl->locate(SourceDepthConstraint::Free, emptyOrigin, false);
}

/*
void NLOpt::locateEarthquake(const Origin &origin)
{
    double depth = pImpl->mOptions.getInitialEarthquakeSearchDepth();
    if (origin.haveDepth()){depth = origin.getDepth();}
    locate(SourceDepthConstraint::Free, depth);
}
*/

/*
void NLOpt::locate(const SourceDepthConstraint sourceConstraint,
                   const double depth) 
{
    if (!pImpl->mHaveObjectiveFunction)
    {
        throw std::runtime_error("Objective function (options) not set");
    }

    ObjectiveFunctions::IObjectiveFunction *objectiveFunction{nullptr};
    if (pImpl->mOptions.getObjectiveFunction() ==
        NLOptOptions::ObjectiveFunction::LeastSquares)
    {
        objectiveFunction = pImpl->mLeastSquaresObjectiveFunction.get();
    }
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::L1)
    {
        objectiveFunction = pImpl->mL1ObjectiveFunction.get();
    }
    else if (pImpl->mOptions.getObjectiveFunction() ==
             NLOptOptions::ObjectiveFunction::LeastSquaresDoubleDifference)
    {
        //objectiveFunction = pImpl->mLeastSquaresDoubleDifferenceObjectiveFunction.get();
    }
    else
    {
        throw std::runtime_error("Unhandled objective function");
    }
    // Set the inversion strategy 
    if (sourceConstraint == SourceDepthConstraint::Free)
    {
        objectiveFunction->setInversionStrategy(
           InversionStrategy::LatitudeLongitudeDepth);
    }
    else if (sourceConstraint == SourceDepthConstraint::Fixed)
    {
        objectiveFunction->setInversionStrategy(
           InversionStrategy::LatitudeLongitudeFixedDepth);
    }
    else if (sourceConstraint == SourceDepthConstraint::FixedToFreeSurface)
    {
        objectiveFunction->setInversionStrategy(
            InversionStrategy::LatitudeLongitudeDepthAtFreeSurface);
        objectiveFunction->setTopographyCallbackFunction(
            pImpl->mTopographyCallbackFunction);
    }
    else
    {
        throw std::runtime_error("Unhandled source depth constraint");
    }
    // Set up some geographic information
    pImpl->mNorth = pImpl->mOptions.isNorthernHemisphere();
    // Set up the optimization scheme
    objectiveFunction->toggleApplyCorrection(true); 
    objectiveFunction->resetCounters();
    auto nParameters = objectiveFunction->getNumberOfModelParameters();
    auto nlOptObjectiveFunction
         = std::bind(&ObjectiveFunctions::IObjectiveFunction::operator(),
                     objectiveFunction,
                     std::placeholders::_1,
                     std::placeholders::_2,
                     std::placeholders::_3);
    // Get the boundaries
    auto latitudeBoundaries  = pImpl->mOptions.getLatitudeBoundaries();
    auto longitudeBoundaries = pImpl->mOptions.getLongitudeBoundaries();
    auto depthBoundaries     = pImpl->mOptions.getDepthBoundaries();
    auto [minimumUTMX, minimumUTMY] = pImpl->toUTM(latitudeBoundaries[0],
                                                   longitudeBoundaries[0]);
    auto [maximumUTMX, maximumUTMY] = pImpl->toUTM(latitudeBoundaries[1],
                                                   longitudeBoundaries[1]);
    // Perform the initial optimization
    nlopt::opt initialOptimizer(nlopt::GN_DIRECT_L, nParameters);
//nlopt::opt optimizer(nlopt::LN_BOBYQA, nParameters);
//nlopt::opt optimizer(nlopt::LD_SLSQP, nParameters);
//nlopt::opt optimizer(nlopt::GD_STOGO, nParameters);
    initialOptimizer.set_min_objective(nlOptObjectiveFunction);
    initialOptimizer.set_maxeval(750); //pImpl->mOptions.getMaximumNumberOfFunctionEvaluations());
    initialOptimizer.set_lower_bounds(
        std::vector<double> {minimumUTMX, minimumUTMY, depthBoundaries[0]});
    initialOptimizer.set_upper_bounds(
        std::vector<double> {maximumUTMX, maximumUTMY, depthBoundaries[1]});
initialOptimizer.set_xtol_abs(1.e-1);
std::vector<double> xLatLon{38.1,    -112.1,       0};
    auto [initialSourceUTMX, initialSourceUTMY] = pImpl->toUTM(xLatLon[0], xLatLon[1]);
    std::vector<double> xInitialLocation{xLatLon};
    xInitialLocation[0] = initialSourceUTMX;
    xInitialLocation[1] = initialSourceUTMY;
    double initialOptimalValue = 0;
    try
    {
        pImpl->mLogger->debug("Performing initial location...");
        initialOptimizer.optimize(xInitialLocation, initialOptimalValue);
    }
    catch (const std::exception &e)
    {
        pImpl->mLogger->error(e.what());
        throw std::runtime_error("Initial optimization failed");
    }
    auto initialDepth = xInitialLocation[2];
    auto [initialLatitude, initialLongitude] = pImpl->toLatLon(xInitialLocation[0],
                                                               xInitialLocation[1]);
    pImpl->mLogger->debug("Initial (lat,lon,depth,time,value) = ("
                       +       std::to_string(initialLatitude)
                       + "," + std::to_string(initialLongitude)
                       + "," + std::to_string(initialDepth)
                       + "," + std::to_string(initialOptimalValue)
                       + ")");
    // Now perform the optimization over a few depths
    objectiveFunction->resetCounters();
    nlopt::opt optimizer(nlopt::LD_MMA, nParameters);
    optimizer.set_min_objective(nlOptObjectiveFunction);
    optimizer.set_maxeval(1500); //pImpl->mOptions.getMaximumNumberOfFunctionEvaluations());
    optimizer.set_xtol_abs(1.e-1);
    std::vector<double> x(nParameters);
    x.at(0) = initialSourceUTMX;
    x.at(1) = initialSourceUTMY;
    x.at(2) = initialDepth;
    auto latitudeRefinement  = pImpl->mOptions.getLatitudeRefinement();
    auto longitudeRefinement = pImpl->mOptions.getLongitudeRefinement();
    auto refinedMinimumUTMX = initialSourceUTMX - latitudeRefinement;
    auto refinedMaximumUTMX = initialSourceUTMX + latitudeRefinement;
    auto refinedMinimumUTMY = initialSourceUTMY - longitudeRefinement;
    auto refinedMaximumUTMY = initialSourceUTMY + longitudeRefinement;
    optimizer.set_lower_bounds(
        std::vector<double> {refinedMinimumUTMX,
                             refinedMinimumUTMY,
                             depthBoundaries[0]});
    optimizer.set_upper_bounds(
        std::vector<double> {refinedMaximumUTMX,
                             refinedMaximumUTMY,
                             depthBoundaries[1]});
 
std::vector<double> trialDepths{-1000, 0, 1000, 5000, 10000, 20000};
    std::vector<double> finalLocation{x};
    double bestOptimalValue{std::numeric_limits<double>::max()};
    for (auto candidateDepth : trialDepths)
    {
        double optimalValue{0};
        x.at(0) = initialSourceUTMX;
        x.at(1) = initialSourceUTMY;
        x.at(2) = candidateDepth;
        try
        {
            pImpl->mLogger->debug(
                "Performing final location search with initial depth: "
              + std::to_string(candidateDepth));
            optimizer.optimize(x, optimalValue);
            pImpl->mLogger->debug("Optimal value: "
                                + std::to_string(optimalValue));
        }
        catch (const std::exception &e)
        {
            pImpl->mLogger->error(e.what());
            throw std::runtime_error("Final optimization failed");
        }
        if (optimalValue < bestOptimalValue &&
            optimalValue < initialOptimalValue)
        {
            bestOptimalValue = optimalValue;
            finalLocation = x;
        }
    } 
    auto [finalLatitude, finalLongitude] = pImpl->toLatLon(xInitialLocation[0],
                                                           xInitialLocation[1]);
    auto finalDepth = finalLocation[2];
    pImpl->mLogger->debug("Final (lat,lon,depth,time,value) = ("
                       +       std::to_string(finalLatitude)
                       + "," + std::to_string(finalLongitude)
                       + "," + std::to_string(finalDepth)
                       + "," + std::to_string(bestOptimalValue)
                       + ")");

}
*/

int NLOpt::getNumberOfObjectiveFunctionEvaluations() const noexcept
{
    return pImpl->mObjectiveFunctionCounter;
}

int NLOpt::getNumberOfGradientEvaluations() const noexcept
{
    return pImpl->mGradientCounter;
}

