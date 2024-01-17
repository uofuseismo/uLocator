#include <string>
#include <vector>
#include <umps/logging/standardOut.hpp>
#include "uLocator/optimizers/optimizer.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"

using namespace ULocator::Optimizers;

class IOptimizer::IOptimizerImpl
{
public:
    IOptimizerImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    ULocator::Origin mOrigin;
    std::vector<ULocator::Arrival> mArrivals;
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::unique_ptr<ULocator::TravelTimeCalculatorMap> mTravelTimeCalculatorMap{nullptr};
    std::unique_ptr<ULocator::Topography::ITopography> mTopography{nullptr};
    std::unique_ptr<ULocator::Position::IGeographicRegion> mGeographicRegion{nullptr};
    bool mHaveOrigin{false};
};

/// Constructor
IOptimizer::IOptimizer() :
    pImpl(std::make_unique<IOptimizerImpl> ())
{
}

/// Constructor
IOptimizer::IOptimizer(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<IOptimizerImpl> (logger))
{
}

/// Set arrivals
void IOptimizer::setArrivals(const std::vector<ULocator::Arrival> &arrivals)
{
    if (arrivals.empty())
    {
        throw std::invalid_argument("No arrivals!");
    }
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator map not set");
    }
    pImpl->mArrivals.clear();
    pImpl->mArrivals.reserve(arrivals.size());
    std::vector<std::string> namePhases;
    namePhases.reserve(arrivals.size());
    for (const auto &arrival : arrivals)
    {
        if (!arrival.haveStation())
        {
            pImpl->mLogger->warn(
                "Arrival's station not set.  Skipping...");
            continue;
        }
        auto station = arrival.getStation();
        auto name = station.getNetwork() + "." + station.getName();
        if (!arrival.havePhase())
        {
            pImpl->mLogger->warn("Arrival's phase not set on " + name
                               + ".  Skipping...");
            continue;
        }
        auto phase = arrival.getPhase();
        auto namePhase = name + "." + phase;
        if (!arrival.haveTime())
        {
            pImpl->mLogger->warn("Arrival time not set on "
                               + namePhase + ".  Skipping...");
            continue;
        }
        if (!arrival.haveStandardError())
        {
            pImpl->mLogger->warn("Standard error not set on "
                              + namePhase + ".  Skipping...");
            continue;
        } 
        if (!pImpl->mTravelTimeCalculatorMap->contains(station, phase))
        {
            pImpl->mLogger->warn("No travel time calculator for "
                               + namePhase
                               + ".  Skipping...");
            continue;
        }
        bool duplicate{false};
        for (const auto &np : namePhases)
        {
            if (namePhase == np)
            {
                pImpl->mLogger->warn(np + " already exists.  Skipping...");
                duplicate = true;
                break;
            }
        }
        if (duplicate){continue;}
        // Add it
        pImpl->mArrivals.push_back(arrival);
        namePhases.push_back(namePhase);
    } 
    pImpl->mLogger->debug("Set " + std::to_string(pImpl->mArrivals.size())
                        + " phase arrivals...");
}

const std::vector<ULocator::Arrival> 
&IOptimizer::getArrivalsReference() const noexcept
{
    return *&pImpl->mArrivals;
}


/// Geographic region
void IOptimizer::setGeographicRegion(
    const ULocator::Position::IGeographicRegion &region)
{
    pImpl->mGeographicRegion = region.clone();
}

std::unique_ptr<ULocator::Position::IGeographicRegion>
IOptimizer::getGeographicRegion() const
{
    if (!haveGeographicRegion())
    {
        throw std::runtime_error("Geographic region not set");
    }
    return pImpl->mGeographicRegion->clone();
}

bool IOptimizer::haveGeographicRegion() const noexcept
{
    return pImpl->mGeographicRegion != nullptr;
}

/// Release topography
void IOptimizer::setTopography(
    std::unique_ptr<ULocator::Topography::ITopography> &&topography)
{
    if (topography == nullptr)
    {
        throw std::invalid_argument("Topography is NULL");
    }
    pImpl->mTopography = std::move(topography);
}

const ULocator::Topography::ITopography*
    IOptimizer::getTopography() const
{
    if (!haveTopography()){throw std::runtime_error("Topography not set");}
    return pImpl->mTopography.get();
}

std::unique_ptr<ULocator::Topography::ITopography>
IOptimizer::releaseTopography()
{
    if (!haveTopography()){throw std::runtime_error("Topography not set");}
    std::unique_ptr<ULocator::Topography::ITopography>
        result{pImpl->mTopography.release()};
    pImpl->mTopography = nullptr;
    return result;
}

bool IOptimizer::haveTopography() const noexcept
{
    return (pImpl->mTopography != nullptr);
}

/// Travel time calculator map
void IOptimizer::setTravelTimeCalculatorMap(
    std::unique_ptr<ULocator::TravelTimeCalculatorMap> &&map)
{
    if (map == nullptr)
    {
        throw std::invalid_argument("Travel time calculator map is NULL");
    }
    pImpl->mTravelTimeCalculatorMap = std::move(map);
}

const ULocator::TravelTimeCalculatorMap 
    *IOptimizer::getTravelTimeCalculatorMap() const
{
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator map not set");
    }
    return pImpl->mTravelTimeCalculatorMap.get();
}

std::unique_ptr<ULocator::TravelTimeCalculatorMap>
IOptimizer::releaseTravelTimeCalculatorMap()
{
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator map not set");
    }
    std::unique_ptr<ULocator::TravelTimeCalculatorMap>
        result{pImpl->mTravelTimeCalculatorMap.release()};
    pImpl->mTravelTimeCalculatorMap = nullptr;
    return result;
}
 
bool IOptimizer::haveTravelTimeCalculatorMap() const noexcept
{
    return (pImpl->mTravelTimeCalculatorMap != nullptr);
}

/// Destructor
IOptimizer::~IOptimizer() = default;

/// The origin
void IOptimizer::setOrigin(const ULocator::Origin &origin)
{
    auto work = origin;
    setOrigin(std::move(work));
}

void IOptimizer::setOrigin(ULocator::Origin &&origin)
{
    if (!origin.haveTime())
    {
        throw std::invalid_argument("Origin time not set");
    }
    if (!origin.haveEpicenter())
    {
        throw std::invalid_argument("Origin epicenter not set");
    }
    if (!origin.haveDepth())
    {
        throw std::invalid_argument("Origin depth not set");
    }
    pImpl->mOrigin = std::move(origin);
}

ULocator::Origin IOptimizer::getOrigin() const
{
    if (!haveOrigin()){throw std::runtime_error("Origin not set");}
    return pImpl->mOrigin;
}

bool IOptimizer::haveOrigin() const noexcept
{
    return pImpl->mHaveOrigin;
}
