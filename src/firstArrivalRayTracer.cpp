#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#ifdef WITH_EIKONALXX
#include <eikonalxx/ray/layerSolver.hpp>
#include <eikonalxx/ray/path2d.hpp>
#include <eikonalxx/ray/segment2d.hpp>
#include <eikonalxx/ray/point2d.hpp>
#else
#include "ray/layerSolver.hpp"
#include "ray/path2d.hpp"
#include "ray/segment2d.hpp"
#include "ray/point2d.hpp"
#endif
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include "uLocator/firstArrivalRayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/staticCorrection.hpp"
#include "uLocator/sourceSpecificStationCorrection.hpp"
#include "uLocator/position/wgs84.hpp"

using namespace ULocator;

class FirstArrivalRayTracer::FirstArrivalRayTracerImpl
{
public:
    FirstArrivalRayTracerImpl(
        std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    double evaluateLayerSolver(const double sourceDepth, const double offset,
                               double *dtdx, double *dtdz) const
    { 
        constexpr bool doReflections{false};
        double travelTime{0};
        if (dtdx != nullptr){*dtdx = 0;} 
        if (dtdz != nullptr){*dtdz = 0;} 
        try
        {
            auto sourceDepthToUse = sourceDepth;
            if (sourceDepth < mMinimumDepth)
            {
                sourceDepthToUse = mMinimumDepth;
                mLogger->warn("Overriding source depth "
                            + std::to_string(sourceDepth)
                            + " to " + std::to_string(mMinimumDepth));
            }
            mSolver.setStationOffsetAndDepth(offset, mStationDepth);
            mSolver.setSourceDepth(sourceDepthToUse);
            mSolver.solve(doReflections);
        }
        catch (const std::exception &e) 
        {
            mLogger->error("Failed to compute ray paths for "
                         + mStation.getNetwork() + "." 
                         + mStation.getName() + "." 
                         + mPhase
                         + ".  Failed with: "
                         + std::string {e.what()}
                         +  ".  Source depth is: " 
                         + std::to_string(sourceDepth));
            return travelTime;
        }
        const auto &rayPaths = mSolver.getRayPaths(); 
        if (!rayPaths.empty())
        {
            travelTime = rayPaths[0].getTravelTime();
            // Decompose segments in first leg into slowness in x and z
            // as these are the sensitivities at the source.
            if (dtdx != nullptr || dtdz != nullptr)
            {
                auto takeOffAngle = rayPaths[0].getTakeOffAngle()*(M_PI/180);
                auto slowness = rayPaths[0].at(0).getSlowness();
                if (dtdx != nullptr){*dtdx = slowness*std::sin(takeOffAngle);}
                if (dtdz != nullptr){*dtdz = slowness*std::cos(takeOffAngle);}
            }
        }
        return travelTime;
    }
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    SourceSpecificStationCorrection mSourceSpecificStationCorrection;
    StaticCorrection mStaticCorrection;
    mutable EikonalXX::Ray::LayerSolver mSolver;
    Station mStation;
    std::string mPhase;
    std::vector<double> mInterfaces;
    std::vector<double> mVelocities;
    double mStationUTMX{0};
    double mStationUTMY{0};
    double mStationDepth{0};
    double mMinimumDepth{0};
    bool mInitialized{false};
};

/// Constructor
FirstArrivalRayTracer::FirstArrivalRayTracer() :
    pImpl(std::make_unique<FirstArrivalRayTracerImpl> ())
{
}

/// Constructor with logger
FirstArrivalRayTracer::FirstArrivalRayTracer(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<FirstArrivalRayTracerImpl> (logger))
{
}

/// Move constructor
FirstArrivalRayTracer::FirstArrivalRayTracer(
    FirstArrivalRayTracer &&rayTracer) noexcept
{
    *this = std::move(rayTracer);
}

/// Reset class
void FirstArrivalRayTracer::clear() noexcept
{
    pImpl->mSolver.clear();
    pImpl->mStation.clear();
    pImpl->mPhase.clear();
    pImpl->mInterfaces.clear();
    pImpl->mVelocities.clear();
    pImpl->mStationDepth = 0;
    pImpl->mMinimumDepth = 0;
    pImpl->mInitialized = false;
}

/// Destructor
FirstArrivalRayTracer::~FirstArrivalRayTracer() = default;

/// Move assignment
FirstArrivalRayTracer&
FirstArrivalRayTracer::operator=(FirstArrivalRayTracer &&rayTracer) noexcept
{
    if (&rayTracer == this){return *this;}
    pImpl = std::move(rayTracer.pImpl);
    return *this;
}

/// Initialized?
bool FirstArrivalRayTracer::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Evaluate
double FirstArrivalRayTracer::evaluate(
    const Position::WGS84 &epicenter,
    const double depth,
    const bool applyCorrection) const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Ray tracer class not initialized");
    }
    if (!epicenter.havePosition())
    {
        throw std::invalid_argument("Latitude/longitude not set");
    }
    double epicentralDistance, distance;
    ULocator::Position::computeDistanceAzimuth(
        epicenter,
        pImpl->mStation.getGeographicPositionReference(),
        &epicentralDistance,
        &distance,
        nullptr,
        nullptr);
    auto travelTime = pImpl->evaluateLayerSolver(depth, distance,
                                                 nullptr, nullptr);
    if (applyCorrection)
    {
        if (pImpl->mStaticCorrection.haveCorrection())
        {
            travelTime = pImpl->mStaticCorrection.evaluate(travelTime);
        }
        if (pImpl->mSourceSpecificStationCorrection.haveModel())
        {
            auto correction
                = pImpl->mSourceSpecificStationCorrection.evaluate(
                     epicenter.getLatitude(),
                     epicenter.getLongitude(),
                     depth,
                     SourceSpecificStationCorrection::
                        EvaluationMethod::InverseDistanceWeighted);
            travelTime = travelTime + correction;
        }
    }
    return travelTime;
}

/*
/// Evaluate
double FirstArrivalRayTracer::evaluate(
    const double utmX,
    const double utmY,
    const double depth,
    double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    if (!isInitialized())
    {   
        throw std::runtime_error("Ray tracer class not initialized");
    }   
    auto dx = pImpl->mStationUTMX - utmX;
    auto dy = pImpl->mStationUTMY - utmY;
    auto azimuth = std::atan2(dy, dx); 
    double distance = std::hypot(dx, dy);
    double dtdr;
    auto travelTime = pImpl->evaluateLayerSolver(depth, distance, &dtdr, dtdz);
    constexpr double degreesToRadians{M_PI/180.};
    *dtdx = dtdr*std::cos(azimuth);
    *dtdy = dtdr*std::sin(azimuth);
    if (applyCorrection)
    {
        if (pImpl->mStaticCorrection.haveCorrection())
        {
            travelTime = pImpl->mStaticCorrection.evaluate(travelTime);
        }
    }
    return travelTime;
}
*/

/// Evaluate
double FirstArrivalRayTracer::evaluate(
    const Position::WGS84 &epicenter,
    const double depth,
    double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Ray tracer class not initialized");
    }
    if (!epicenter.havePosition())
    {   
        throw std::invalid_argument("Latitude/longitude not set");
    }   
    double epicentralDistance, distance, azimuth;
    ULocator::Position::computeDistanceAzimuth(
        epicenter,
        pImpl->mStation.getGeographicPositionReference(),
        &epicentralDistance,
        &distance,
        &azimuth,
        nullptr);
    double dtdr;
    auto travelTime = pImpl->evaluateLayerSolver(depth, distance, &dtdr, dtdz);
    constexpr double degreesToRadians{M_PI/180.};
    if (dtdx != nullptr){*dtdx = dtdr*std::sin(azimuth*degreesToRadians);}
    if (dtdy != nullptr){*dtdy = dtdr*std::cos(azimuth*degreesToRadians);}
    if (applyCorrection)
    {
        if (pImpl->mStaticCorrection.haveCorrection())
        {
            travelTime = pImpl->mStaticCorrection.evaluate(travelTime);
        }
        if (pImpl->mSourceSpecificStationCorrection.haveModel())
        {
            auto correction
                = pImpl->mSourceSpecificStationCorrection.evaluate(
                     epicenter.getLatitude(),
                     epicenter.getLongitude(),
                     depth,
                     SourceSpecificStationCorrection::
                        EvaluationMethod::InverseDistanceWeighted);
            travelTime = travelTime + correction;
        }
    }
    return travelTime;
}

/// Predefined models
void FirstArrivalRayTracer::initialize(
    const Station &station,
    const std::string &phase,
    const std::vector<double> &interfaces,
    const std::vector<double> &velocities)
{
    clear();
    if (phase.empty()){throw std::invalid_argument("Phase is empty");}
    if (!station.haveGeographicPosition())
    {
        throw std::invalid_argument("Station position not set");
    }
    if (!station.haveElevation())
    {
        throw std::invalid_argument("Station elevation not set");
    }
    if (interfaces.empty()){throw std::invalid_argument("Interfaces is empty");}
    if (velocities.empty()){throw std::invalid_argument("Velocities is empty");}
    pImpl->mSolver.setVelocityModel(interfaces, velocities); 
    pImpl->mMinimumDepth = interfaces.at(0);
    pImpl->mStationDepth =-station.getElevation();
    pImpl->mStation = station;
    pImpl->mStationUTMX = station.getGeographicPosition().getEasting();
    pImpl->mStationUTMY = station.getGeographicPosition().getNorthing();
    pImpl->mPhase = phase;
    pImpl->mInterfaces = interfaces;
    pImpl->mVelocities = velocities;
    pImpl->mInitialized = true;
}

void FirstArrivalRayTracer::initialize(
    const Station &station,
    const std::string &phase,
    const std::vector<double> &interfaces,
    const std::vector<double> &velocities,
    const StaticCorrection &staticCorrection)
{
    initialize(station, phase, interfaces, velocities);
    if (staticCorrection.haveCorrection())
    {
        pImpl->mStaticCorrection = staticCorrection;
    }
}

void FirstArrivalRayTracer::initialize(
    const Station &station,
    const std::string &phase,
    const std::vector<double> &interfaces,
    const std::vector<double> &velocities,
    const SourceSpecificStationCorrection &sssc)
{
    initialize(station, phase, interfaces, velocities);
    if (sssc.haveModel())
    {   
        pImpl->mSourceSpecificStationCorrection = sssc;
    }   
}


void FirstArrivalRayTracer::initialize(
    const Station &station,
    const std::string &phase,
    const std::vector<double> &interfaces,
    const std::vector<double> &velocities,
    const StaticCorrection &staticCorrection,
    const SourceSpecificStationCorrection &sssc)
{
    initialize(station, phase, interfaces, velocities);
    if (staticCorrection.haveCorrection())
    {   
        pImpl->mStaticCorrection = staticCorrection;
    }
    if (sssc.haveModel())
    {   
        pImpl->mSourceSpecificStationCorrection = sssc;
    }   
}

void FirstArrivalRayTracer::initialize(
    const Station &station,
    const std::string &phase,
    const std::vector<double> &interfaces,
    const std::vector<double> &velocities,
    StaticCorrection &&staticCorrection,
    SourceSpecificStationCorrection &&sssc)
{
    initialize(station, phase, interfaces, velocities);
    if (staticCorrection.haveCorrection())
    {   
        pImpl->mStaticCorrection = std::move(staticCorrection);
    }   
    if (sssc.haveModel())
    {   
        pImpl->mSourceSpecificStationCorrection = std::move(sssc);
    }   
}


void FirstArrivalRayTracer::initializeUtahP(const Station &station,
                                            const bool useAlternateModel)
{
    const std::vector<double> interfaces{-4500,   40,  15600, 26500, 40500};
    std::vector<double> velocities{       3400, 5900,  6400,  7500,  7900};
    if (useAlternateModel)
    {
        //velocities = std::vector<double> {3400, 5950, 6450,  7550,  7900};
        velocities = std::vector<double> {3685, 5972, 6553, 7106, 7672};
    }
    initialize(station, "P", interfaces, velocities);
}

void FirstArrivalRayTracer::initializeUtahP(
    const Station &station,
    StaticCorrection &&staticCorrection,
    SourceSpecificStationCorrection &&sssc,
    const bool useAlternateModel)
{
    initializeUtahP(station, useAlternateModel);
    if (staticCorrection.haveCorrection())
    {
        pImpl->mStaticCorrection = std::move(staticCorrection);
    }   
    if (sssc.haveModel())
    {
        pImpl->mSourceSpecificStationCorrection = std::move(sssc);
    }
}


void FirstArrivalRayTracer::initializeUtahS(const Station &station,
                                            const bool useAlternateModel)
{
    const std::vector<double> interfaces{-4500,  40,  15600, 26500, 40500};
    std::vector<double> velocities{1950, 3390, 3680,   4310,  4540};
    if (useAlternateModel)
    {
        //velocities = std::vector<double> {2000, 3425, 3700,   4400,  4550};
        velocities = std::vector<double> {2169, 3443, 3653, 4486, 5022};
    }
    initialize(station, "S", interfaces, velocities);
}

void FirstArrivalRayTracer::initializeUtahS(
    const Station &station,
    StaticCorrection &&staticCorrection,
    SourceSpecificStationCorrection &&sssc,
    const bool useAlternateModel)
{
    initializeUtahS(station, useAlternateModel);
    if (staticCorrection.haveCorrection())
    {   
        pImpl->mStaticCorrection = std::move(staticCorrection);
    }   
    if (sssc.haveModel())
    {   
        pImpl->mSourceSpecificStationCorrection = std::move(sssc);
    }   
}

void FirstArrivalRayTracer::initializeYellowstoneP(const Station &station,
                                                   const bool useAlternateModel)
{
    const std::vector<double> interfaces{-4500, -1000,  2000,  5000,  8000,
                                         12000, 16000, 21000, 50000};
    std::vector<double> velocities{2720, 2790, 5210, 5560, 5770,
                                   6070, 6330, 6630, 8000};
    if (useAlternateModel)
    {
        velocities = std::vector<double> {2512, 3398, 4689, 5456, 5674,
                                          6250, 6398, 6574, 8200};
    }
    initialize(station, "P", interfaces, velocities);
}

void FirstArrivalRayTracer::initializeYellowstoneP(
    const Station &station,
    StaticCorrection &&staticCorrection,
    SourceSpecificStationCorrection &&sssc,
    const bool useAlternateModel)
{
    initializeYellowstoneP(station, useAlternateModel);
    if (staticCorrection.haveCorrection())
    {   
        pImpl->mStaticCorrection = std::move(staticCorrection);
    }   
    if (sssc.haveModel())
    {   
        pImpl->mSourceSpecificStationCorrection = std::move(sssc);
    }   
}

void FirstArrivalRayTracer::initializeYellowstoneS(const Station &station,
                                                   const bool useAlternateModel)
{
    const std::vector<double> interfaces{-4500, -1000,  2000,  5000,  8000,
                                         12000, 16000, 21000, 50000};
    std::vector<double> velocities{1950, 2000, 3400, 3420, 3490,
                                   3680, 3780, 4000, 4850};
    if (useAlternateModel)
    {   
        velocities = std::vector<double> {1725, 2343, 3064, 3425, 3569,
                                          3690, 3705, 3975, 4950};
    }
    initialize(station, "S", interfaces, velocities);
}

void FirstArrivalRayTracer::initializeYellowstoneS(
    const Station &station,
    StaticCorrection &&staticCorrection,
    SourceSpecificStationCorrection &&sssc,
    const bool useAlternateModel)
{
    initializeYellowstoneS(station, useAlternateModel);
    if (staticCorrection.haveCorrection())
    {
        pImpl->mStaticCorrection = std::move(staticCorrection);
    }
    if (sssc.haveModel())
    {
        pImpl->mSourceSpecificStationCorrection = std::move(sssc);
    }
}

std::vector<double> FirstArrivalRayTracer::getVelocites() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Ray tracer not initialized");
    }
    return pImpl->mVelocities;
}

std::vector<double> FirstArrivalRayTracer::getInterfaces() const
{
    if (!isInitialized())
    {
        throw std::runtime_error("Ray tracer not initialized");
    }
    return pImpl->mInterfaces;
} 
