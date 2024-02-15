#include <iostream>
#include <vector>
#include <mutex>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include "uLocator/uussRayTracer.hpp"
#include "uLocator/rayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/corrections/static.hpp"
#include "uLocator/corrections/sourceSpecific.hpp"

using namespace ULocator;

class UUSSRayTracer::UUSSRayTracerImpl
{
public:
    UUSSRayTracerImpl(const Station &station,
                      const UUSSRayTracer::Phase phase,
                      const UUSSRayTracer::Region region,
                      Corrections::Static &&staticCorrection,
                      Corrections::SourceSpecific &&sourceSpecific,
                      std::shared_ptr<UMPS::Logging::ILog> &logger)
    {
        if (!station.haveNetwork())
        {
            throw std::invalid_argument("Station's network not set");
        }
        if (!station.haveName())
        {
            throw std::invalid_argument("Station name not set");
        }
        if (!station.haveGeographicPosition())
        {
            throw std::invalid_argument("Station geographic position not set");
        }
        if (!station.haveElevation())
        {
            throw std::invalid_argument("Station elevation not set");
        }
        auto stationDepth =-station.getElevation();
        auto interfaces = UUSSRayTracer::getInterfaces(region);
        interfaces[0] = std::min(stationDepth - 1.e-8, interfaces.at(0));
        std::vector<double> velocities;
        std::string stringPhase;
        if (phase == UUSSRayTracer::Phase::P)
        {
            if (logger != nullptr)
            {
                if (region == UUSSRayTracer::Region::Utah)
                {
                    logger->debug("Initializing Utah P layer solver...");
                }
                else
                {
                    logger->debug("Initializing YNP P layer solver...");
                }
            }
            velocities = UUSSRayTracer::getPVelocities(region);
            stringPhase = "P";
        }
        else
        {
            if (logger != nullptr)
            {
                if (region == UUSSRayTracer::Region::Utah)
                {
                    logger->debug("Initializing Utah S layer solver...");
                }
                else
                {
                    logger->debug("Initializing YNP S layer solver...");
                }
            }
            velocities = UUSSRayTracer::getSVelocities(region);
            stringPhase = "S"; 
        }
#ifndef NDEBUG
        assert(interfaces.size() == velocities.size());
#endif
        mRayTracer
            = std::make_unique<ULocator::RayTracer>(station,
                                                    stringPhase,
                                                    interfaces,
                                                    velocities,
                                                    std::move(staticCorrection),
                                                    std::move(sourceSpecific),
                                                    logger);
    } 
    std::unique_ptr<ULocator::RayTracer> mRayTracer{nullptr};
};

UUSSRayTracer::UUSSRayTracer(const Station &station,
                             const UUSSRayTracer::Phase phase,
                             const UUSSRayTracer::Region region,
                             std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<UUSSRayTracerImpl> (
        station,
        phase,
        region,
        std::move(ULocator::Corrections::Static {}),
        std::move(ULocator::Corrections::SourceSpecific {}),
        logger)
    )
{
}

/*
UUSSRayTracer::UUSSRayTracer(const Station &station,
                             const UUSSRayTracer::Phase phase,
                             const UUSSRayTracer::Region region,
                             ULocator::Corrections::Static &&staticCorrection,
                             std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<UUSSRayTracerImpl> (
        station,
        phase,
        region,
        std::move(staticCorrection),
        std::move(ULocator::Corrections::SourceSpecific {}),
        logger)
    )   
{
}
*/

UUSSRayTracer::UUSSRayTracer(const Station &station,
                             const UUSSRayTracer::Phase phase,
                             const UUSSRayTracer::Region region,
                             ULocator::Corrections::Static &&staticCorrection,
                             ULocator::Corrections::SourceSpecific &&sssc,
                             std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<UUSSRayTracerImpl> (
        station,
        phase,
        region,
        std::move(staticCorrection),
        std::move(sssc),
        logger)
    )   
{
}

/*
UUSSRayTracer::UUSSRayTracer(const Station &station,
                             const UUSSRayTracer::Phase phase,
                             const UUSSRayTracer::Region region,
                             ULocator::Corrections::SourceSpecific &&sssc,
                             std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<UUSSRayTracerImpl> (
        station,
        phase,
        region,
        std::move(ULocator::Corrections::Static {}),
        std::move(sssc),
        logger)
    )
{
}
*/

std::vector<double> UUSSRayTracer::getInterfaces(
    const UUSSRayTracer::Region region)
{
    if (region == UUSSRayTracer::Region::Utah)
    {
        return std::vector<double> {-4500,  40,  15600, 26500, 40500};
    }
    else
    {
        return std::vector<double> {-4500, -1000,  2000,  5000,  8000,
                                    12000, 16000, 21000, 50000};
    }
}

std::vector<double> UUSSRayTracer::getPVelocities(
    const UUSSRayTracer::Region region)
{
    if (region == UUSSRayTracer::Region::Utah)
    {
        return std::vector<double> {3664, 5965, 6566, 7100, 7856};
    }
    else
    {
        return std::vector<double> {2713, 2787, 5207, 5565, 5812,
                                    6150, 6290, 6611, 7990};
    }
}

std::vector<double> UUSSRayTracer::getSVelocities(
    const UUSSRayTracer::Region region)
{
    if (region == UUSSRayTracer::Region::Utah)
    {
        return std::vector<double> {2004, 3480, 3858, 4168, 5006};
    }
    else
    {
        return std::vector<double> {1745, 2316, 3066, 3375, 3529,
                                    3650, 3713, 4014, 5054};
    }
}



double UUSSRayTracer::evaluate(
    const double t0, const double x, const double y, const double z,
    double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    if (pImpl->mRayTracer == nullptr)
    {
        throw std::runtime_error("UUSS ray tracer not initialized");
    }
    return pImpl->mRayTracer->evaluate(t0, x, y, z,
                                       dtdt0, dtdx, dtdy, dtdz,
                                       applyCorrection);
}

UUSSRayTracer::~UUSSRayTracer() = default;
