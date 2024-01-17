#include <vector>
#include <mutex>
#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include <umps/logging/log.hpp>
#include <umps/logging/standardOut.hpp>
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
#include "uLocator/uussRayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/geographicRegion.hpp"

using namespace ULocator;

class UUSSRayTracer::UUSSRayTracerImpl
{
public:
    UUSSRayTracerImpl(const Station &station,
                      const UUSSRayTracer::Phase phase,
                      const UUSSRayTracer::Region region,
                      std::shared_ptr<UMPS::Logging::ILog> &logger) :
        mLogger(logger),
        mStation(station)
    {
        if (!mStation.haveNetwork())
        {
            throw std::invalid_argument("Station's network not set");
        }
        if (!mStation.haveName())
        {
            throw std::invalid_argument("Station name not set");
        }
        if (!mStation.haveGeographicPosition())
        {
            throw std::invalid_argument("Station geographic position not set");
        }
        if (!mStation.haveElevation())
        {
            throw std::invalid_argument("Station elevation not set");
        }
        mStationDepth = -mStation.getElevation();
        mStationX = mStation.getLocalCoordinates().first;
        mStationY = mStation.getLocalCoordinates().second;
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
        if (region == UUSSRayTracer::Region::Utah)
        {
            std::vector<double> interfaces{-4500,  40,  15600, 26500, 40500};
            interfaces[0] = std::min(mStationDepth - 1.e-8, interfaces[0]);
            mMinimumDepth = interfaces.at(0);
            if (phase == UUSSRayTracer::Phase::P) 
            {
                mLogger->debug("Initializing Utah P layer solver...");
                std::vector<double> velocities{3685, 5972, 6553, 7106, 7672};
                mLayerSolver.setVelocityModel(interfaces, velocities);
                mPhase = "P";
            }
            else if (phase == UUSSRayTracer::Phase::S)
            {
                mLogger->debug("Initializing Utah S layer solver...");
                std::vector<double> velocities{2169, 3443, 3653, 4486, 5022};
                mLayerSolver.setVelocityModel(interfaces, velocities);
                mPhase = "S";
            }
#ifndef NDEBUG
            else
            {
               assert(false);
            }
#endif
        }
        else if (region == UUSSRayTracer::Region::YNP)
        {
            std::vector<double> interfaces{-4500, -1000,  2000,  5000,  8000,
                                           12000, 16000, 21000, 50000};
            interfaces[0] = std::min(mStationDepth - 1.e-8, interfaces[0]);
            mMinimumDepth = interfaces.at(0);
            if (phase == UUSSRayTracer::Phase::P)
            {
                mLogger->debug("Initializing YNP P layer solver...");
                std::vector<double> velocities{2512, 3398, 4689, 5456, 5674,
                                               6250, 6398, 6574, 8200};
                mLayerSolver.setVelocityModel(interfaces, velocities);
                mPhase = "P";
            }
            else if (phase == UUSSRayTracer::Phase::S)
            {
                mLogger->debug("Initializing YNP S layer solver...");
                std::vector<double> velocities{1725, 2343, 3064, 3425, 3569,
                                               3690, 3705, 3975, 4950};
                mLayerSolver.setVelocityModel(interfaces, velocities);
                mPhase = "S";
            }
       }
#ifndef NDEBUG
       else
       {
            assert(false);
       }
       assert(mLayerSolver.haveVelocityModel());
#endif        

    } 
    mutable std::mutex mMutex;
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    Station mStation;
    mutable EikonalXX::Ray::LayerSolver mLayerSolver;
    std::string mPhase;
    double mStationX{0};
    double mStationY{0};
    double mStationDepth{0};
    double mMinimumDepth{0};
};

UUSSRayTracer::UUSSRayTracer(const Station &station,
                             const UUSSRayTracer::Phase phase,
                             const UUSSRayTracer::Region region,
                             std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<UUSSRayTracerImpl> (station, phase, region, logger))
{
}

double UUSSRayTracer::evaluate(
    const double t0, const double x, const double y, const double z,
    double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    double dx = pImpl->mStationX - x; // Ray points source to receiver
    double dy = pImpl->mStationY - y;
    double sourceDepthToUse = z;
    if (sourceDepthToUse < pImpl->mMinimumDepth)
    {
        auto warningMessage = "Event depth " + std::to_string(sourceDepthToUse)
                            + " being overriden to " 
                            + std::to_string(pImpl->mMinimumDepth);
        pImpl->mLogger->warn(warningMessage);
        sourceDepthToUse = pImpl->mMinimumDepth;
    }
    auto offset = std::hypot(dx, dy);

    // Evaluate solver - in case something parallel calls this we need to lock
    // access to our solver.  Since this is the only implementation in this
    // source this is all we need to do.
    double travelTime{0};
    double takeOffAngleRadians{0};
    double slowness{0};
    {
    std::lock_guard lock(pImpl->mMutex);
    try
    {
        pImpl->mLayerSolver.setStationOffsetAndDepth(offset,
                                                     pImpl->mStationDepth);
        pImpl->mLayerSolver.setSourceDepth(sourceDepthToUse);
        pImpl->mLayerSolver.solve();
    }
    catch (const std::exception &e)
    {
         auto errorMessage = "Failed to compute ray paths for "
                           + pImpl->mStation.getNetwork() + "." 
                           + pImpl->mStation.getName() + "." 
                           + pImpl->mPhase
                           + ".  Failed with: "
                           + std::string {e.what()}
                           +  ".  Source depth is: " 
                           + std::to_string(sourceDepthToUse);
         pImpl->mLogger->error(errorMessage);
         throw std::runtime_error(errorMessage);
    }
    const auto &rayPaths = pImpl->mLayerSolver.getRayPaths();
    if (rayPaths.empty())
    {
        throw std::runtime_error("Ray paths are empty!");
    }
    travelTime = rayPaths[0].getTravelTime();
    takeOffAngleRadians = rayPaths[0].getTakeOffAngle()*(M_PI/180);
    slowness = rayPaths[0].at(0).getSlowness();
    } // End mutex

    // Derivatives
    if (dtdt0 != nullptr){*dtdt0 = 1;}
    if (dtdx != nullptr || dtdy != nullptr || dtdz != nullptr)
    {
        // Decompose initial ray path into slowness in radial and z (+down)
        double dtdr = slowness*std::sin(takeOffAngleRadians);
        if (dtdz != nullptr){*dtdz = slowness*std::cos(takeOffAngleRadians);}
        // A few more considerations for going from radial to cartesian
        // coordinates:
        if (dtdx != nullptr || dtdy != nullptr)
        {
            // We have computed a travel time of the form T(r,z) so
            // {dT/dr}   = [   cos(theta)   sin(theta) 0 ] {dT/dx}
            // {dT/dphi} = [-r sin(theta) r cos(theta) 0 ] {dT/dy}
            // {dT/dz}   = [     0              0      1 ] {dT/dz}
            // The inverse relations are 
            // {dT/dx}   = [   cos(theta)  -sin(theta)/r 0 ] {dT/dr}
            // {dT/dz}   = [   sin(theta)   cos(theta)/r 0 ] {dT/dtheta}
            // {dT/dz}   = ]      0             0        1 ] {dT/dz}
            // Note, the layer tracer models are 1D so dT/dtheta = 0.  Thus:
            //   dT/dx = cos(theta) dT/dr
            //   dT/dy = sin(theta) dT/dr
            //   dT/dz = dT/dz
            auto thetaRadians = std::atan2(dy, dx);
            if (dtdx != nullptr){*dtdx = dtdr*std::cos(thetaRadians);}
            if (dtdy != nullptr){*dtdy = dtdr*std::sin(thetaRadians);}
        } 
    }
    if (applyCorrection)
    {
    }
    return t0 + travelTime; 
}

UUSSRayTracer::~UUSSRayTracer() = default;
