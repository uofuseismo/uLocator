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
#include "uLocator/rayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/corrections/static.hpp"
#include "uLocator/corrections/sourceSpecific.hpp"

using namespace ULocator;

class RayTracer::RayTracerImpl
{
public:
    RayTracerImpl(const Station &station,
                  const std::string &phase,
                  const std::vector<double> &interfaces,
                  const std::vector<double> &velocities,
                  std::shared_ptr<UMPS::Logging::ILog> &logger) :
        mLogger(logger),
        mStation(station),
        mPhase(phase)
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
        if (mPhase.empty()){throw std::invalid_argument("Phase is empty");}
        if (interfaces.empty())
        {
            throw std::invalid_argument("Interfaces is empty");
        }
        if (velocities.size() != interfaces.size())
        {
            throw std::invalid_argument(
                "Size of velocity != size of interfaces");
        }
        mStationDepth =-mStation.getElevation();
        mStationX = mStation.getLocalCoordinates().first;
        mStationY = mStation.getLocalCoordinates().second;
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
        setVelocityModel(interfaces, velocities);
    }
    void setVelocityModel(const std::vector<double> &interfacesIn,
                          const std::vector<double> &velocities)
    {
        auto interfaces = interfacesIn;
        interfaces.at(0) = std::min(mStationDepth - 1.e-8, interfaces.at(0));
        mMinimumDepth = interfaces.at(0);
        mLayerSolver.setVelocityModel(interfaces, velocities);
#ifndef NDEBUG
        assert(mLayerSolver.haveVelocityModel());
#endif
    } 
    RayTracerImpl(const Station &station,
                  const std::string &phase,
                  const std::vector<double> &interfaces,
                  const std::vector<double> &velocities,
                  ULocator::Corrections::Static &&staticCorrection,
                  ULocator::Corrections::SourceSpecific &&sourceSpecific,
                  std::shared_ptr<UMPS::Logging::ILog> &logger) :
        RayTracerImpl(station, phase, interfaces, velocities, logger)
    {
        if (staticCorrection.haveCorrection())
        {
            mStaticCorrection = std::move(staticCorrection);
        }
        if (sourceSpecific.haveModel())
        {
            mSourceSpecific = std::move(sourceSpecific);
        }
    }
    mutable std::mutex mMutex;
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    Corrections::Static mStaticCorrection;
    Corrections::SourceSpecific mSourceSpecific;
    Station mStation;
    mutable EikonalXX::Ray::LayerSolver mLayerSolver;
    std::string mPhase;
    double mStationX{0};
    double mStationY{0};
    double mStationDepth{0};
    double mMinimumDepth{0};
};

RayTracer::RayTracer(const Station &station,
                     const std::string &phase,
                     const std::vector<double> &interfaces,
                     const std::vector<double> &velocities,
                     std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<RayTracerImpl> (station, phase, 
                                           interfaces, velocities,
                                           logger))
{
}

RayTracer::RayTracer(const Station &station,
                     const std::string &phase,
                     const std::vector<double> &interfaces,
                     const std::vector<double> &velocities,
                     ULocator::Corrections::Static &&staticCorrection,
                     ULocator::Corrections::SourceSpecific &&sourceSpecific,
                     std::shared_ptr<UMPS::Logging::ILog> logger) :
    pImpl(std::make_unique<RayTracerImpl> (station, phase, 
                                           interfaces, velocities,
                                           std::move(staticCorrection),
                                           std::move(sourceSpecific),
                                           logger))
{
}

/// Source-receiver distance
double RayTracer::computeDistance(const double x, const double y) const
{
    if (!pImpl->mLayerSolver.haveVelocityModel())
    {   
        throw std::runtime_error("Solver not initialized");
    }   
    double dx = pImpl->mStationX - x; // Ray points source to receiver
    double dy = pImpl->mStationY - y;
    return std::hypot(dx, dy);
}

double RayTracer::computeDistance(
    const double x, const double y, const double z) const
{
    if (!pImpl->mLayerSolver.haveVelocityModel())
    {   
        throw std::runtime_error("Solver not initialized");
    }   
    double dx = pImpl->mStationX - x; // Ray points source to receiver
    double dy = pImpl->mStationY - y;
    double dz = pImpl->mStationDepth - z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/// Apply the solver
double RayTracer::evaluate(
    const double t0, const double x, const double y, const double z,
    double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    if (!pImpl->mLayerSolver.haveVelocityModel())
    {
        throw std::runtime_error("Solver not initialized");
    }
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
        constexpr bool doReflections{false};
        pImpl->mLayerSolver.solve(doReflections);
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
    travelTime = rayPaths.at(0).getTravelTime();
    takeOffAngleRadians = rayPaths.at(0).getTakeOffAngle()*(M_PI/180);
    slowness = rayPaths.at(0).at(0).getSlowness();
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
        constexpr double zeroTime{0};
        if (dtdx != nullptr or dtdy != nullptr or dtdz != nullptr) 
        {
            double dcdt, dcdx, dcdy, dcdz;
            if (pImpl->mStaticCorrection.haveCorrection())
            {
                auto staticCorrection
                    = pImpl->mStaticCorrection.evaluate(&dcdt,
                                                        &dcdx,
                                                        &dcdy,
                                                        &dcdz,
                                                        zeroTime);
                if (travelTime + staticCorrection >= 0)
                {
                    travelTime = travelTime + staticCorrection;
                    if (dtdt0 != nullptr){*dtdt0 = *dtdt0 + dcdt;}
                    if (dtdx  != nullptr){*dtdx  = *dtdx  + dcdx;}
                    if (dtdy  != nullptr){*dtdy  = *dtdy  + dcdy;}
                    if (dtdz  != nullptr){*dtdz  = *dtdz  + dcdz;}
                }
            }
            if (pImpl->mSourceSpecific.haveModel())
            {
                auto sourceSpecificCorrection
                    = pImpl->mSourceSpecific.evaluate(x, y, z,
                                                      &dcdt,
                                                      &dcdx,
                                                      &dcdy,
                                                      &dcdz,
                                                      zeroTime);
                if (travelTime + sourceSpecificCorrection >= 0)
                {
                    travelTime = travelTime + sourceSpecificCorrection;
                    if (dtdt0 != nullptr){*dtdt0 = *dtdt0 + dcdt;}
                    if (dtdx  != nullptr){*dtdx  = *dtdx  + dcdx;}
                    if (dtdy  != nullptr){*dtdy  = *dtdy  + dcdy;}
                    if (dtdz  != nullptr){*dtdz  = *dtdz  + dcdz;}
                }
            }
        }
        else
        {
            if (pImpl->mStaticCorrection.haveCorrection())
            {
                auto staticCorrection
                    = pImpl->mStaticCorrection.evaluate(zeroTime);
                if (travelTime + staticCorrection >= 0)
                {
                    travelTime = travelTime + staticCorrection;
                }
            }
            if (pImpl->mSourceSpecific.haveModel())
            {
                auto sourceSpecificCorrection
                    = pImpl->mSourceSpecific.evaluate(x, y, z, zeroTime);
                if (travelTime + sourceSpecificCorrection >= 0)
                {
                    travelTime = travelTime + sourceSpecificCorrection;
                }
            }
        }
    }
    return t0 + travelTime; 
}

RayTracer::~RayTracer() = default;
