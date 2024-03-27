#ifndef ULOCATOR_RAY_TRACER_HPP
#define ULOCATOR_RAY_TRACER_HPP
#include <memory>
#include <vector>
#include <string>
#include <uLocator/travelTimeCalculator.hpp>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator
{
 class Station;
 namespace Corrections
 {
  class Static;
  class SourceSpecific;
 }
}
namespace ULocator
{
/// @class RayTracer "rayTracer.hpp" "uLocator/rayTracer.hpp"
/// @brief This defines a 1D layer-cake ray-tracer.
/// @bug This cannot handle velocity inversions.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class RayTracer : public ITravelTimeCalculator
{
public:
    /// @brief Constructor.
    /// @param[in] station     Holds the station name and location.
    /// @param[in] phase       Defines the seismic phase.
    /// @param[in] interfaces  The interfaces of the velocity model w.r.t.
    ///                        sea-level in meters.  This increases +z down
    ///                        into the earth - i.e., the first interface
    ///                        is the free-surface and interfaces.back()
    ///                        is the top of the underlying half-space.
    /// @param[in] velocities  The velocities in each layer in m/s.  Here
    ///                        velocities.front is the velocity in the top-most
    ///                        layer and velocities.back is the velocity in the
    ///                        underlying half-space.
    /// @param[in] logger   The logging facility.
    RayTracer(const Station &station,
              const std::string &phase,
              const std::vector<double> &interfaces,
              const std::vector<double> &velocities,
              std::shared_ptr<UMPS::Logging::ILog> logger = nullptr);
    /// @brief Constructor.
    /// @param[in] station     Holds the station name and location.
    /// @param[in] phase       Defines the seismic phase.
    /// @param[in] interfaces  The interfaces of the velocity model w.r.t.
    ///                        sea-level in meters.  This increases +z down
    ///                        into the earth - i.e., the first interface
    ///                        is the free-surface and interfaces.back()
    ///                        is the top of the underlying half-space.
    /// @param[in] velocities  The velocities in each layer in m/s.  Here
    ///                        velocities.front is the velocity in the top-most
    ///                        layer and velocities.back is the velocity in the
    ///                        underlying half-space.
    /// @param[in,out] staticCorrection  Optionally, the static correction.
    ///                                  If stationCorrection.haveCorrection()
    ///                                  is true then the memory from
    ///                                  staticCorrection will be moved to this
    ///                                  and its behavior is undefined on exit. 
    /// @param[in,out] sourceSpecific  Optionally, the source-specific 
    ///                                correction.  If sourceSpecific.haveModel()
    ///                                is true then the memory from
    ///                                sourceSpecific will be moved to this
    ///                                and its behavior is undefined on exit. 
    /// @param[in] logger   The logging facility.
    RayTracer(const Station &station,
              const std::string &phase,
              const std::vector<double> &interfaces,
              const std::vector<double> &velocities,
              ULocator::Corrections::Static &&staticCorrection,
              ULocator::Corrections::SourceSpecific &&sourceSpecific,
              std::shared_ptr<UMPS::Logging::ILog> logger = nullptr);
    /// @param[in] t0   The origin time in seconds.
    /// @param[in] x    The x position of the source in meters.
    /// @param[in] y    The y position of the source in meters.
    /// @param[in] z    The z position of the source in meters.
    /// @param[in,out] dtdt0  If this is not NULL then this is the derivative
    ///                       of the travel time w.r.t. to the origin time.
    ///                       This is unitless.
    /// @param[in,out] dtdx   If this i snot NULL then this is the derivative
    ///                       of the travel time w.r.t. to the x source
    ///                       position.  This has units of seconds per meter. 
    /// @param[in,out] dtdy   If this i snot NULL then this is the derivative
    ///                       of the travel time w.r.t. to the y source
    ///                       position.  This has units of seconds per meter. 
    /// @param[in,out] dtdz   If this i snot NULL then this is the derivative
    ///                       of the travel time w.r.t. to the z source
    ///                       position.  This has units of seconds per meter. 
    /// @param[in] applyCorrection  True indicates the static and/or source
    ///                             specific corrections should be applied
    ///                             provided they were set.
    [[nodiscard]] double evaluate(double t0, double x, double y, double z,
                                  double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
                                  bool applyCorrection) const final;
    /// @result The source-receiver epicentral distance in meters.
    [[nodiscard]] double computeDistance(double x, double y) const final;
    /// @result The distance between the source and receiver in meters.
    [[nodiscard]] double computeDistance(double x, double y, double z) const final;
    /// @brief Destructor. 
    ~RayTracer() override; 
    RayTracer() = delete;

    RayTracer(const RayTracer &) = delete;
    RayTracer(RayTracer &&) noexcept = delete;
    RayTracer& operator=(const RayTracer &) = delete;
    RayTracer& operator=(RayTracer &&) noexcept = delete;
private:
    class RayTracerImpl;
    std::unique_ptr<RayTracerImpl> pImpl;
};
}
#endif
