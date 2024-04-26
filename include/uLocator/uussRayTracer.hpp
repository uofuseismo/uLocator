#ifndef ULOCATOR_UUSS_RAY_TRACER_HPP
#define ULOCATOR_UUSS_RAY_TRACER_HPP
#include <memory>
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
/// @class UUSSRayTracer "uussRayTracer.hpp" "uLocator/uussRayTracer.hpp"
/// @brief This defines a 1D layer-cake ray-tracer that has hard-wired constants
///        that are appropriate for the UUSS monitoring regions.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class UUSSRayTracer : public ITravelTimeCalculator
{
public:
    /// @brief Defines the monitoring region.
    enum class Region
    {
        Utah, /*!< The Utah region. */
        YNP   /*!< Yellowstone National Park. */
    };
    /// @brief The phase type.
    enum class Phase
    {
        P,    /*!< Predicts first-arriving compressional waves. */
        S     /*!< Predicts first-arriving shear waves. */
    };
public:
    /// @brief Constructor.
    /// @param[in] station  Holds the station information - i.e., the name
    ///                     and geographic position.
    /// @param[in] phase    The phase - i.e., P or S.
    /// @param[in] region   The region - i.e., Utah or Yellowstone.
    /// @param[in] logger   The logging facility.
    UUSSRayTracer(const Station &station,
                  Phase phase,
                  Region region,
                  std::shared_ptr<UMPS::Logging::ILog> logger = nullptr);
    /// @brief Constructor.
    /// @param[in] station  Holds the station information - i.e., the name
    ///                     and geographic position.
    /// @param[in] phase    The phase - i.e., P or S.
    /// @param[in] region   The region - i.e., Utah or Yellowstone.
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
    UUSSRayTracer(const Station &station,
                  const Phase phase,
                  const Region region,
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
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    [[nodiscard]] double computeDistance(double x, double y) const final;
    /// @result The source-receiver distance in meters.
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    /// @param[in] z       The z position of the source in meters.
    [[nodiscard]] double computeDistance(double x, double y, double z) const final;
    /// @result The layer interfaces in the specified 1D model.
    [[nodiscard]] static std::vector<double> getInterfaces(Region region);
    /// @result The P velocity (m/s) in each layer of the specified 1D model.
    [[nodiscard]] static std::vector<double> getPVelocities(Region region);
    /// @result The S velocity (m/s) in each layer of the specified 1D model.
    [[nodiscard]] static std::vector<double> getSVelocities(Region region); 
    /// @brief Destructor. 
    ~UUSSRayTracer() override; 
    UUSSRayTracer() = delete;

    UUSSRayTracer(const UUSSRayTracer &) = delete;
    UUSSRayTracer(UUSSRayTracer &&) noexcept = delete;
    UUSSRayTracer& operator=(const UUSSRayTracer &) = delete;
    UUSSRayTracer& operator=(UUSSRayTracer &&) noexcept = delete;
private:
    class UUSSRayTracerImpl;
    std::unique_ptr<UUSSRayTracerImpl> pImpl;
};
}
#endif
