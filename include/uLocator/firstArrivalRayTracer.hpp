#ifndef ULOCATOR_FIRST_ARRIVAL_RAY_TRACER_HPP
#define ULOCATOR_FIRST_ARRIVAL_RAY_TRACER_HPP
#include <memory>
//#include <umps/logging/log.hpp>
#include "uLocator/travelTimeCalculator.hpp"
namespace ULocator
{
class FirstArrivalRayTracer : public ITravelTimeCalculator
{
public:
    FirstArrivalRayTracer();
    explicit FirstArrivalRayTracer(std::shared_ptr<UMPS::Logging::ILog> &logger);
    FirstArrivalRayTracer(FirstArrivalRayTracer &&rayTracer) noexcept;

    void initialize(const Station &station,
                    const std::string &phase,
                    const std::vector<double> &interfaces,
                    const std::vector<double> &velocities);

    void initializeUtahP(const Station &station, const bool useAlternateModel = true);
    void initializeUtahS(const Station &station, const bool useAlternateModel = true);

    void initializeYellowstoneP(const Station &station);
    void initializeYellowstoneS(const Station &station);
    [[nodiscard]] bool isInitialized() const noexcept;


/*
    [[nodiscard]] double evaluate(double utmX, double utmY, double depth, bool applyCorrection = true) const override;
    [[nodiscard]] double evaluate(double utmX, double utmY, double depth, 
                                  double *dtdx, double *dtdy, double *dtdz,
                                  bool applyCorrection = true) const override;
*/

    /// @result The first-arrival travel time from the given event location.
    /// @param[in] epicenter        The event's epicenter.
    /// @param[in] depth            The event's depth in meters.
    /// @param[in] applyCorrection  If true then this will apply travel time
    ///                             corrections.
    /// @result The source-to-receiver travel time in seconds.
    /// @throws std::invalid_argument if the epicenter's latitude and longitude
    ///         is not set or the depth makes this an air quake.
    [[nodiscard]] double evaluate(const Position::WGS84 &epicenter,
                                  double depth, bool applyCorrection = true) const override;
    /// @param[in] epicenter        The event's epicenter.
    /// @param[in] depth            The event's depth in meters.
    /// @param[in] applyCorrection  If true then this will apply travel time
    ///                             corrections.
    /// @param[out] dtdx   The derivative of the source-to-receiver travel time
    ///                    with respect to x.  This has units s/m.
    /// @param[out] dtdy   The derivative of the source-to-receiver travel time
    ///                    with respect to y.  This has units s/m.
    /// @param[out] dtdz   The derivative of the source-to-receiver travel time
    ///                    with respect to z.  This has units s/m.
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] double evaluate(const Position::WGS84 &epicenter, double depth,
                                  double *dtdx, double *dtdy, double *dtdz,
                                  bool applyCorrection = true) const override;

    /// @brief Resets the class and releases memory.
    void clear() noexcept; 
    /// @brief Destructor.
    ~FirstArrivalRayTracer() override;
    /// @result The memory from rayTracer moved to this.
    FirstArrivalRayTracer& operator=(FirstArrivalRayTracer &&rayTracer) noexcept;
private:
    class FirstArrivalRayTracerImpl;
    std::unique_ptr<FirstArrivalRayTracerImpl> pImpl;
};
}
#endif
