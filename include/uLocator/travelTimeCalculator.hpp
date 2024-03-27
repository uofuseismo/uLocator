#ifndef ULOCATOR_TRAVEL_TIME_CALCULATOR_HPP
#define ULOCATOR_TRAVEL_TIME_CALCULATOR_HPP
#include <string>
namespace ULocator
{
/// @class ITravelTimeCalculator "travelTimeCalculator.hpp" "uLocator/travelTimeCalculator.hpp"
/// @brief Defines the minimum functionality of a travel time calculator.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class ITravelTimeCalculator
{
public:
    /// @brief Destructor.
    virtual ~ITravelTimeCalculator() = default;
    /// @param[in] t0  The origin time in seconds.
    /// @param[in] x   The x position of the source in meters.
    /// @param[in] y   The y position of the source in meters.
    /// @param[in] z   The z position of the source in meters.
    /// @result The theoretical arrival time of the phase at the station
    ///         in seconds.
    [[nodiscard]] virtual double evaluate(double t0, double x, double y, double z,
                                          bool applyCorrection = true) const;
    /// @param[in] x   The x position of the source in meters.
    /// @param[in] y   The y position of the source in meters.
    /// @param[in] z   The z position of the source in meters.
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] virtual double evaluate(double x, double y, double z,
                                          bool applyCorrection = true) const;
    /// @param[in] t0      The origin time in seconds.
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    /// @param[in] z       The z position of the source in meters.
    /// @param[out] dtdt0  The derivative of the origin time.  This is unitless
    ///                    and will likely be unity.  If this is a nullptr then
    ///                    it will not be accessed.
    /// @param[out] dtdx   The derivative of the source-to-receiver travel time
    ///                    with respect to x.  This has units s/m.  If this is a
    ///                    nullptr then it will not be accessed.
    /// @param[out] dtdy   The derivative of the source-to-receiver travel time
    ///                    with respect to y.  This has units s/m.  If this is a
    ///                    nullptr then it will not be accessed.
    /// @param[out] dtdz   The derivative of the source-to-receiver travel time
    ///                    with respect to z.  This has units s/m.  If this is a
    ///                    nullptr then it will not be accessed.
    /// @result The theoretical arrival time of the phase in seconds.
    [[nodiscard]] virtual double evaluate(double t0, double x, double y, double z,
                                          double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
                                          bool applyCorrection = true) const = 0; 
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    /// @param[in] z       The z position of the source in meters.
    /// @param[out] dtdx   The derivative of the source-to-receiver travel time
    ///                    with respect to x.  This has units s/m.
    /// @param[out] dtdy   The derivative of the source-to-receiver travel time
    ///                    with respect to y.  This has units s/m.
    /// @param[out] dtdz   The derivative of the source-to-receiver travel time
    ///                    with respect to z.  This has units s/m.
    /// @param[in] applyCorrection  True indicates the static and 
    ///                             source-specific corrections should be
    ///                             applied; assuming they were set.
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] virtual double evaluate(double x, double y, double z,
                                          double *dtdx, double *dtdy, double *dtdz,
                                          bool applyCorrection = true) const;
    /// @param[in] t0      The origin time in seconds.
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    /// @param[in] z       The z position of the source in meters.
    /// @param[in] applyCorrection  True indicates the static and 
    ///                             source-specific corrections should be
    ///                             applied; assuming they were set.
    /// @result The theoretical arrival time of the phase at the station in
    ///         seconds.
    [[nodiscard]] virtual double operator()(double t0, double x, double y, double z, bool applyCorrection = true) const;
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    /// @param[in] z       The z position of the source in meters.
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] virtual double operator()(double x, double y, double z, bool applyCorrection = true) const;
    /// @result The source-receiver epicentral distance in meters.
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    [[nodiscard]] virtual double computeDistance(double x, double y) const = 0;
    /// @result The source-receiver distance in meters.
    /// @param[in] x       The x position of the source in meters.
    /// @param[in] y       The y position of the source in meters.
    /// @param[in] z       The z position of the source in meters.
    [[nodiscard]] virtual double computeDistance(double x, double y, double z) const = 0;

};
}
#endif
