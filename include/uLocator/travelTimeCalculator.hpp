#ifndef ULOCATOR_TRAVEL_TIME_CALCULATOR_HPP
#define ULOCATOR_TRAVEL_TIME_CALCULATOR_HPP
#include <string>
namespace ULocator
{
 namespace Position
 {
  class WGS84;
 }
 class Station;
}
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
    /// @result The source-to-receiver travel time in seconds.
/*
    [[nodiscard]] virtual double evaluate(double utmX, double utmY, double depth, bool applyCorrection = true) const = 0;
    [[nodiscard]] virtual double evaluate(double utmX, double utmY, double depth,
                                          double *dtdx, double *dtdy, double *dtdz,
                                          bool applyCorrection = true) const = 0;
*/

    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] virtual double evaluate(const Position::WGS84 &epicenter, double depth, bool applyCorrection = true) const = 0;
    /// @param[out] dtdx   The derivative of the source-to-receiver travel time
    ///                    with respect to x.  This has units s/m.
    /// @param[out] dtdy   The derivative of the source-to-receiver travel time
    ///                    with respect to y.  This has units s/m.
    /// @param[out] dtdz   The derivative of the source-to-receiver travel time
    ///                    with respect to z.  This has units s/m.
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] virtual double evaluate(const Position::WGS84 &epicenter, double depth,
                                          double *dtdx, double *dtdy, double *dtdz,
                                          bool applyCorrection = true) const = 0;
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] virtual double operator()(const Position::WGS84 &epicenter, double depth, bool applyCorrection = true) const
    {
        return evaluate(epicenter, depth);
    }
};
}
#endif
