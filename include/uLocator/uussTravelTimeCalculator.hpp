#ifndef ULOCATOR_UUSS_TRAVEL_TIME_CALCULATOR_HPP
#define ULOCATOR_UUSS_TRAVEL_TIME_CALCULATOR_HPP
#include <memory>
#include <uLocator/travelTimeCalculator.hpp>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator
{
 class Station;
}
namespace ULocator
{
/// @class UUSSTravelTimeCalculator "uussTravelTimeCalculator.hpp" "uLocator/uussTravelTimeCalculator.hpp"
/// @brief The travel time calculator which works with stations and UUSS
///        station corrections.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class UUSSTravelTimeCalculator : public ITravelTimeCalculator
{
public:
    /// @name Construtors
    /// @{

    UUSSTravelTimeCalculator();
    explicit UUSSTravelTimeCalculator(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @}

    /// @brief Loads the travel time field.
    void load(const std::string &fileName,
              const Station &station, const std::string &phase, const std::string &region = "utah");
    /// @brief Loads the static correction.
    void loadStaticCorrection(const std::string &fileName);
    /// @brief Loads the source-specific station corrections.
    void loadSourceSpecificStationCorrections(const std::string &fileName);

    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] double evaluate(const Position::WGS84 &epicenter, double depth,
                                  bool evaluateCorrection = true) const override;
    /// @param[out] dtdx   The derivative of the source-to-receiver travel time
    ///                    with respect to x.  This has units s/m.
    /// @param[out] dtdy   The derivative of the source-to-receiver travel time
    ///                    with respect to y.  This has units s/m.
    /// @param[out] dtdz   The derivative of the source-to-receiver travel time
    ///                    with respect to z.  This has units s/m.
    /// @result The source-to-receiver travel time in seconds.
    [[nodiscard]] double evaluate(const Position::WGS84 &epicenter, double depth,
                                  double *dtdx, double *dtdy, double *dtdz,
                                  bool evaluateCorrection = true) const override;

    /// @name Destructors
    /// @{
 
    /// Release memory and reset class.
    void clear() noexcept;
    /// @brief Destructor.
    ~UUSSTravelTimeCalculator() override;
    /// @}

    UUSSTravelTimeCalculator(const UUSSTravelTimeCalculator &) = delete;
    UUSSTravelTimeCalculator& operator=(const UUSSTravelTimeCalculator &) = delete;
private:
    class UUSSTravelTimeCalculatorImpl;
    std::unique_ptr<UUSSTravelTimeCalculatorImpl> pImpl;
};
}
#endif
