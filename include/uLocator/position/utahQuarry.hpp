#ifndef ULOCATOR_POSITION_UTAH_QUARRY_HPP
#define ULOCATOR_POSITION_UTAH_QUARRY_HPP
#include <memory>
#include <string>
#include <uLocator/position/geographicPoint.hpp>
namespace ULocator::Position
{
/// @class UtahQuarry "utahQuarry.hpp" "uLocator/position/utahQuarry.hpp"
/// @brief This is a quarry in the Utah geographic region.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class UtahQuarry : public IGeographicPoint
{
public:
    /// @name Constructors
    /// @{

    /// @brief Default constructor.
    UtahQuarry();
    /// @brief Copy constructor.
    /// @param[in] quarry  The quarry from which to initialize this class.
    UtahQuarry(const UtahQuarry &quarry);
    /// @brief Move constructor.
    /// @param[in,out] quarry  The quarry from which to initialize this class.
    ///                        On exit, quarry's behavior is undefined.
    UtahQuarry(UtahQuarry &&quarry) noexcept;
    /// @brief Fully constructs the quarry.
    /// @param[in] latitude   The latitude of the quarry in degrees.
    /// @param[in] longitude  The longitude of the quarry in degrees.
    /// @param[in] name       The name of the quarry.
    /// @throws std::invalid_argument if the latitude is out of range.
    UtahQuarry(const double latitude,
               const double longitude,
               const std::string &name);
    /// @brief Fully constructs the quarry.
    /// @param[in] latitude   The latitude of the quarry in degrees.
    /// @param[in] longitude  The longitude of the quarry in degrees.
    /// @param[in] elevation  The elevation above sea-level in meters of
    ///                       the quarry.
    /// @param[in] name       The name of the quarry.
    /// @throws std::invalid_argument if the latitude is out of range.
    UtahQuarry(double latitude,
               double longitude,
               double elevation,
               const std::string &name);
    /// @}

    /// @name Quarry Name
    /// @{

    /// @param[in] name  The name of the quarry.
    void setName(const std::string &name);
    /// @result The name of the quarry.
    [[nodiscard]] std::string getName() const;

    /// @brief Set elevation.
    /// @param[in] elevation  The elevation  above sea-level in meters of
    ///                       the quarry. 
    void setElevation(double elevation);
    /// @result The elevation, in meters, above sea level of the quarry.
    [[nodiscard]] double getElevation() const noexcept;
    /// @}

    /// @result A copy of this class.
    [[nodiscard]] std::unique_ptr<IGeographicPoint> clone() const override;

    /// @name Destructors
    /// {

    /// @brief Destructor.
    ~UtahQuarry() override;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] quarry   The quarry to copy to this.
    /// @result A deep copy of the input quarry.
    UtahQuarry& operator=(const UtahQuarry &quarry);
    /// @brief Move assignment operator.
    /// @result[in,out] quarry  The quarry whose memory will be moved to this.
    /// @result The memory from quarry moved to this.
    UtahQuarry& operator=(UtahQuarry &&quarry) noexcept;
    /// @}
private:
    class UtahQuarryImpl;
    std::unique_ptr<UtahQuarryImpl> pImpl;
};
}
#endif
