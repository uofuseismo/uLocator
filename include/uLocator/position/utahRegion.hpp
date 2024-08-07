#ifndef ULOCATOR_POSITION_UTAH_REGION_HPP
#define ULOCATOR_POSITION_UTAH_REGION_HPP
#include <memory>
#include <uLocator/position/geographicRegion.hpp>
namespace ULocator::Position
{
/// @brief Defines the Utah geographic region for optimization.
///        This region is approximately -114.75 to -108.25 degrees
///        longitude and approximately 36.25 to 43 degrees latitude.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class UtahRegion : public IGeographicRegion
{
public:
    /// @name Constructors.
    /// @{

    /// @brief Constructor.
    UtahRegion();
    /// @brief Copy constructor.
    /// @param[in] utah   The Utah region from which to initialize this class.
    UtahRegion(const UtahRegion &utah);
    /// @brief Move constructor.
    /// @param[in,out] utah  The Utah region from which to initialize this
    ///                      class.  On exit, utah's behavior is undefined.
    UtahRegion(UtahRegion &&utah) noexcept;
    /// @}

    /// @name Destructors.
    /// @{

    /// @brief Destructor.
    ~UtahRegion();
    /// @}

    /// @result Converts the given latitude, longitude pair to a local x, y coordinate in meters.
    [[nodiscard]] std::pair<double, double> geographicToLocalCoordinates(double latitude, double longitude) const final;
    /// @result Converts the given local (x, y) coordinates in metres to latitude and longitude in degrees.
    [[nodiscard]] std::pair<double, double> localToGeographicCoordinates(double x, double y) const final;

    /// @result The minimum and maximum x extent of the Utah search region in 
    ///         meters in the local coordinate system.
    [[nodiscard]] std::pair<double, double> getExtentInX() const noexcept final;
    /// @result The minimum and maximum y extent of the Utah search region in
    ///         meters in the local coordinate system. 
    [[nodiscard]] std::pair<double, double> getExtentInY() const noexcept final;

    /// @name Operators
    /// @{

    /// @brief Copy assignment. 
    /// @param[in] utah  The Utah region to copy to this.
    /// @result A deep copy of the Utah boundaries.
    UtahRegion& operator=(const UtahRegion &utah);
    /// @brief Move assignment.
    /// @param[in,out] utah  The Utah region whose memory will be moved to this.
    ///                      On exit, utah's behavior is undefined.
    /// @result The memory from the Utah boundaries moved to this.
    UtahRegion& operator=(UtahRegion &&utah) noexcept;
    /// @result A copy of the Utah region.
    [[nodiscard]] std::unique_ptr<IGeographicRegion> clone() const final; 
    /// @}
private:
    class UtahRegionImpl;
    std::unique_ptr<UtahRegionImpl> pImpl;
};
}
#endif
