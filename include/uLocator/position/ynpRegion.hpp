#ifndef ULOCATOR_POSITION_YNP_REGION_HPP
#define ULOCATOR_POSITION_YNP_REGION_HPP
#include <memory>
#include <uLocator/position/geographicRegion.hpp>
namespace ULocator::Position
{
/// @brief Defines the YNP geographic region for optimization.
///        This region is approximately -114.75 to -108.25 degrees
///        longitude and approximately 36.25 to 43 degrees latitude.
/// @copyright Ben Baker (University of YNP) distributed under the MIT license.
class YNPRegion : public IGeographicRegion
{
public:
    /// @name Constructors.
    /// @{

    /// @brief Constructor.
    YNPRegion();
    /// @brief Copy constructor.
    /// @param[in] ynp  The YNP region from which to initialize this class.
    YNPRegion(const YNPRegion &ynp);
    /// @brief Move constructor.
    /// @param[in,out] ynp  The YNP region from which to initialize this class.
    ///                     On exit, YNP's behavior is undefined.
    YNPRegion(YNPRegion &&ynp) noexcept;
    /// @}

    /// @name Destructors.
    /// @{

    /// @brief Destructor.
    ~YNPRegion();
    /// @}

    /// @result Converts the given latitude, longitude pair to a local x, y coordinate in meters.
    [[nodiscard]] std::pair<double, double> geographicToLocalCoordinates(double latitude, double longitude) const final;
    /// @result Converts the given local (x, y) coordinates in metres to latitude and longitude in degrees.
    [[nodiscard]] std::pair<double, double> localToGeographicCoordinates(double x, double y) const final;

    /// @result The minimum and maximum x extent of the YNP search region in 
    ///         meters in the local coordinate system.
    [[nodiscard]] std::pair<double, double> getExtentInX() const noexcept final;
    /// @result The minimum and maximum y extent of the YNP search region in
    ///         meters in the local coordinate system. 
    [[nodiscard]] std::pair<double, double> getExtentInY() const noexcept final;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] ynp  The YNP region to copy to this.
    /// @result A deep copy of the YNP boundaries.
    YNPRegion& operator=(const YNPRegion &ynp);
    /// @brief Move assignment.
    /// @param[in,out] ynp  The YNP region whose memory will be moved to this.
    ///                     On exit, YNP's behavior is undefined.
    /// @result The memory from the YNP boundaries moved to this.
    YNPRegion& operator=(YNPRegion &&ynp) noexcept;
    /// @result A copy of the YNP region.
    [[nodiscard]] std::unique_ptr<IGeographicRegion> clone() const final; 
    /// @}
private:
    class YNPRegionImpl;
    std::unique_ptr<YNPRegionImpl> pImpl;
};
}
#endif
