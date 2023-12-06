#ifndef ULOCATOR_POSITION_GEOGRAPHIC_REGION_HPP
#define ULOCATOR_POSITION_GEOGRAPHIC_REGION_HPP
#include <utility>
namespace ULocator::Position
{
/// @class IGeographicRegion "geographicRegion.hpp" "uLocator/postion/geographicRegion.hpp"
/// @brief A geographic region defines the region that will be searched by
///        an optimization scheme.  In the local coordinates Cartesian
///        coordinates the region must be square - this is for the benefit
///        of the optimization scheme which is more easily defined in 
///        (x, y, z, t) for network seismology tasks.  Of course, near the poles
///        there will be over-sampling in longitude at the extreme latitudes.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class IGeographicRegion
{
public:
    /// @brief Destructor.
    virtual ~IGeographicRegion() = default;
    /// @brief From the local WGS84 geodetic coordinates this returns
    ///        the local Cartesian (x,y) coordinates.
    /// @param[in] latitude   The latitude in degrees.
    /// @param[in] longitude  The longitude in degrees.
    /// @result result.first is the x position and result.second is the y
    ///         position; both in meters.
    /// @throws std::invalid_argument if the latitude is not the range [-90,90].
    [[nodiscard]] virtual std::pair<double, double> geodeticToLocalCoordinates(double latitude, double longitude) const = 0;
    /// @brief From the local Cartesian (x,y) this returns the WGS84
    ///        latitude and longitude.
    /// @param[in] x  The x local coordinate in meters.
    /// @param[in] y  The y local coordiante in meters.
    /// @result result.first is the latitude and result.second is the longitude;
    ///         both in degrees.
    [[nodiscard]] virtual std::pair<double, double> localToGeodeticCoordinates(double x, double y) const = 0;
    /// @result result.first is the smallest x in the region and result.second
    ///         is the largest x; both in meters. 
    [[nodiscard]] virtual std::pair<double, double> getExtentInX() const noexcept = 0;
    /// @result result.first is the smallest y in the region and result.second
    ///         is the largest y; both in meters.
    [[nodiscard]] virtual std::pair<double, double> getExtentInY() const noexcept = 0;
};
}
#endif
