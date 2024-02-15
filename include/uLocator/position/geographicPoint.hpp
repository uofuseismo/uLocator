#ifndef ULOCATOR_POSITION_GEOGRAPHIC_POINT_HPP
#define ULOCATOR_POSITION_GEOGRAPHIC_POINT_HPP
#include <memory>
#include <utility>
#include <uLocator/position/geographicRegion.hpp>
namespace ULocator::Position
{
/// @class IGeographicPoint "geographicPoint.hpp" "uLocator/postion/geographicPoint.hpp"
/// @brief A geographic point is a latitude/longitude point in a
//         geographic region.  
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class IGeographicPoint
{
public:
    /// @brief Sets the geographic position.
    explicit IGeographicPoint(const ULocator::Position::IGeographicRegion &region);
    /// @brief Destructor.
    virtual ~IGeographicPoint();
    /// @brief Sets the WGS84 geographic coordiantes.
    /// @param[in] latitude   The latitude of the point in degrees.
    /// @param[in] longitude  The longitude of the point in degrees.
    /// @throws std::invalid_argument if the latitude is not the range [-90,90].
    virtual void setGeographicCoordinates(double latitude, double longitude);
    /// @result True indicates the geographic coordaintes were set.
    [[nodiscard]] virtual bool haveGeographicCoordinates() const noexcept;
    /// @brief From the local WGS84 geographic coordinates this returns
    ///        the local Cartesian (x,y) coordinates.
    /// @result result.first is the x position and result.second is the y
    ///         position; both in meters.
    [[nodiscard]] virtual std::pair<double, double> getLocalCoordinates() const;
    /// @brief From the local Cartesian (x,y) this returns the WGS84
    ///        latitude and longitude.
    /// @result result.first is the latitude and result.second is the longitude;
    ///         both in degrees.
    [[nodiscard]] virtual std::pair<double, double> getGeographicCoordinates() const;
    /// @result The geographic region of this point.
    [[nodiscard]] std::unique_ptr<IGeographicRegion> getGeographicRegion() const;
    /// @result A copy of the class.
    [[nodiscard]] virtual std::unique_ptr<IGeographicPoint> clone() const = 0;
private:
    class IGeographicPointImpl;
    std::unique_ptr<IGeographicPointImpl> pImpl; 
};
}
#endif
