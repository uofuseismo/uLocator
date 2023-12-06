#ifndef ULOCATOR_POSITION_WGS84_HPP
#define ULOCATOR_POSITION_WGS84_HPP
#include <memory>
namespace ULocator::Position
{
/// @class WGS84 "wgs84.hpp" "uLocator/position/wgs84.hpp"
/// @brief Defines a geographic position in the WGS84 projection.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class WGS84
{
public:
    /// @name Constructor
    /// @{

    /// @brief Constructor.
    /// @param[in] utmZone  The UTM zone.  By default this will be gleaned
    ///                     from the latitude and longitude. 
    WGS84();//const int utmZone =-1);
    /// @brief Copy constructor.
    /// @param[in] position  The position class from which to initialize
    ///                      this class.
    WGS84(const WGS84 &position);
    /// @brief Move constructor.
    /// @param[in,out] position  The position class from which to initialize
    ///                          this class.  On exit, position's behavior
    ///                          is undefined.
    WGS84(WGS84 &&position) noexcept;
    /// @brief Constructs a WGS84 position from the given latitude and
    ///        longitude.
    /// @param[in] latitude       The latitude in degrees.
    /// @param[in] longitude      The longitude in degrees.
    /// @param[in] alternateZone  The alternate UTM zone.  By default this
    ///                           will be set from the latitude and longitude.
    /// @sa \c setLatitudeAndLongitude(), \c setAlternateZone()
    WGS84(double latitude, double longitude, int alternateZone =-1);
    /// @brief Constructs a WGS84 position from an easting and northing.
    /// @param[in] zone      The UTM zone number.
    /// @param[in] isNorth   True indicates this is in the northern hemisphere.
    /// @param[in] easting   The UTM easting in meters.
    /// @param[in] northing  The UTM northing in meters.
    WGS84(int zone, bool isNorth, double easting, double northing);
    /// @brief

    /// @}

    /// @name Latitude and Longitude
    /// @{

    /// @brief Sets the latitude and longitude.
    /// @param[in] latitude  The latitude in degrees where positive is north.
    /// @param[in] longitude  The longitude in degrees where positive is east.
    /// @throws std::runtime_error if the latitude is not in range [-90,90].
    void setLatitudeAndLongitude(double latitude, double longitude);
    /// @result The latitude in degrees.
    /// @throws std::runtime_error if \c havePosition() is false.
    [[nodiscard]] double getLatitude() const;

    /// @result The longitude in degrees.  This will be in the range [-180,180).
    /// @throws std::runtime_error if \c havePosition() is false.
    [[nodiscard]] double getLongitude() const;

    /// @brief Defines an alternate UTM zone.
    /// @param[in] zone  The alternate UTM zone.  -1 indicates to use
    ///                  the default from the latitude and longitude.
    ///                  0 is for UPS.  Otherwise this will be in the
    ///                  range of 1 to 60.
    /// @throws std::invalid_argument if this is out of bounds. 
    void setAlternateUTMZone(int zone);
    /// @result The alternate UTM zone.  By default this is -1 indicating it
    ///         is undefined.
    [[nodiscard]] int getAlternateUTMZone() const noexcept;
    /// @}

    /// @name UTM Easting and Northing
    /// @{

    /// @brief Sets the position as a northing and easting.
    /// @param[in] zone      The UTM zone number.  This must be in the
    ///                      range [0,60] where 0 indicates UPS (useful for
    ///                      high latitudes).
    /// @param[in] isNorth   True indicates this is in the northern hemisphere.
    /// @param[in] easting   The easting in meters.
    /// @param[in] northing  The northing in meters.
    void setEastingAndNorthing(int zone, bool isNorth, double easting, double northing);
    /// @result The UTM easting in meters.
    /// @throws std::runtime_error if \c havePosition() is false.
    [[nodiscard]] double getEasting() const;
    /// @result The UTM northing in meters.
    /// @throws std::runtime_error if \c havePosition() is false. 
    [[nodiscard]] double getNorthing() const; 
    /// @result The UTM zone.
    /// @throws std::runtime_error if \c havePosition() is false.
    [[nodiscard]] int getUTMZone() const;
    /// @result True indicates this is in the northern hemisphere.
    /// @throws std::runtime_error if \c havePosition() is false.
    [[nodiscard]] bool isNorth() const;
    /// @}

    /// @result True indicates the position was set.
    [[nodiscard]] bool havePosition() const noexcept;

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] position  The position to copy to this.
    /// @result A deep copy of the input position.
    WGS84& operator=(const WGS84 &position);
    /// @brief Move assignment operator.
    /// @param[in,out] position  The position whose memory will be moved to
    ///                          this.  On exit, position's behavior is
    ///                          undefined.
    /// @result The memory from position moved to this.
    WGS84& operator=(WGS84 &&position) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~WGS84();
    /// @}
private:
    class WGS84Impl;
    std::unique_ptr<WGS84Impl> pImpl;
};

/// @param[in] source   The source position.
/// @param[in] station  The station position.
/// @param[out] greatCircleDistance  The great circle distance between the
///                                  sourcce and station in degrees.
/// @param[out] distance             The source-receiver distance in meters.
/// @param[out] azimuth              The source-to-station azimuth measured
///                                  positive east of north in degrees.
/// @param[out] backAzimuth          The station-to-source azimuth measured
///                                  positive east of north in degrees.
/// @throws std::invalid_argument if source.havePosition() is
///         station.havePosition() is false.
void computeDistanceAzimuth(const WGS84 &source,
                            const WGS84 &station,
                            double *greatCircleDistance,
                            double *distance,
                            double *azimuth,
                            double *backAzimuth);
 
}
#endif
