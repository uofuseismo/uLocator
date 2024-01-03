#ifndef ULOCATOR_ORIGIN_HPP
#define ULOCATOR_ORIGIN_HPP
#include <memory>
#include <vector>
namespace ULocator
{
 class Arrival;
 namespace Position
 {
  class WGS84;
 }
}
namespace ULocator
{
/// @class Origin "origin.hpp" "uLocator/origin.hpp"
/// @brief Defines an event origin.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Origin
{
public:
    enum class EventType
    {
        Unknown = 0,
        Earthquake,
        QuarryBlast
    };
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Origin();
    /// @brief Copy constructor.
    /// @param[in] origin  The origin from which to initialize this class.
    Origin(const Origin &origin);
    /// @brief Move constructor.
    /// @param[in,out] origin  The origin from which to initialize this class.
    ///                        On exit, origin's behavior is undefined.
    Origin(Origin &&origin) noexcept;
    /// @}

    /// @name Origin Time
    /// @{

    /// @brief Sets the origin time.
    /// @param[in] time  The origin time (UTC) measured in seconds since the
    ///                  epoch.
    void setTime(double time) noexcept;
    /// @result The origin time. 
    [[nodiscard]] double getTime() const;
    /// @result True indicates the origin time was set.
    [[nodiscard]] bool haveTime() const noexcept;
    /// @}

    /// @name Epicenter
    /// @{

    /// @brief Sets the event epicenter.
    /// @throws std::invalid_argument if the latitude or longitude is not set.
    void setEpicenter(const ULocator::Position::WGS84 &epicenter);
    /// @result The epicenter.
    /// @throws std::runtime_error if \c haveEpicenter() is false.
    [[nodiscard]] ULocator::Position::WGS84 getEpicenter() const;
    /// @result True indicates the epicenter was set.
    [[nodiscard]] bool haveEpicenter() const noexcept;
    /// @}

    /// @name Depth
    /// @{

    /// @brief Sets the event depth.
    /// @param[in] depth  The event depth in meters measured on the WGS84
    ///                   reference ellipsoid.
    /// @param[in] isFixed  True indicates this depth was not a free parameter
    ///                     in the estimation and was set.
    /// @throws std::invalid_argument if this exceeds the 6,400,000 
    ///         which is approximately the radius of the earth.
    void setDepth(double depth, bool isFixed = false);
    /// @result The event depth in meters.  This uses the WGS84 reference
    ///         ellipsoid so 0 is sea level.
    [[nodiscard]] double getDepth() const;
    /// @result True indicates the event depth was set.
    [[nodiscard]] bool haveDepth() const noexcept;
    /// @result True indicates the depth is fixed.
    [[nodiscard]] bool isFixedDepth() const noexcept;
    /// @}

    /// @name Event Identifier
    /// @{

    /// @brief Sets the origin identifier.
    void setIdentifier(int64_t identifier) noexcept;
    /// @result The origin identifier.
    [[nodiscard]] int64_t getIdentifier() const noexcept;
    /// @}

    /// @name Arrivals
    /// @{

    /// @brief Sets the arrivals afixed to this origin.
    /// @param[in] arrivals  The arrivals.  This must have the station
    ///                      name, phase, and arrival time.
    void setArrivals(const std::vector<Arrival> &arrivals);
    /// @result The arrivals.
    [[nodiscard]] std::vector<Arrival> getArrivals() const;
    /// @result A reference to the underlying arrivals.
    [[nodiscard]] const std::vector<Arrival> &getArrivalsReference() const;
    /// @}

    /// @name Event Type
    /// @{

    /// @brief Sets the event type.
    void setEventType(EventType eventType) noexcept;
    /// @result The event type.  By default this is unknown.
    [[nodiscard]] EventType getEventType() const noexcept;
    /// @}

    /// @name Derivative Products
    /// @{

    /// @result For the given arrivals and epicenter this is the largest 
    ///         azimuthal gap in station coverage in degrees.
    /// @throws std::runtime_error if \c haveAzimuthalGap() is false.
    [[nodiscard]] double getAzimuthalGap() const;
    /// @result True indicates the azimuthal gap was computed.
    [[nodiscard]] bool haveAzimuthalGap() const noexcept;
    /// @result The smallest source-station distance in meters.
    /// @throws std::runtime_error if \c haveNearestStationDistance() is false.
    [[nodiscard]] double getNearestStationDistance() const;
    /// @result True indicates
    [[nodiscard]] bool haveNearestStationDistance() const noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] origin  The origin to copy to this.
    /// @result A deep copy of the input origin.
    Origin& operator=(const Origin &origin);
    /// @brief Move assignment.
    /// @param[in,out] origin  The origin whose memory will be moved to this.
    ///                        On exit, origin's behavior is undefined.
    /// @result The memory from origin moved to this.
    Origin& operator=(Origin &&origin) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Origin();
    /// @}
private:
    class OriginImpl;
    std::unique_ptr<OriginImpl> pImpl;
};
}
#endif
