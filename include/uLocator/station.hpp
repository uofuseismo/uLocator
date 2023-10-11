#ifndef ULOCATOR_STATION_HPP
#define ULOCATOR_STATION_HPP
#include <memory>
namespace ULocator::Position
{
 class WGS84;
}
namespace ULocator
{
class Station
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Station();
    /// @brief Copy constructor.
    /// @param[in] station  The station class from which to initialize
    ///                     this class.
    Station(const Station &station);
    /// @brief Move constructor.
    /// @param[in,out] station  The station class from which to initialize
    ///                         this class.  On exit, station's behavior
    ///                         is undefined.
    Station(Station &&station) noexcept;
    /// @}

    /// @brief Sets the network code.
    /// @param[in] network  The station's network code - e.g., UU.
    void setNetwork(const std::string &network);
    /// @result The network code.
    /// @throws std::runtime_error if \c haveNetwork() is false.
    [[nodiscard]] std::string getNetwork() const;
    /// @result True indicates the network code was set.
    [[nodiscard]] bool haveNetwork() const noexcept;

    /// @brief Sets the station name.
    /// @param[in] name   The station name - e.g., FORK.
    /// @throws std::invalid_argument if the name is empty.
    void setName(const std::string &name);
    /// @result The station name. 
    /// @throws std::runtime_error if \c haveName() is false.
    [[nodiscard]] std::string getName() const;
    /// @result True indicates the station name was set.
    [[nodiscard]] bool haveName() const noexcept;
    /// @result The hash for this station.
    /// @throws std::runtime_error if \c haveName() or \c haveNetwork()
    ///         is false.
    [[nodiscard]] size_t getHash() const;

    /// @brief Sets the latitude and longitude of the station.
    /// @param[in] position  The geographic position of the station.
    void setGeographicPosition(const Position::WGS84 &position);
    /// @result The geographic position of the station.
    [[nodiscard]] Position::WGS84 getGeographicPosition() const;
    /// @result The geographic position of the station.
    /// @note This exists for performance reasons.  You should prefer
    ///       \c getGeographicPosition().
    [[nodiscard]] const Position::WGS84& getGeographicPositionReference() const;
    /// @result True indicates the geographic position was set. 
    [[nodiscard]] bool haveGeographicPosition() const noexcept; 

    /// @brief Sets the station's elevation.
    /// @param[in] elevation  The station's elevation in meters with
    ///                       respect to sea-level.
    /// @throws std::invalid_argument if this is the center of the Earth or
    ///         exceeds the elevation of Mt. Everest.
    void setElevation(double elevation);
    /// @result The station's elevation.
    [[nodiscard]] double getElevation() const;
    /// @result True indicates the station's elevation was set.
    [[nodiscard]] bool haveElevation() const noexcept;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] station  The station class to copy to this.
    /// @result A deep copy of the station class.
    Station& operator=(const Station &station);
    /// @brief Move assignment.
    /// @param[in,out] station  The station class whose memory will be moved
    ///                         to this.  On exit, station's behavior is
    ///                         undefined.
    /// @result The memory from station moved to this.
    Station& operator=(Station &&station) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Station();
    /// @}
private:
    class StationImpl;
    std::unique_ptr<StationImpl> pImpl;
};
}
#endif
