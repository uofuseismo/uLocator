#ifndef ULOCATOR_ARRIVAL_HPP
#define ULOCATOR_ARRIVAL_HPP
#include <string>
#include <memory>
namespace ULocator
{
 class Station;
}
namespace ULocator
{
/// @class Arrival "arrival.hpp" "uLocator/arrival.hpp"
/// @brief Defines an arrival time.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Arrival
{
public:
    /// @brief A convenience utility for defining P and S phase types.
    enum class PhaseType
    {
        P = 0,
        S = 1
    };
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Arrival();
    /// @brief Copy constructor.
    /// @param[in] arrival  The arrival class from which to copy this pick.
    Arrival(const Arrival &arrival);
    /// @brief Move constructor.
    /// @param[in,out] arrival  The arrival class from which to initialize this
    ///                         class.  On exit, arrival's behavior is
    ///                         undefined.
    Arrival(Arrival &&arrival) noexcept;
    /// @}

    /// @name Required Parameters
    /// @{

    /// @brief Sets the observation's arrival time.
    /// @param[in] time  The arrival time (UTC) measured in seconds since
    ///                  the epoch (Jan 1 1970) of the phase. 
    void setTime(double time) noexcept;
    /// @result The arrival time.
    [[nodiscard]] double getTime() const;
    /// @result True indicates the arrival time was set.
    [[nodiscard]] bool haveTime() const noexcept;

    /// @brief Sets the standard error for the arrival time.
    /// @param[in] standardError  Assuming a Gaussian model for the pick
    ///                           time uncertainty, this defines the +/- one
    ///                           standard deviation in arrival time in seconds.
    /// @throws std::invalid_argument if this is not positive.
    void setStandardError(double standardError);
    /// @result The standard error deviation  
    [[nodiscard]] double getStandardError() const;
    /// @result True indicates the standard error was set.
    [[nodiscard]] bool haveStandardError() const noexcept;

    /// @brief Sets the phase arrival type.
    /// @param[in] phase  The phase type.
    void setPhase(const PhaseType phase)  noexcept;
    /// @brief Sets the phase arrival type.
    /// @param[in] phase  The phase type - e.g., P or S.
    /// @throw std::invalid_argument if the phase type is empty.
    void setPhase(const std::string &phase);
    /// @result The phsae of the arrival.
    [[nodiscard]] std::string getPhase() const;
    /// @result True indicates the phase type was set.
    [[nodiscard]] bool havePhase() const noexcept;

    /// @brief Sets the station information on which the arrival was observed.
    /// @param[in] station  The station information.
    /// @throws std::invalid_argument if the network code or station name was
    ///         not set.
    void setStation(const Station &station);
    /// @result The station information.
    [[nodiscard]] Station getStation() const;
    /// @result A reference to the station information.
    /// @note This exists for optimization reasons and \c getStation()
    ///       should be preferred. 
    [[nodiscard]] const Station& getStationReference() const;
    /// @result True indicates the station information was set.
    [[nodiscard]] bool haveStation() const noexcept;
    /// @}

    /// @name Optional Parameters
    /// @{

    /// @brief Sets the travel time residual which is observed - predicted.
    /// @param[in] residual  The travel time residual in seconds.
    void setResidual(double residual) noexcept;
    /// @result The travel time residual (observed - predicted) in seconds.
    [[nodiscard]] double getResidual() const;
    /// @result True indicates the travel time residual was set.
    [[nodiscard]] bool haveResidual() const noexcept;

    /// @brief Sets the unique arrival identifier.
    /// @param[in] identifier  The arrival identifier.
    void setIdentifier(const int64_t identifier) noexcept;
    /// @result The arrival identifier.
    [[nodiscard]] int64_t getIdentifier() const noexcept;

    /// @brief Sets the source-receiver distance.
    /// @param[in] distance  The source-receiver distance in meters.
    /// @throws std::invalid_argument if distance is negative.
    void setDistance(double distance);
    /// @result The source-receiver distance in meters.
    [[nodiscard]] double getDistance() const;
    /// @result True indicates the source-receiver distance was set.
    [[nodiscard]] bool haveDistance() const noexcept;

    /// @brief Sets the source-to-receiver azimuth.
    /// @param[in] azimuth   The source-to-receiver azimuth in degrees measured
    ///                      positive-east of north.
    /// @throws std::invalid_argument if azimuth is not in the range [0,360).
    void setAzimuth(double azimuth);
    /// @result The source-to-receiver azimuth in degrees.
    [[nodiscard]] double getAzimuth() const;
    /// @result True indicates the source-to-receiver azimuth was not set.
    [[nodiscard]] bool haveAzimuth() const noexcept;

    /// @brief Sets the receiver-to-source azimuth.
    /// @param[in] backAzimuth  The receiver-to-source azimuth in degrees
    ///                         measured positive-east of north.
    /// @throws std::invalid_argument if backAzimuthis not in the range [0,360).
    void setBackAzimuth(double backAzimuth);
    /// @result The receiver-to-source azimuth in degrees.
    [[nodiscard]] double getBackAzimuth() const;
    /// @result True indicates the receiver-to-source azimuth was not set.
    [[nodiscard]] bool haveBackAzimuth() const noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @param[in] arrival  The arrival to copy to this.
    /// @result A deep copy of the input arrival.
    Arrival& operator=(const Arrival &arrival);
    /// @param[in,out] arrival  The arrival whose memory will be moved to this.
    ///                         On exit, arrival's behavior is undefined. 
    /// @result The memory from arrival moved to this.
    Arrival& operator=(Arrival &&arrival) noexcept;
    /// @}

    /// @name Destructors 
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Arrival();
    /// @} 
private:
    class ArrivalImpl;
    std::unique_ptr<ArrivalImpl> pImpl;
};
}
#endif
