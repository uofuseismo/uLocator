#ifndef ULOCATOR_POSITION_KNOWN_YNP_EVENT_HPP
#define ULOCATOR_POSITION_KNOWN_YNP_EVENT_HPP
#include <vector>
#include <memory>
#include <uLocator/position/knownLocalLocation.hpp>
namespace ULocator::Position
{
/// @class KnownYNPEvent "knownYNPEvent.hpp" "uLocator/position/knownYNPEvent.hpp"
/// @brief A known YNP provides the local coordinates of a place to search
///        during an optimization.  For example, you may define such a place
///        as a swarm centroid.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class KnownYNPEvent : public IKnownLocalLocation
{
public:
    /// @brief Constructor
    /// @param[in] latitude   The event latitude in degrees.
    /// @param[in] longitude  The event longitude in degrees.
    /// @param[in] depth      The event depth in meters.
    /// @throws std::invalid_argument if the latitude is not in the
    ///         range [-90,90] or the depth is less than -8500 meters
    ///         or greater than 800000 meters.
    KnownYNPEvent(const double latitude,
                   const double longitude,
                   const double depth);
    /// @brief Copy constructor.
    /// @param[in] event  The event from which to initialize this class.
    KnownYNPEvent(const KnownYNPEvent &event);
    /// @brief Move constructor.
    /// @param[in,out] event  The event from which to initialize this class.
    ///                       On exit, event's behavior is undefined.
    KnownYNPEvent(KnownYNPEvent &&event) noexcept;
    /// @result The local x position in meters of the event.
    [[nodiscard]] double x() const override;
    /// @result The local y position in meters of the event.
    [[nodiscard]] double y() const override;
    /// @result The local z position in meters of the event.
    [[nodiscard]] double z() const override;
    /// @result A copy of the class.
    [[nodiscard]] std::unique_ptr<IKnownLocalLocation> clone() const override;
    /// @brief Copy assignment operator.
    /// @param[in] event  The event to copy to this.
    /// @result A deep copy of the input event.
    KnownYNPEvent& operator=(const KnownYNPEvent &event);
    /// @brief Move assingment operator.
    /// @param[in,out] event  The event whose memory will be moved to this.
    ///                       On exit, event's behavior is undefined.
    /// @result The memory from event moved to this.
    KnownYNPEvent& operator=(KnownYNPEvent &&event) noexcept;
    
    /// @brief destructor.
    ~KnownYNPEvent() override;
    KnownYNPEvent() = delete;
private:
    class KnownYNPEventImpl;
    std::unique_ptr<KnownYNPEventImpl> pImpl;
};
/// @result The current default event search locations in YNP.
std::vector<ULocator::Position::KnownYNPEvent> getKnownYNPEvents();
}
#endif
