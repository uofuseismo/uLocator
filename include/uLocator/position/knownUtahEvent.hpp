#ifndef ULOCATOR_POSITION_KNOWN_UTAH_EVENT_HPP
#define ULOCATOR_POSITION_KNOWN_UTAH_EVENT_HPP
#include <vector>
#include <memory>
#include <uLocator/position/knownLocalLocation.hpp>
namespace ULocator::Position
{
/// @class KnownUtahEvent "knownUtahEvent.hpp" "uLocator/position/knownUtahEvent.hpp"
/// @brief A known Utah provides the local coordinates of a place to search
///        during an optimization.  For example, you may define such a place
///        as an aftershock centroid.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class KnownUtahEvent : public IKnownLocalLocation
{
public:
    /// @brief Constructor
    /// @param[in] latitude   The event latitude in degrees.
    /// @param[in] longitude  The event longitude in degrees.
    /// @param[in] depth      The event depth in meters.
    /// @throws std::invalid_argument if the latitude is not in the
    ///         range [-90,90] or the depth is less than -8500 meters
    ///         or greater than 800000 meters.
    KnownUtahEvent(const double latitude,
                   const double longitude,
                   const double depth);
    /// @brief Copy constructor.
    /// @param[in] event  The event from which to initialize this class.
    KnownUtahEvent(const KnownUtahEvent &event);
    /// @brief Move constructor.
    /// @param[in,out] event  The event from which to initialize this class.
    ///                       On exit, event's behavior is undefined.
    KnownUtahEvent(KnownUtahEvent &&event) noexcept;
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
    KnownUtahEvent& operator=(const KnownUtahEvent &event);
    /// @brief Move assingment operator.
    /// @param[in,out] event  The event whose memory will be moved to this.
    ///                       On exit, event's behavior is undefined.
    /// @result The memory from event moved to this.
    KnownUtahEvent& operator=(KnownUtahEvent &&event) noexcept;
    
    /// @brief destructor.
    ~KnownUtahEvent() override;
    KnownUtahEvent() = delete;
private:
    class KnownUtahEventImpl;
    std::unique_ptr<KnownUtahEventImpl> pImpl;
};
/// @result The current default event search locations in Utah.
std::vector<ULocator::Position::KnownUtahEvent> getKnownUtahEvents();
/*
{
    std::vector<ULocator::Position::KnownUtahEvent> events{
          KnownUtahEvent (37.461412, -114.035346, 5692.475248),
          KnownUtahEvent {36.624736, -112.928499, 10421.710526},
          KnownUtahEvent {37.873043, -110.507263, 6413.653846},
          KnownUtahEvent {37.604445, -113.467926, 5037.192982},
          KnownUtahEvent {39.425405, -111.956333, 6537.510917},
          KnownUtahEvent {37.109262, -112.853806, 16126.457399},
          KnownUtahEvent {38.622589, -112.154460, 6378.744186},
          KnownUtahEvent {36.960758, -112.308144, 16690.987654},
          KnownUtahEvent {39.104179, -109.079286, 4593.200000},
          KnownUtahEvent {41.885767, -112.510553, 4003.215078},
          KnownUtahEvent {37.733097, -112.461298, 6198.433180},
          KnownUtahEvent {39.934295, -111.697597, 5262.568807},
          KnownUtahEvent {41.679496, -109.739289, 4795.172414},
          KnownUtahEvent {40.851306, -111.577731, 9097.089947},
          KnownUtahEvent {39.453065, -111.245483, 3270.021598},
          KnownUtahEvent {41.452097, -112.820454, 5088.750000},
          KnownUtahEvent {38.131472, -112.324373, 5459.079498},
          KnownUtahEvent {42.352765, -111.392682, 5288.333333},
          KnownUtahEvent {38.285265, -113.027094, 4811.232227},
          KnownUtahEvent {42.041059, -112.021268, 4437.246377},
          KnownUtahEvent {36.833859, -113.974749, 7014.936709},
          KnownUtahEvent {36.974799, -113.544430, 9758.940397},
          KnownUtahEvent {38.927739, -111.416057, 3340.490998},
          KnownUtahEvent {41.359725, -111.718352, 8720.059880},
          KnownUtahEvent {41.753682, -111.630500, 8369.726776},
          KnownUtahEvent {40.427553, -111.249713, 7976.808511},
          KnownUtahEvent {39.712275, -110.753325, -1744.260870},
          KnownUtahEvent {39.066326, -110.938887, 6542.880000},
          KnownUtahEvent {40.288048, -113.434374, 7364.166667},
          KnownUtahEvent {38.618886, -112.584403, 3942.283951},
          KnownUtahEvent {38.200110, -111.281470, 7854.000000},
          KnownUtahEvent {41.077816, -110.645927, 15837.500000},
          KnownUtahEvent {37.216025, -109.159245, 6680.000000},
          KnownUtahEvent {40.481416, -112.003325, 5898.465608},
          KnownUtahEvent {40.747322, -112.053724, 8168.579161},
          KnownUtahEvent {37.492067, -113.041498, 7294.243697},
          KnownUtahEvent {38.472293, -112.831202, 1329.542396},
          KnownUtahEvent {37.797463, -113.113250, 5969.433198},
          KnownUtahEvent {40.181620, -108.852801, 6799.090909},
          KnownUtahEvent {39.423903, -110.309202, -1709.524398},
          KnownUtahEvent {38.255689, -108.899164, 3095.952381},
          KnownUtahEvent {38.901065, -114.055206, 3992.500000},
          KnownUtahEvent {38.190056, -112.656680, 5713.505155},
          KnownUtahEvent {39.706717, -111.223339, 1222.222222},
          KnownUtahEvent {41.537563, -112.265330, 3749.936709},
          KnownUtahEvent {40.971500, -109.883167, 58110.000000}
    };
    return events;
}
*/

/*
std::vector<ULocator::Position::KnownYNPEvent> getKnownYNPEvents()
{
    std::vector<KnownYNPEvent> events{
          KnownYNPEvent {44.042996, -110.297501, 8447.500000},
          KnownYNPEvent {44.615943, -111.001421, 8805.620253},
          KnownYNPEvent {44.598103, -110.740705, 5737.287615},
          KnownYNPEvent {44.781297, -110.478683, 5067.831325},
          KnownYNPEvent {44.775547, -110.793679, 5308.191230},
          KnownYNPEvent {44.752074, -111.139571, 9971.768868},
          KnownYNPEvent {44.424962, -110.990331, 7659.076305},
          KnownYNPEvent {44.252810, -110.719965, 5639.319149},
          KnownYNPEvent {44.377787, -110.728424, 3730.296610},
          KnownYNPEvent {44.523714, -111.096710, 14262.535211},
          KnownYNPEvent {43.798107, -110.957131, 8371.851852},
          KnownYNPEvent {44.447925, -110.361887, 4800.767386},
          KnownYNPEvent {44.834802, -111.416688, 10724.375000},
          KnownYNPEvent {44.329277, -110.498721, 4987.519380},
          KnownYNPEvent {44.615620, -110.409933, 3963.893805},
          KnownYNPEvent {44.788535, -111.025461, 7755.848075},
          KnownYNPEvent {44.450988, -110.568044, 3137.515152},
          KnownYNPEvent {44.065834, -110.637226, 9194.404762},
          KnownYNPEvent {44.675678, -110.022057, 10811.473684},
          KnownYNPEvent {44.757390, -110.909975, 6440.319716},
          KnownYNPEvent {44.736921, -110.997092, 6580.154440},
          KnownYNPEvent {44.755156, -110.683070, 3646.256158}
    };
    return events;
}
*/

}
#endif
