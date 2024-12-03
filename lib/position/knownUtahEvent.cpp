#include <vector>
#include <memory>
#include "uLocator/position/knownUtahEvent.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/position/knownLocalLocation.hpp"

using namespace ULocator::Position;

class KnownUtahEvent::KnownUtahEventImpl
{
public:
    double mX{0};
    double mY{0};
    double mDepth{0};
    bool mInitialized{false};
};

KnownUtahEvent::KnownUtahEvent(const double latitude,
                               const double longitude,
                               const double depth) :
    pImpl(std::make_unique<KnownUtahEventImpl> ())
{
    ULocator::Position::UtahRegion region;
    auto localCoordinates
        = region.geographicToLocalCoordinates(latitude, longitude);
    if (depth < -8500 || depth > 800000)
    {
        throw std::invalid_argument(
             "Depth not between [-8,500 and 800,000] m");
    }
    pImpl->mX = localCoordinates.first;
    pImpl->mY = localCoordinates.second;
    pImpl->mDepth = depth;
    pImpl->mInitialized = true;
}

KnownUtahEvent::KnownUtahEvent(const KnownUtahEvent &event)
{
    *this = event;
}

KnownUtahEvent::KnownUtahEvent(KnownUtahEvent &&event) noexcept
{
    *this = std::move(event);
}

double KnownUtahEvent::x() const
{
    if (!pImpl->mInitialized){throw std::runtime_error("Class not constructed");}
    return pImpl->mX;
}

double KnownUtahEvent::y() const
{
    if (!pImpl->mInitialized){throw std::runtime_error("Class not constructed");}
    return pImpl->mY;
}

double KnownUtahEvent::z() const
{
    if (!pImpl->mInitialized){throw std::runtime_error("Class not constructed");}
    return pImpl->mDepth;
}

std::unique_ptr<IKnownLocalLocation> KnownUtahEvent::clone() const
{
    std::unique_ptr<IKnownLocalLocation> result
        = std::make_unique<KnownUtahEvent> (*this);
    return result;
}

KnownUtahEvent& KnownUtahEvent::operator=(const KnownUtahEvent &event)
{
    if (&event == this){return *this;}
    pImpl = std::make_unique<KnownUtahEventImpl> (*event.pImpl);
    return *this;
}

KnownUtahEvent& KnownUtahEvent::operator=(KnownUtahEvent &&event) noexcept
{
    if (&event == this){return *this;}
    pImpl = std::move(event.pImpl);
    return *this;
}
 
KnownUtahEvent::~KnownUtahEvent() = default;

std::vector<ULocator::Position::KnownUtahEvent> 
    ULocator::Position::getKnownUtahEvents()
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
          KnownUtahEvent {40.971500, -109.883167, 58110.000000},
          KnownUtahEvent {38.500000, -112.890000, 3200.000000} // FERVO
    };
    return events;
}

