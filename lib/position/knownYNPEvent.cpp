#include <vector>
#include <memory>
#include "uLocator/position/knownYNPEvent.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/position/knownLocalLocation.hpp"

using namespace ULocator::Position;

class KnownYNPEvent::KnownYNPEventImpl
{
public:
    double mX{0};
    double mY{0};
    double mDepth{0};
    bool mInitialized{false};
};

KnownYNPEvent::KnownYNPEvent(const double latitude,
                               const double longitude,
                               const double depth) :
    pImpl(std::make_unique<KnownYNPEventImpl> ())
{
    ULocator::Position::YNPRegion region;
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

KnownYNPEvent::KnownYNPEvent(const KnownYNPEvent &event)
{
    *this = event;
}

KnownYNPEvent::KnownYNPEvent(KnownYNPEvent &&event) noexcept
{
    *this = std::move(event);
}

double KnownYNPEvent::x() const
{
    if (!pImpl->mInitialized){throw std::runtime_error("Class not constructed");}
    return pImpl->mX;
}

double KnownYNPEvent::y() const
{
    if (!pImpl->mInitialized){throw std::runtime_error("Class not constructed");}
    return pImpl->mY;
}

double KnownYNPEvent::z() const
{
    if (!pImpl->mInitialized){throw std::runtime_error("Class not constructed");}
    return pImpl->mDepth;
}

std::unique_ptr<IKnownLocalLocation> KnownYNPEvent::clone() const
{
    std::unique_ptr<IKnownLocalLocation> result
        = std::make_unique<KnownYNPEvent> (*this);
    return result;
}

KnownYNPEvent& KnownYNPEvent::operator=(const KnownYNPEvent &event)
{
    if (&event == this){return *this;}
    pImpl = std::make_unique<KnownYNPEventImpl> (*event.pImpl);
    return *this;
}

KnownYNPEvent& KnownYNPEvent::operator=(KnownYNPEvent &&event) noexcept
{
    if (&event == this){return *this;}
    pImpl = std::move(event.pImpl);
    return *this;
}
 
KnownYNPEvent::~KnownYNPEvent() = default;

std::vector<ULocator::Position::KnownYNPEvent> 
    ULocator::Position::getKnownYNPEvents()
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
