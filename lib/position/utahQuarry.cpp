#include <string>
#include "uLocator/position/utahQuarry.hpp"
#include "uLocator/position/utah.hpp"

using namespace ULocator::Position;

class UtahQuarry::UtahQuarryImpl
{
public:
    std::string mName;
};

/// Constructor
UtahQuarry::UtahQuarry() :
    IGeographicPoint(Utah {}),
    pImpl(std::make_unique<UtahQuarryImpl> ())
{
}

/// Copy constructor
UtahQuarry::UtahQuarry(const UtahQuarry &quarry) :
    IGeographicPoint(Utah {})
{
    *this = quarry;
}

/// Move constructor
UtahQuarry::UtahQuarry(UtahQuarry &&quarry) noexcept :
    IGeographicPoint(Utah {})
{
    *this = std::move(quarry);
}

/// Destructor
UtahQuarry::~UtahQuarry() = default;

/// Copy assignment
UtahQuarry& UtahQuarry::operator=(const UtahQuarry &quarry)
{
    if (&quarry == this){return *this;}
    pImpl = std::make_unique<UtahQuarryImpl> (*quarry.pImpl);
    if (quarry.haveGeographicCoordinates())
    {
        auto [latitude, longitude] = quarry.getGeographicCoordinates();
        setGeographicCoordinates(latitude, longitude);
    }
    return *this;
}

/// Move assignment
UtahQuarry& UtahQuarry::operator=(UtahQuarry &&quarry) noexcept
{
    if (&quarry == this){return *this;}
    if (quarry.haveGeographicCoordinates())
    {   
        auto [latitude, longitude] = quarry.getGeographicCoordinates();
        setGeographicCoordinates(latitude, longitude);
    }   
    pImpl = std::move(quarry.pImpl);
    return *this;
}

/// Name
void UtahQuarry::setName(const std::string &name)
{
    pImpl->mName = name;
}

std::string UtahQuarry::getName() const
{
    return pImpl->mName;
}

/// Clone operator
std::unique_ptr<IGeographicPoint> UtahQuarry::clone() const
{
    std::unique_ptr<IGeographicPoint> point
        = std::make_unique<UtahQuarry> (*this);
    return point;
}
