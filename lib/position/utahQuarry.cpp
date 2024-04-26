#include <string>
#include "uLocator/position/utahQuarry.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "utahQuarries.hpp"

using namespace ULocator::Position;

class UtahQuarry::UtahQuarryImpl
{
public:
    std::string mName;
    double mElevation{0};
};

/// Constructor
UtahQuarry::UtahQuarry() :
    IGeographicPoint(UtahRegion {}),
    pImpl(std::make_unique<UtahQuarryImpl> ())
{
}

/// Copy constructor
UtahQuarry::UtahQuarry(const UtahQuarry &quarry) :
    IGeographicPoint(UtahRegion {})
{
    *this = quarry;
}

/// Move constructor
UtahQuarry::UtahQuarry(UtahQuarry &&quarry) noexcept :
    IGeographicPoint(UtahRegion {})
{
    *this = std::move(quarry);
}

/// Constructor
UtahQuarry::UtahQuarry(const double latitude,
                       const double longitude,
                       const std::string &name) :
    UtahQuarry()
{
    setName(name);
    setGeographicCoordinates(latitude, longitude);
}

/// Constructor
UtahQuarry::UtahQuarry(const double latitude,
                       const double longitude,
                       const double elevation,
                       const std::string &name) :
    UtahQuarry()
{
    setName(name);
    setGeographicCoordinates(latitude, longitude);
    setElevation(elevation);
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

/// Elevation
void UtahQuarry::setElevation(const double elevation)
{
    if (elevation >= 8850)
    {
        throw std::invalid_argument("Elevation must be less than 8,850 m");
    }
    if (elevation <= -10000)
    {
        throw std::invalid_argument("Elevation must be greater than -10,000 m");
    }
    pImpl->mElevation = elevation;
}

double UtahQuarry::getElevation() const noexcept
{
    return pImpl->mElevation;
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
