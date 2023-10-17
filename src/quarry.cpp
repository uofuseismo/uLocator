#include <string>
#include "uLocator/quarry.hpp"
#include "uLocator/position/wgs84.hpp"

using namespace ULocator;

class Quarry::QuarryImpl
{
public:
    Position::WGS84 mPosition;
    std::string mName;
};

/// Constructor
Quarry::Quarry() :
    pImpl(std::make_unique<QuarryImpl> ())
{
}

/// Constructor
Quarry::Quarry(const Position::WGS84 &position, const std::string &name) :
    pImpl(std::make_unique<QuarryImpl> ())
{
    setGeographicPosition(position);
    setName(name);
}

/// Copy constructor
Quarry::Quarry(const Quarry &quarry)
{
    *this = quarry;
}

/// Move constructor
Quarry::Quarry(Quarry &&quarry) noexcept
{
    *this = std::move(quarry);
}

/// Copy assignment
Quarry& Quarry::operator=(const Quarry &quarry)
{
    if (&quarry == this){return *this;}
    pImpl = std::make_unique<QuarryImpl> (*quarry.pImpl);
    return *this;
}

/// Move assignment
Quarry& Quarry::operator=(Quarry &&quarry) noexcept
{
    if (&quarry == this){return *this;}
    pImpl = std::move(quarry.pImpl);
    return *this;
}

/// Reset class
void Quarry::clear() noexcept
{
    pImpl = std::make_unique<QuarryImpl> ();
}

/// Destructor
Quarry::~Quarry() = default;

/// Name
void Quarry::setName(const std::string &name)
{
    pImpl->mName = name;
}

std::string Quarry::getName() const noexcept
{
    return pImpl->mName;
}

/// Position
void Quarry::setGeographicPosition(const Position::WGS84 &position)
{
    if (!position.havePosition())
    {
        throw std::invalid_argument("Quarry position not set");
    }
    pImpl->mPosition = position;
}

Position::WGS84 Quarry::getGeographicPosition() const
{
    if (!pImpl->mPosition.havePosition())
    {
        throw std::invalid_argument("Quarry position not set");
    }
    return pImpl->mPosition;
}

bool Quarry::haveGeographicPosition() const noexcept
{
    return pImpl->mPosition.havePosition();
}
