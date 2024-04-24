#include "uLocator/position/knownUtahQuarry.hpp"
#include "uLocator/position/utahQuarry.hpp"

using namespace ULocator::Position;

class KnownUtahQuarry::KnownUtahQuarryImpl
{
public:
    KnownUtahQuarryImpl(const ULocator::Position::UtahQuarry &quarry) :
        mQuarry(quarry)
    {
    }
    ULocator::Position::UtahQuarry mQuarry;
};

KnownUtahQuarry::KnownUtahQuarry(const ULocator::Position::UtahQuarry &quarry) :
    pImpl(std::make_unique<KnownUtahQuarryImpl> (quarry))
{    
}

KnownUtahQuarry::KnownUtahQuarry(const KnownUtahQuarry &quarry)
{
    *this = quarry;
}

KnownUtahQuarry::KnownUtahQuarry(KnownUtahQuarry &&quarry) noexcept
{
    *this = std::move(quarry);
}
    
double KnownUtahQuarry::x() const
{
    return pImpl->mQuarry.getLocalCoordinates().first;
}
    
double KnownUtahQuarry::y() const
{
    return pImpl->mQuarry.getLocalCoordinates().second;
}
    
double KnownUtahQuarry::z() const
{
    return -pImpl->mQuarry.getElevation();
}

KnownUtahQuarry& KnownUtahQuarry::operator=(const KnownUtahQuarry &quarry)
{
    if (&quarry == this){return *this;}
    pImpl = std::make_unique<KnownUtahQuarryImpl> (*quarry.pImpl);
    return *this;
}

KnownUtahQuarry& KnownUtahQuarry::operator=(KnownUtahQuarry &&quarry) noexcept
{
    if (&quarry == this){return *this;}
    pImpl = std::move(quarry.pImpl);
    return *this;
}
    
KnownUtahQuarry::~KnownUtahQuarry() = default;


/// The quarry
ULocator::Position::UtahQuarry KnownUtahQuarry::getUtahQuarry() const noexcept
{
    return pImpl->mQuarry;
}

const ULocator::Position::UtahQuarry&
 KnownUtahQuarry::getUtahQuarryReference() const noexcept
{
    return *&pImpl->mQuarry;
}

/// Clone
std::unique_ptr<IKnownLocalLocation> KnownUtahQuarry::clone() const
{
    std::unique_ptr<IKnownLocalLocation> region
        = std::make_unique<KnownUtahQuarry> (*this);
    return region;
}

