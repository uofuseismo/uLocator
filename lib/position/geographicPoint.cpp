#include "uLocator/position/geographicPoint.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "lonTo180.hpp"

using namespace ULocator::Position;

class IGeographicPoint::IGeographicPointImpl
{
public:
    IGeographicPointImpl(const IGeographicRegion &region) :
        mRegion(region.clone())
    {
    }
    IGeographicPointImpl(const IGeographicPointImpl &impl)
    {
        if (impl.mRegion != nullptr)
        {
            mRegion = impl.mRegion->clone();
        }
        mLatitude = impl.mLatitude;
        mLongitude = impl.mLongitude;
        mX = impl.mX;
        mY = impl.mY;
        mHaveGeographicCoordinates = impl.mHaveGeographicCoordinates;
    }
    std::unique_ptr<IGeographicRegion> mRegion{nullptr};
    double mLatitude{0};
    double mLongitude{0};
    double mX{0};
    double mY{0};
    bool mHaveGeographicCoordinates{false};
};

/// Constructor
IGeographicPoint::IGeographicPoint(
    const IGeographicRegion &region) :
    pImpl(std::make_unique<IGeographicPointImpl> (region))
{
}

/// Destructor
IGeographicPoint::~IGeographicPoint() = default;

/// Geographic coordinates set?
bool IGeographicPoint::haveGeographicCoordinates() const noexcept
{
    return pImpl->mHaveGeographicCoordinates;
}

/// Local coordinates
std::pair<double, double> IGeographicPoint::getLocalCoordinates() const
{
    if (!haveGeographicCoordinates())
    {
        throw std::runtime_error("Geographic coordinates not set");
    }
    return std::pair {pImpl->mX, pImpl->mY};
}

/// Geographic coordinates
std::pair<double, double> IGeographicPoint::getGeographicCoordinates() const
{
    if (!haveGeographicCoordinates())
    {
        throw std::runtime_error("Geographic coordinates not set");
    }
    return std::pair {pImpl->mLatitude, pImpl->mLongitude};
}

/// Set geographic point
void IGeographicPoint::setGeographicCoordinates(
    const double latitude, const double longitudeIn)
{
    if (latitude < -90 || latitude > 90)
    {
        throw std::invalid_argument("Latitude must be in range [-90,90]");
    }
    auto longitude = ::lonTo180(longitudeIn);
    auto [x, y]
        = pImpl->mRegion->geographicToLocalCoordinates(latitude, longitude);
    pImpl->mLatitude = latitude;
    pImpl->mLongitude = longitude;
    pImpl->mX = x;
    pImpl->mY = y;
    pImpl->mHaveGeographicCoordinates = true;
}

std::unique_ptr<IGeographicRegion> IGeographicPoint::getGeographicRegion() const
{
    return pImpl->mRegion->clone();
}

/*
/// Clone
std::unique_ptr<IGeographicPoint> IGeographicPoint::clone() const
{
    std::unique_ptr<IGeographicPoint> point
        = std::make_unique<IGeographicPoint> (*this);
    return point;
}
*/
