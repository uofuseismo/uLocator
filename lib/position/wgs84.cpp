#include <cmath>
#ifndef NDEBUG
#include <cassert>
#endif
#include <GeographicLib/UTMUPS.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include "uLocator/position/wgs84.hpp"
#include "geodetic.hpp"
#include "lonTo180.hpp"

using namespace ULocator::Position;

class WGS84::WGS84Impl
{
public:
    WGS84Impl() = default;
    void forward()
    {
        try
        {
            ::GeographicLib::UTMUPS::Forward(mLatitude, mLongitude,
                                             mZone,
                                             mNorth,
                                             mEasting,
                                             mNorthing,
                                             mAlternateZone);
        }
        catch (...)
        {
            if (mAlternateZone !=-1)
            {
                mZone = mAlternateZone;
            }
            else
            {
                mZone
                   = GeographicLib::UTMUPS::StandardZone(mLatitude, mLongitude);
            }
            mNorth = true;
            if (mLatitude < 0){mNorth = false;}
            auto zone = mZone;
            if (!mNorth){zone =-zone;} 
            ::latitudeLongitudeToUTM(mLatitude, mLongitude, zone,
                                     &mEasting, &mNorthing);
         }
    }
    void reverse()
    {
        try
        {
            double gamma, k;
            GeographicLib::UTMUPS::Reverse(mZone, mNorth,
                                           mEasting, mNorthing,
                                           mLatitude, mLongitude,
                                           gamma, k);
        }
        catch (const std::exception &e)
        {
            auto zone = mZone;
            if (!mNorth){zone =-zone;}
            ::utmToLatitudeLongitude(mEasting, mNorthing, zone,
                                     &mLatitude, &mLongitude);
        }
    }
    double mLatitude{0};
    double mLongitude{0};
    double mEasting{0};
    double mNorthing{0}; 
    int mAlternateZone{-1};
    int mZone{-1};
    bool mHavePosition{false};
    bool mHaveLongitude{false};
    bool mNorth{true};
};

/// Constructor
WGS84::WGS84() :
    pImpl(std::make_unique<WGS84Impl> ())
{
}

/// Copy constructor
WGS84::WGS84(const WGS84 &position)
{
    *this = position; 
}

/// Move constructor
WGS84::WGS84(WGS84 &&position) noexcept
{
    *this = std::move(position);
}

/// Copy assignment
WGS84& WGS84::operator=(const WGS84 &position)
{
    if (&position == this){return *this;}
    pImpl = std::make_unique<WGS84Impl> (*position.pImpl);
    return *this;
}

/// Move assignment
WGS84& WGS84::operator=(WGS84 &&position) noexcept
{
    if (&position == this){return *this;}
    pImpl = std::move(position.pImpl);
    return *this;
}

/// Latitude and longitude
WGS84::WGS84(const double latitude, const double longitude,
             const int alternateZone) :
    pImpl(std::make_unique<WGS84Impl> ())
{
    pImpl->mAlternateZone = alternateZone;
    setLatitudeAndLongitude(latitude, longitude);
}

WGS84::WGS84(const int zone, const bool isNorth,
             const double easting, const double northing) :
    pImpl(std::make_unique<WGS84Impl> ())
{
    setEastingAndNorthing(zone, isNorth, easting, northing);
}

void WGS84::setLatitudeAndLongitude(const double latitude,
                                    const double longitude)
{
    if (latitude <-90 || latitude > 90)
    {
        throw std::runtime_error("Latitude not in bounds [-90,90]");
    }
    pImpl->mLatitude = latitude;
    pImpl->mLongitude = ::lonTo180(longitude);
    pImpl->mHavePosition = true;
    pImpl->forward();
}

void WGS84::setAlternateUTMZone(const int zone)
{
    if (zone < 1 || zone > 60)
    {
        if (zone != -1)
        {
            throw std::invalid_argument(
               "Zone must be -1 or in range [1,60]");
        }
    }
    pImpl->mAlternateZone = zone;
    if (havePosition())
    {
        pImpl->forward();
    }
}

double WGS84::getLatitude() const
{
    if (!havePosition()){throw std::runtime_error("Position not set");}
    return pImpl->mLatitude;
}

double WGS84::getLongitude() const
{
    if (!havePosition()){throw std::runtime_error("Position not set");}
    return pImpl->mLongitude;
}

bool WGS84::havePosition() const noexcept
{
    return pImpl->mHavePosition;
}

/// Northing and easting
void WGS84::setEastingAndNorthing(const int zone,
                                  const bool isNorth,
                                  const double easting,
                                  const double northing)
{
    if (zone <-1 || zone > 60)
    {
        throw std::invalid_argument("Zone must be in range [-1,60]");
    }
    pImpl->mNorthing = northing;
    pImpl->mEasting = easting; 
    pImpl->mZone = zone;
    pImpl->mNorth = isNorth;
    pImpl->reverse(); 
    pImpl->mHavePosition = true;
}

double WGS84::getEasting() const
{
    if (!havePosition()){throw std::runtime_error("Position not set");}
    return pImpl->mEasting;
}

double WGS84::getNorthing() const
{
    if (!havePosition()){throw std::runtime_error("Position not set");}
    return pImpl->mNorthing;
}

bool WGS84::isNorth() const
{
    if (!havePosition()){throw std::runtime_error("Position not set");}
    return pImpl->mNorth;
}

int WGS84::getAlternateUTMZone() const noexcept
{
    return pImpl->mAlternateZone;
}

int WGS84::getUTMZone() const
{
    if (!havePosition()){throw std::runtime_error("Position not set");}
    return pImpl->mZone;
}

/// Reset class
void WGS84::clear() noexcept
{
    pImpl = std::make_unique<WGS84Impl> ();
}

/// Destructor
WGS84::~WGS84() = default;

/// Distance between two wgs84 points
void ULocator::Position::computeDistanceAzimuth(
    const WGS84 &source,
    const WGS84 &station,
    double *mGreatCircleDistance,
    double *mDistance,
    double *mAzimuth,
    double *mBackAzimuth)
{
    if (!source.havePosition())
    {
        throw std::invalid_argument("Source position not set");
    } 
    if (!station.havePosition())
    {
        throw std::invalid_argument("Station position not set");
    }
    if (mGreatCircleDistance == nullptr &&
        mDistance == nullptr &&
        mAzimuth == nullptr &&
        mBackAzimuth == nullptr)
    {
        return;
    }
    GeographicLib::Geodesic geodesic{GeographicLib::Constants::WGS84_a(),
                                     GeographicLib::Constants::WGS84_f()};

    auto sourceLatitude   = source.getLatitude();
    auto sourceLongitude  = source.getLongitude(); 
    auto stationLatitude  = station.getLatitude();
    auto stationLongitude = station.getLongitude();
    double distance, azimuth, backAzimuth;
    double greatCircleDistance
         = geodesic.Inverse(sourceLatitude,  sourceLongitude,
                            stationLatitude, stationLongitude,
                            distance, azimuth, backAzimuth);
    // Translate azimuth from [-180,180] to [0,360].
    if (azimuth < 0){azimuth = azimuth + 360;}
    // Translate azimuth from [-180,180] to [0,360] then convert to a 
    // back-azimuth by subtracting 180, i.e., +180.
    backAzimuth = backAzimuth + 180;
    if (mGreatCircleDistance != nullptr)
    {
        *mGreatCircleDistance = greatCircleDistance;
    }
    if (mDistance != nullptr)
    {
        *mDistance = distance;
    }
    if (mAzimuth != nullptr)
    {
        *mAzimuth = azimuth;
    }
    if (mBackAzimuth != nullptr)
    {
        // Deal with issue if azimuth = 180
        if (backAzimuth >= 360)
        {
            backAzimuth = backAzimuth - 360;
        }
        *mBackAzimuth = backAzimuth;
    }
}
