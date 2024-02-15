#include <cmath>
#include <string>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include "uLocator/position/utahRegion.hpp"

using namespace ULocator::Position;

#define MINIMUM_LATITUDE 36.25
#define MAXIMUM_LATITUDE 43.00
#define CENTER_LATITUDE  39.625
#define MINIMUM_LONGITUDE -114.75
#define MAXIMUM_LONGITUDE -108.25
#define CENTER_LONGITUDE  -111.5
#define MINIMUM_X -264866.0565037383
#define MAXIMUM_X  264866.0565037381
#define MINIMUM_Y -369110.6893637062
#define MAXIMUM_Y  379402.8900178153

class UtahRegion::UtahRegionImpl
{
public:
    /*
    UtahRegionImpl()
    {
        double x, y, z, lat, lon, h;
        mProjection.Forward(MINIMUM_LATITUDE, MINIMUM_LONGITUDE, 0, x, y, z);
        mMinimumX = x;
        mMinimumY = y;
        mProjection.Forward(MAXIMUM_LATITUDE, MINIMUM_LONGITUDE, 0, x, y, z);
        mMinimumX = std::max(mMinimumX, x);
        mMaximumY = y;
        
        mProjection.Forward(MINIMUM_LATITUDE, MAXIMUM_LONGITUDE, 0, x, y, z);
        mMaximumX = x;
        mMinimumY = std::max(mMinimumY, y);
        mProjection.Forward(MAXIMUM_LATITUDE, MAXIMUM_LONGITUDE, 0, x, y, z); 
        mMaximumX = std::min(mMaximumX, x);
        mMaximumY = std::min(mMaximumY, y);

        std::cout << std::setprecision(16) << mMinimumX << " " << mMaximumX << " " << mMinimumY << " " << mMaximumY << std::endl;
         
        mProjection.Reverse(mMinimumX, mMinimumY, 0, lat, lon, h);
        std::cout << lat << " " << lon << std::endl;
        mProjection.Reverse(mMaximumX, mMinimumY, 0, lat, lon, h);
        std::cout << lat << " " << lon << std::endl;
        mProjection.Reverse(mMinimumX, mMaximumY, 0, lat, lon, h); 
        std::cout << lat << " " << lon << std::endl;
        mProjection.Reverse(mMaximumX, mMaximumY, 0, lat, lon, h); 
        std::cout << lat << " " << lon << std::endl;
    }
    */
    GeographicLib::Geocentric mEarth{GeographicLib::Constants::WGS84_a(),
                                     GeographicLib::Constants::WGS84_f()};
    GeographicLib::LocalCartesian mProjection{CENTER_LATITUDE,
                                              CENTER_LONGITUDE,
                                              0.0, // Height above ellipsoid at origin (meters)
                                              mEarth};
    double mMinimumX{MINIMUM_X};
    double mMaximumX{MAXIMUM_X};
    double mMinimumY{MINIMUM_Y};
    double mMaximumY{MAXIMUM_Y};
};

/// Constructor.
UtahRegion::UtahRegion() :
    pImpl(std::make_unique<UtahRegionImpl> ())
{
}

/// Copy constructor
UtahRegion::UtahRegion(const UtahRegion &utah)
{
    *this = utah;
}

/// Move constructor
UtahRegion::UtahRegion(UtahRegion &&utah) noexcept
{
    *this = std::move(utah);
}

/// Copy assignment
UtahRegion& UtahRegion::operator=(const UtahRegion &utah)
{
    if (&utah == this){return *this;}
    pImpl = std::make_unique<UtahRegionImpl> (*utah.pImpl);
    return *this;
}

/// Move assignment
UtahRegion& UtahRegion::operator=(UtahRegion &&utah) noexcept
{
    if (&utah == this){return *this;}
    pImpl = std::move(utah.pImpl);
    return *this;
}

/// X-extent
std::pair<double, double> UtahRegion::getExtentInX() const noexcept
{
    return std::pair {MINIMUM_X, MAXIMUM_X}; 
}

/// Y-extent
std::pair<double, double> UtahRegion::getExtentInY() const noexcept
{
    return std::pair {MINIMUM_Y, MAXIMUM_Y};
}

/// Forward transformation
std::pair<double, double> 
UtahRegion::geographicToLocalCoordinates(const double latitude,
                                         const double longitude) const
{
    if (latitude < -90 || latitude > 90)
    {
        throw std::invalid_argument("Latitude " + std::to_string(latitude)
                                  + " must be in range [-90,90]");
    }
    double x, y, z;
    pImpl->mProjection.Forward(latitude, longitude, 0.0, x, y, z); 
    return std::pair {x, y};
}


/// Reverse transformation
std::pair<double, double> 
UtahRegion::localToGeographicCoordinates(const double x, const double y) const
{
    double latitude, longitude, h;
    pImpl->mProjection.Reverse(x, y, 0.0, latitude, longitude, h);
    return std::pair {latitude, longitude}; 
}


/// Destructor.
UtahRegion::~UtahRegion() = default;

/// Clone
std::unique_ptr<IGeographicRegion> UtahRegion::clone() const
{
    std::unique_ptr<IGeographicRegion> region
        = std::make_unique<UtahRegion> (*this);
    return region;
}

/*
#include <GeographicLib/Geodesic.hpp>
int main()
{
    Utah utah;
    double stGeorgeLatitude{37.0965};
    double stGeorgeLongitude{-113.5684};
    double slcLatitude{40.7608};
    double slcLongitude{-111.8910};
    auto [xsg, ysg] = utah.geographicToLocalCoordinates(stGeorgeLatitude, stGeorgeLongitude);
    auto [xslc, yslc] = utah.geographicToLocalCoordinates(slcLatitude, slcLongitude);
    auto distance = std::hypot(xslc - xsg, yslc - ysg);
    std::cout << std::setprecision(12) << distance << std::endl;

    GeographicLib::Geodesic geodesic{GeographicLib::Constants::WGS84_a(),
                                     GeographicLib::Constants::WGS84_f()};
    double distanceCorrect, azimuth, backAzimuth;
    double greatCircleDistance
         = geodesic.Inverse(slcLatitude,      slcLongitude,
                            stGeorgeLatitude, stGeorgeLongitude,
                            distanceCorrect, azimuth, backAzimuth);
    if (azimuth < 0){azimuth = azimuth + 360;}
    backAzimuth = backAzimuth + 180;

    std::cout << distance << " " << distanceCorrect << std::endl; 
    std::cout << azimuth << " " << 90 - std::atan2(ysg - yslc, xsg - xslc)*180/M_PI << std::endl; // target - origin
    std::cout << backAzimuth << " " << 90 - std::atan2(yslc - ysg, xslc - xsg)*180/M_PI << std::endl;
}
*/
