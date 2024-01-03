#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#include "uLocator/position/ynp.hpp"

using namespace ULocator::Position;

#define MINIMUM_LATITUDE 43.5
#define MAXIMUM_LATITUDE 45.5
#define CENTER_LATITUDE  44.5
#define MINIMUM_LONGITUDE -111.75
#define MAXIMUM_LONGITUDE -109.75
#define CENTER_LONGITUDE  -110.75
#define MINIMUM_X -78154.09571331975
#define MAXIMUM_X  78154.09571331974
#define MINIMUM_Y -110611.925792941
#define MAXIMUM_Y  111604.1838556079

class YNP::YNPImpl
{
public:
/*
    YNPImpl()
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
YNP::YNP() :
    pImpl(std::make_unique<YNPImpl> ())
{
}

/// Copy constructor
YNP::YNP(const YNP &ynp)
{
    *this = ynp;
}

/// Move constructor
YNP::YNP(YNP &&ynp) noexcept
{
    *this = std::move(ynp);
}

/// Copy assignment
YNP& YNP::operator=(const YNP &ynp)
{
    if (&ynp == this){return *this;}
    pImpl = std::make_unique<YNPImpl> (*ynp.pImpl);
    return *this;
}

/// Move assignment
YNP& YNP::operator=(YNP &&ynp) noexcept
{
    if (&ynp == this){return *this;}
    pImpl = std::move(ynp.pImpl);
    return *this;
}

/// X-extent
std::pair<double, double> YNP::getExtentInX() const noexcept
{
    return std::pair {MINIMUM_X, MAXIMUM_X}; 
}

/// Y-extent
std::pair<double, double> YNP::getExtentInY() const noexcept
{
    return std::pair {MINIMUM_Y, MAXIMUM_Y};
}

/// Forward transformation
std::pair<double, double> 
YNP::geographicToLocalCoordinates(const double latitude,
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
YNP::localToGeographicCoordinates(const double x, const double y) const
{
    double latitude, longitude, h;
    pImpl->mProjection.Reverse(x, y, 0.0, latitude, longitude, h);
    return std::pair {latitude, longitude}; 
}


/// Destructor.
YNP::~YNP() = default;

/// Clone
std::unique_ptr<IGeographicRegion> YNP::clone() const
{
    std::unique_ptr<IGeographicRegion> region
        = std::make_unique<YNP> (*this);
    return region;
}
