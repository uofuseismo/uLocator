#include "uLocator/position/wgs84.hpp"
#include <gtest/gtest.h>

using namespace ULocator;

namespace
{

TEST(ULocatorPosition, WGS84)
{
    // http://rcn.montana.edu/resources/Converter.aspx
    const double latitude{40.7608};
    const double longitude{-111.8910 + 360};
    const double easting{424795};
    const double northing{4512586};
    const int utmZone{12};
    const bool isNorth{true};
 
    Position::WGS84 position(latitude, longitude);
    EXPECT_TRUE(position.havePosition());
    EXPECT_NEAR(position.getLatitude(), latitude, 1.e-8);
    EXPECT_NEAR(position.getLongitude(), longitude - 360, 1.e-8);
    EXPECT_TRUE(position.isNorth());
    EXPECT_EQ(position.getUTMZone(), utmZone);
    EXPECT_NEAR(position.getEasting(),  easting,  1);
    EXPECT_NEAR(position.getNorthing(), northing, 1);

    
    Position::WGS84 utmPosition(utmZone, isNorth, easting, northing);
    EXPECT_TRUE(utmPosition.havePosition());
    EXPECT_NEAR(utmPosition.getEasting(),  easting,  1.e-6);
    EXPECT_NEAR(utmPosition.getNorthing(), northing, 1.e-6);
    EXPECT_NEAR(utmPosition.getLatitude(),    40.76080047521302, 1.e-5);
    EXPECT_NEAR(utmPosition.getLongitude(), -111.89099802607673, 1.e-5);
}

TEST(ULocatorPosition, DistanceAzimuthBackAzimuth)
{
    constexpr double azimuthRef{-101.69839545926366 + 360};
    constexpr double backAzimuthRef{-99.84855154050379 + 180};
    constexpr double gcDistanceRef{124.17888035704068};
    constexpr double distanceRef{13779227.281877534};
    Position::WGS84 source{19, 32.5};
    Position::WGS84 station{-20, 272.4};

    double greatCircleDistance, distance, azimuth, backAzimuth;
    Position::computeDistanceAzimuth(source, station,
                                     &greatCircleDistance,
                                     &distance,
                                     &azimuth,
                                     &backAzimuth);
    EXPECT_NEAR(greatCircleDistance, gcDistanceRef, 1.e-8);
    EXPECT_NEAR(distance, distanceRef, 1.e-7);
    EXPECT_NEAR(azimuth, azimuthRef, 1.e-10);
    EXPECT_NEAR(backAzimuth, backAzimuthRef, 1.e-10);
}

}
