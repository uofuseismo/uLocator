#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/position/utahQuarry.hpp"
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

TEST(ULocatorPosition, UtahRegion)
{
    Position::UtahRegion utah;
    const double centerLatitude{39.625};
    const double centerLongitude{-111.5};
    auto localCoordinates
        = utah.geographicToLocalCoordinates(centerLatitude, centerLongitude);
    EXPECT_NEAR(localCoordinates.first,  0, 1.e-10);
    EXPECT_NEAR(localCoordinates.second, 0, 1.e-10);
    auto geographicCoordinates
        = utah.localToGeographicCoordinates(localCoordinates.first,
                                            localCoordinates.second); 
    EXPECT_NEAR(geographicCoordinates.first,  centerLatitude,  1.e-10);
    EXPECT_NEAR(geographicCoordinates.second, centerLongitude, 1.e-10);
    auto minMaxX = utah.getExtentInX();
    auto minMaxY = utah.getExtentInY();
    EXPECT_NEAR(minMaxX.first,  -264866.0565037383, 1.e-3);
    EXPECT_NEAR(minMaxX.second,  264866.0565037381, 1.e-3);
    EXPECT_NEAR(minMaxY.first,  -369110.6893637062, 1.e-3);
    EXPECT_NEAR(minMaxY.second,  379402.8900178153, 1.e-3);
}

TEST(ULocatorPosition, YNPRegion)
{
    Position::YNPRegion ynp;
    const double centerLatitude{44.5};
    const double centerLongitude{-110.75};
    auto localCoordinates
        = ynp.geographicToLocalCoordinates(centerLatitude, centerLongitude);
    EXPECT_NEAR(localCoordinates.first,  0, 1.e-10);
    EXPECT_NEAR(localCoordinates.second, 0, 1.e-10);
    auto geographicCoordinates
        = ynp.localToGeographicCoordinates(localCoordinates.first,
                                           localCoordinates.second); 
    EXPECT_NEAR(geographicCoordinates.first,  centerLatitude,  1.e-10);
    EXPECT_NEAR(geographicCoordinates.second, centerLongitude, 1.e-10);
    auto minMaxX = ynp.getExtentInX();
    auto minMaxY = ynp.getExtentInY();
    EXPECT_NEAR(minMaxX.first,  -78154.09571331975, 1.e-3);
    EXPECT_NEAR(minMaxX.second,  78154.09571331974, 1.e-3);
    EXPECT_NEAR(minMaxY.first,  -110611.925792941,  1.e-3);
    EXPECT_NEAR(minMaxY.second,  111604.1838556079, 1.e-3);
}

TEST(ULocatorPosition, UtahQuarry)
{
    const std::string name{"Bingham"};
    constexpr double latitude{40.529800};
    constexpr double longitude{-112.147700};
    constexpr double elevation{1654};
    Position::UtahRegion utah;
    auto localCoordinates  
         = utah.geographicToLocalCoordinates(latitude, longitude);
 
    Position::UtahQuarry bingham;
    bingham.setName(name);
    bingham.setElevation(elevation);
    bingham.setGeographicCoordinates(latitude, longitude);

    auto copy = bingham;
    EXPECT_NEAR(copy.getGeographicCoordinates().first,  latitude,  1.e-10);
    EXPECT_NEAR(copy.getGeographicCoordinates().second, longitude, 1.e-10);
    EXPECT_NEAR(copy.getElevation(), elevation, 1.e-10);
    EXPECT_EQ(copy.getName(), name);
    EXPECT_NEAR(copy.getLocalCoordinates().first,
                localCoordinates.first, 1.e-10);
    EXPECT_NEAR(copy.getLocalCoordinates().second,
                localCoordinates.second, 1.e-10);

    auto move = std::move(bingham);
    EXPECT_NEAR(move.getGeographicCoordinates().first,  latitude,  1.e-10);
    EXPECT_NEAR(move.getGeographicCoordinates().second, longitude, 1.e-10);
    EXPECT_NEAR(move.getElevation(), elevation, 1.e-10);
    EXPECT_EQ(move.getName(), name);

    auto clone = move.clone();
    EXPECT_TRUE(clone->haveGeographicCoordinates());
    EXPECT_NEAR(clone->getGeographicCoordinates().first,  latitude,  1.e-10);
    EXPECT_NEAR(clone->getGeographicCoordinates().second, longitude, 1.e-10);
    EXPECT_NEAR(clone->getLocalCoordinates().first,
                localCoordinates.first, 1.e-10);
    EXPECT_NEAR(clone->getLocalCoordinates().second,
                localCoordinates.second, 1.e-10);
}

}
