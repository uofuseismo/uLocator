#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/position/utahQuarry.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace ULocator;

TEST_CASE("ULocator::Position", "[position]")
{
    // http://rcn.montana.edu/resources/Converter.aspx
    const double latitude{40.7608};
    const double longitude{-111.8910 + 360};
    const double easting{424795};
    const double northing{4512586};
    const int utmZone{12};
    const bool isNorth{true};
 
    Position::WGS84 position(latitude, longitude);
    REQUIRE(position.havePosition());
    REQUIRE_THAT(position.getLatitude(),
                 Catch::Matchers::WithinAbs(latitude, 1.e-8));
    REQUIRE_THAT(position.getLongitude(),
                 Catch::Matchers::WithinAbs(longitude - 360, 1.e-8));
    REQUIRE(position.isNorth());
    REQUIRE(position.getUTMZone() == utmZone);
    REQUIRE_THAT(position.getEasting(),
                 Catch::Matchers::WithinAbs(easting,  1));
    REQUIRE_THAT(position.getNorthing(),
                 Catch::Matchers::WithinAbs(northing, 1));

    
    Position::WGS84 utmPosition(utmZone, isNorth, easting, northing);
    REQUIRE(utmPosition.havePosition());
    REQUIRE_THAT(utmPosition.getEasting(),  
                 Catch::Matchers::WithinAbs(easting,  1.e-6));
    REQUIRE_THAT(utmPosition.getNorthing(),
                 Catch::Matchers::WithinAbs(northing, 1.e-6));
    REQUIRE_THAT(utmPosition.getLatitude(),
                 Catch::Matchers::WithinAbs(40.76080047521302, 1.e-5));
    REQUIRE_THAT(utmPosition.getLongitude(),
                 Catch::Matchers::WithinAbs(-111.89099802607673, 1.e-5));
}

TEST_CASE("ULocator::Position", "[distanceAzimuthBackazimuth]")
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
    REQUIRE_THAT(greatCircleDistance,
                 Catch::Matchers::WithinAbs(gcDistanceRef, 1.e-8));
    REQUIRE_THAT(distance,
                 Catch::Matchers::WithinAbs(distanceRef, 1.e-7));
    REQUIRE_THAT(azimuth,
                 Catch::Matchers::WithinAbs(azimuthRef, 1.e-10));
    REQUIRE_THAT(backAzimuth,
                 Catch::Matchers::WithinAbs(backAzimuthRef, 1.e-10));
}

TEST_CASE("ULocator::Position", "[utahRegion]")
{
    Position::UtahRegion utah;
    const double centerLatitude{39.625};
    const double centerLongitude{-111.5};
    auto localCoordinates
        = utah.geographicToLocalCoordinates(centerLatitude, centerLongitude);
    REQUIRE_THAT(localCoordinates.first,
                 Catch::Matchers::WithinAbs(0, 1.e-10));
    REQUIRE_THAT(localCoordinates.second,
                 Catch::Matchers::WithinAbs(0, 1.e-10));
    auto geographicCoordinates
        = utah.localToGeographicCoordinates(localCoordinates.first,
                                            localCoordinates.second); 
    REQUIRE_THAT(geographicCoordinates.first,
                 Catch::Matchers::WithinAbs(centerLatitude,  1.e-10));
    REQUIRE_THAT(geographicCoordinates.second,
                 Catch::Matchers::WithinAbs(centerLongitude, 1.e-10));
    auto minMaxX = utah.getExtentInX();
    auto minMaxY = utah.getExtentInY();
    REQUIRE_THAT(minMaxX.first,
                 Catch::Matchers::WithinAbs( -264866.0565037383, 1.e-3));
    REQUIRE_THAT(minMaxX.second,
                 Catch::Matchers::WithinAbs(  264866.0565037381, 1.e-3));
    REQUIRE_THAT(minMaxY.first,
                 Catch::Matchers::WithinAbs(  -369110.6893637062, 1.e-3));
    REQUIRE_THAT(minMaxY.second,
                 Catch::Matchers::WithinAbs(  379402.8900178153, 1.e-3));
}

TEST_CASE("ULocator::Position", "[ynpRegion]")
{
    Position::YNPRegion ynp;
    const double centerLatitude{44.5};
    const double centerLongitude{-110.75};
    auto localCoordinates
        = ynp.geographicToLocalCoordinates(centerLatitude, centerLongitude);
    REQUIRE_THAT(localCoordinates.first,
                 Catch::Matchers::WithinAbs(0, 1.e-10));
    REQUIRE_THAT(localCoordinates.second,
                 Catch::Matchers::WithinAbs(0, 1.e-10));
    auto geographicCoordinates
        = ynp.localToGeographicCoordinates(localCoordinates.first,
                                           localCoordinates.second); 
    REQUIRE_THAT(geographicCoordinates.first,
                 Catch::Matchers::WithinAbs(centerLatitude,  1.e-10));
    REQUIRE_THAT(geographicCoordinates.second,
                 Catch::Matchers::WithinAbs(centerLongitude, 1.e-10));
    auto minMaxX = ynp.getExtentInX();
    auto minMaxY = ynp.getExtentInY();
    REQUIRE_THAT(minMaxX.first,
                 Catch::Matchers::WithinAbs(-78154.09571331975, 1.e-3));
    REQUIRE_THAT(minMaxX.second,
                 Catch::Matchers::WithinAbs( 78154.09571331974, 1.e-3));
    REQUIRE_THAT(minMaxY.first,
                 Catch::Matchers::WithinAbs(-110611.925792941,  1.e-3));
    REQUIRE_THAT(minMaxY.second,
                 Catch::Matchers::WithinAbs( 111604.1838556079, 1.e-3));
}

TEST_CASE("ULocator::Position", "[utahQuarry]")
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
    REQUIRE_THAT(copy.getGeographicCoordinates().first,
                 Catch::Matchers::WithinAbs(latitude,  1.e-10));
    REQUIRE_THAT(copy.getGeographicCoordinates().second,
                 Catch::Matchers::WithinAbs(longitude, 1.e-10));
    REQUIRE_THAT(copy.getElevation(),
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    REQUIRE(copy.getName() == name);
    REQUIRE_THAT(copy.getLocalCoordinates().first,
                 Catch::Matchers::WithinAbs(localCoordinates.first, 1.e-10));
    REQUIRE_THAT(copy.getLocalCoordinates().second,
                 Catch::Matchers::WithinAbs(localCoordinates.second, 1.e-10));

    SECTION("move assignment")
    {
    auto move = std::move(bingham);
    REQUIRE_THAT(move.getGeographicCoordinates().first,
                 Catch::Matchers::WithinAbs(latitude,  1.e-10));
    REQUIRE_THAT(move.getGeographicCoordinates().second,
                 Catch::Matchers::WithinAbs(longitude, 1.e-10));
    REQUIRE_THAT(move.getElevation(),
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    REQUIRE(move.getName() == name);
    }

    SECTION("clone")
    {
    auto clone = bingham.clone();
    REQUIRE(clone->haveGeographicCoordinates());
    REQUIRE_THAT(clone->getGeographicCoordinates().first,
                 Catch::Matchers::WithinAbs(latitude,  1.e-10));
    REQUIRE_THAT(clone->getGeographicCoordinates().second,
                 Catch::Matchers::WithinAbs(longitude, 1.e-10));
    REQUIRE_THAT(clone->getLocalCoordinates().first,
                 Catch::Matchers::WithinAbs(localCoordinates.first, 1.e-10));
    REQUIRE_THAT(clone->getLocalCoordinates().second,
                 Catch::Matchers::WithinAbs(localCoordinates.second, 1.e-10));
    }

}
