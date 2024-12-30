#include <string>
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/utahRegion.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace ULocator;

TEST_CASE("ULocator", "[station]")
{
    Position::UtahRegion utah;
    const std::string network{"UU"};
    const std::string name{"FPU"};
    constexpr double latitude{41.02633};
    constexpr double longitude{-111.83683};
    constexpr double elevation{2816};
    Position::WGS84 position{latitude, longitude};
    Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(position, utah);
    station.setElevation(elevation);

    Station sCopy(station);
    REQUIRE(sCopy.getNetwork() == network);
    REQUIRE(sCopy.getName() == name);
    REQUIRE_THAT(sCopy.getElevation(),
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    auto positionBack = sCopy.getGeographicPosition();
    REQUIRE_THAT(positionBack.getLatitude(),
                 Catch::Matchers::WithinAbs(latitude, 1.e-8));
    REQUIRE_THAT(positionBack.getLongitude(),
                 Catch::Matchers::WithinAbs(longitude, 1.e-8));
    auto localCoordinates = sCopy.getLocalCoordinates();
    REQUIRE_THAT(localCoordinates.first,
                 Catch::Matchers::WithinAbs(-28327.81065000718, 1.e-3));
    REQUIRE_THAT(localCoordinates.second,
                 Catch::Matchers::WithinAbs(  155642.5444867647, 1.e-3));
}

