#include <string>
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace ULocator;

namespace
{

TEST_CASE("ULocator", "[arrival]")
{
    const std::string network{"WY"};
    const std::string stationName{"YML"};
    constexpr double stationLatitude{44.60533};
    constexpr double stationLongitude{-110.64317};
    constexpr double stationElevation{2653.0};
    Position::YNPRegion ynp;
    Station station;
    station.setNetwork(network);
    station.setName(stationName);
    station.setGeographicPosition(Position::WGS84{stationLatitude, stationLongitude},  ynp);
    station.setElevation(stationElevation);
         
    constexpr double time{24};
    constexpr double residual{-1};
    constexpr double standardError{0.12};
    constexpr double azimuth{230};
    constexpr double backAzimuth{50};
    constexpr double distance{3000};
    constexpr int64_t identifier{3838};
    const Arrival::PhaseType phaseType{Arrival::PhaseType::S};
    Arrival arrival;

    arrival.setStation(station);
    arrival.setTime(time);
    arrival.setStandardError(standardError);
    arrival.setIdentifier(identifier);
    arrival.setPhase(phaseType);
    arrival.setResidual(residual);
    arrival.setDistance(distance);
    arrival.setAzimuth(azimuth);
    arrival.setBackAzimuth(backAzimuth);
    
    Arrival copy(arrival);
    REQUIRE_THAT(copy.getTime(),
                 Catch::Matchers::WithinAbs(time, 1.e-14));
    REQUIRE_THAT(copy.getStandardError(),
                 Catch::Matchers::WithinAbs(standardError, 1.e-14));
    REQUIRE_THAT(copy.getResidual(),
                 Catch::Matchers::WithinAbs(residual, 1.e-14));
    REQUIRE_THAT(copy.getDistance(),
                 Catch::Matchers::WithinAbs(distance, 1.e-12));
    REQUIRE_THAT(copy.getAzimuth(),
                 Catch::Matchers::WithinAbs(azimuth, 1.e-14));
    REQUIRE_THAT(copy.getBackAzimuth(),
                 Catch::Matchers::WithinAbs(backAzimuth, 1.e-14));
    REQUIRE(copy.getIdentifier() == identifier);
    REQUIRE(copy.getPhase() == "S");
    REQUIRE(copy.getStationReference().getNetwork() == network);
    REQUIRE(copy.haveStation());
    REQUIRE_THAT(copy.getStationReference().getElevation(),
                 Catch::Matchers::WithinAbs(stationElevation, 1.e-8));
    REQUIRE_THAT(copy.getStationReference().getLocalCoordinates().first,
                 Catch::Matchers::WithinAbs(8480.827504103971, 1.e-4));
    REQUIRE_THAT(copy.getStationReference().getLocalCoordinates().second,
                 Catch::Matchers::WithinAbs(11710.12450692187, 1.e-4));
    REQUIRE(copy.getStation().getName() == stationName);

    SECTION("clear")
    {
    arrival.clear();
    REQUIRE_FALSE(arrival.haveTime());
    }
}

}
