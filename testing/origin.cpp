#include <string>
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/utahRegion.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace ULocator;

namespace
{

ULocator::Arrival toArrival(const std::string &network,
                            const std::string &name,
                            const std::string &phase,
                            const double latitude,
                            const double longitude,
                            const double elevation,
                            const double arrivalTime,
                            const double standardError, 
                            const int utmZone = 12) 
{
    ULocator::Position::UtahRegion utah;
    ULocator::Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        ULocator::Position::WGS84 {latitude, longitude, utmZone}, utah);
    station.setElevation(elevation);

    ULocator::Arrival arrival;
    arrival.setStation(station);
    arrival.setPhase(phase);
    arrival.setTime(arrivalTime);
    arrival.setStandardError(standardError);
    if (phase == "P")
    {
        arrival.setResidual(0.1);
    }
    else
    {
        arrival.setResidual(0.2);
    }
    return arrival;
}

/// Magna event uu60484462
std::vector<ULocator::Arrival> createArrivals()
{
    std::vector<ULocator::Arrival> arrivals;
    arrivals.push_back(::toArrival("UU", "ASU7", "P", 38.552595, -112.220772, 1609, 1646479476.6776488,   0.03));
    arrivals.push_back(::toArrival("LB", "MVU" , "P", 38.5037,   -112.212303, 2239, 1646479477.4650624,   0.03));
    arrivals.push_back(::toArrival("UU", "ASU7", "S", 38.552595, -112.220772, 1609, 1646479477.5781794,   0.15));
    arrivals.push_back(::toArrival("UU", "TCRU", "P", 38.6095,   -112.4472,   2293, 1646479479.779858,    0.03));
    arrivals.push_back(::toArrival("UU", "TCRU", "S", 38.6095,   -112.4472,   2293, 1646479483.044294,    0.15));
    arrivals.push_back(::toArrival("UU", "WCU",  "P", 38.96467,  -112.09067,  2673, 1646479483.8797174,   0.06));
    arrivals.push_back(::toArrival("UU", "NMU",  "P", 38.5165,   -112.85,     1853, 1646479485.68,        0.06));
    arrivals.push_back(::toArrival("UU", "ECUT", "P", 39.171697, -112.133201, 2136, 1646479487.7271104,   0.15));
    arrivals.push_back(::toArrival("UU", "DWU",  "P", 38.10534,  -112.9975,   2270, 1646479490.553661,    0.15));
    arrivals.push_back(::toArrival("UU", "NMU",  "S", 38.5165,   -112.85,     1853, 1646479492.619098,    0.06));
    arrivals.push_back(::toArrival("UU", "CVRU", "P", 38.9176,  -111.1716,    1912, 1646479492.699132,    0.15));
    arrivals.push_back(::toArrival("UU", "SWUT", "P", 39.3286,  -113.1954,    1644, 1646479496.371336,    0.15));
    return arrivals;
}

}

TEST_CASE("ULocator", "[origin]")
{
    Origin origin;
    constexpr int64_t identifier{60484462};
    constexpr double latitude{38.562};
    constexpr double longitude{-112.220};
    constexpr double depth{2000};
    constexpr double time{1646479475};
    auto arrivals = ::createArrivals();
 
    REQUIRE_NOTHROW(origin.setEpicenter(Position::WGS84 {latitude, longitude}));
    REQUIRE_NOTHROW(origin.setDepth(depth));
    REQUIRE_NOTHROW(origin.setTime(time));
    REQUIRE_NOTHROW(origin.setArrivals(arrivals));
    origin.setIdentifier(identifier);
    origin.setEventType(Origin::EventType::Earthquake);
 
    SECTION("copy constructor")
    {
    Origin copy(origin);
    REQUIRE_THAT(copy.getEpicenter().getLatitude(),
                 Catch::Matchers::WithinAbs(latitude,  1.e-10));
    REQUIRE_THAT(copy.getEpicenter().getLongitude(),
                 Catch::Matchers::WithinAbs(longitude, 1.e-10));
    REQUIRE_THAT(copy.getDepth(),
                 Catch::Matchers::WithinAbs(depth, 1.e-10));
    REQUIRE(copy.getIdentifier() == identifier);
    REQUIRE(copy.getEventType() == Origin::EventType::Earthquake);
    REQUIRE_THAT(copy.getTime(),
                 Catch::Matchers::WithinAbs(time, 1.e-4));
    REQUIRE(copy.getArrivals().size() == arrivals.size());
    REQUIRE_THAT(copy.getWeightedRootMeanSquaredError(),
                 Catch::Matchers::WithinAbs(0.12139539573337683, 1.e-12));
    REQUIRE_THAT(copy.getAzimuthalGap(),
                 Catch::Matchers::WithinAbs(107.820850789, 1.e-5));
    REQUIRE_THAT(copy.getNearestStationDistance(),
                 Catch::Matchers::WithinAbs(1046.187573, 1.e-5));
    }

    SECTION("clear")
    {
    auto copy = origin;
    copy.clear();
    REQUIRE_FALSE(copy.haveEpicenter());
    REQUIRE_FALSE(copy.haveDepth());
    REQUIRE_FALSE(copy.haveTime());
    REQUIRE(copy.getArrivals().empty());
    REQUIRE_FALSE(copy.haveAzimuthalGap());
    REQUIRE_FALSE(copy.haveNearestStationDistance());
    }
}

