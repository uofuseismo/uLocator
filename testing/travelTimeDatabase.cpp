#include <vector>
#include "travelTimeDatabase.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/position/wgs84.hpp"
#include "originTime.hpp"
#include <gtest/gtest.h>

namespace
{

// Create a list of arrivals for this event
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
    return arrival;
}

/// Quarry blast uu60020257
ULocator::Origin createQuarryBlast60020257()
{
    std::vector<ULocator::Arrival> arrivals;
    arrivals.push_back(::toArrival("UU", "CWU",  "P", 40.44584,-112.10217, 1945.0, 1366569002.4887397, 0.03));
    arrivals.push_back(::toArrival("UU", "MID",  "P", 40.51733,-112.25467, 1722.0, 1366569002.8587527, 0.06));
    arrivals.push_back(::toArrival("UU", "NOQ",  "P", 40.6525, -112.12033, 1622.0, 1366569004.0445964, 0.03));
    arrivals.push_back(::toArrival("UU", "WTU",  "P", 40.45483,-111.9535,  1552.0, 1366569004.5588472, 0.03));
    arrivals.push_back(::toArrival("UU", "GMU",  "P", 40.5755, -111.76317, 1829.0, 1366569006.866271,  0.15));
    arrivals.push_back(::toArrival("UU", "MHD",  "P", 40.66067,-111.80083, 1597.0, 1366569007.6222286, 0.30));
    arrivals.push_back(::toArrival("UU", "SAIU", "P", 40.85483,-112.1815,  1384.0, 1366569007.8094156, 0.06));
    arrivals.push_back(::toArrival("UU", "RBU",  "P", 40.78083,-111.80833, 1676.0, 1366569008.9965146, 0.15));
    arrivals.push_back(::toArrival("UU", "SNUT", "P", 40.88567,-112.509,   1652.0, 1366569009.8162184, 0.15));
    arrivals.push_back(::toArrival("UU", "NAIU", "P", 41.01617,-112.228,   1472.0, 1366569010.6781356, 0.15));
    arrivals.push_back(::toArrival("UU", "NLU",  "P", 39.95483,-112.075,   2036.0, 1366569011.983079,  0.30));
    arrivals.push_back(::toArrival("US", "DUG",  "P", 40.195,  -112.8133,  1477.0, 1366569012.6228065, 0.30));
    arrivals.push_back(::toArrival("UU", "SPU",  "P", 41.30867,-112.44917, 2086.0, 1366569016.7019002, 0.15)); // 13

    ULocator::Origin origin;
    origin.setTime(1366569000.31);
    origin.setEpicenter(ULocator::Position::WGS84 {40.5023333,-112.1658333});
    origin.setDepth(-2000);
    origin.setArrivals(arrivals);
    origin.setEventType(ULocator::Origin::EventType::QuarryBlast);
    return origin;
}

/// Quarry blast uu60484742
ULocator::Origin createQuarryBlast60484742()
{
    std::vector<ULocator::Arrival> arrivals;
    arrivals.push_back(::toArrival("UU", "CWU",  "P", 40.44584,-112.10217, 1945.0, 1646684879.079617,   0.06));
    arrivals.push_back(::toArrival("UU", "MID",  "P", 40.51733,-112.25467, 1722.0, 1646684879.3381984,  0.03));
    arrivals.push_back(::toArrival("UU", "NOQ",  "P", 40.6525, -112.12033, 1622.0, 1646684880.4258294,  0.03));
    arrivals.push_back(::toArrival("UU", "WTU",  "P", 40.45483,-111.9535,  1552.0, 1646684880.785694,   0.15));
    arrivals.push_back(::toArrival("UU", "GMU",  "P", 40.5755, -111.76317, 1829.0, 1646684883.5544589,  0.15));
    arrivals.push_back(::toArrival("UU", "SAIU", "P", 40.85483,-112.1815,  1384.0, 1646684883.9708252,  0.03));
    arrivals.push_back(::toArrival("UU", "RBU",  "P", 40.78083,-111.80833, 1676.0, 1646684885.032984,   0.15));
    arrivals.push_back(::toArrival("UU", "SNUT", "P", 40.88567,-112.509,   1652.0, 1646684886.6219738,  0.15));
    arrivals.push_back(::toArrival("UU", "NAIU", "P", 41.01617,-112.228,   1472.0, 1646684887.0673397,  0.15));
    arrivals.push_back(::toArrival("UU", "NLU",  "P", 39.95483,-112.075,   2036.0, 1646684888.7377942,  0.15));
    arrivals.push_back(::toArrival("UU", "MOUT", "P", 41.19894,-111.87953, 2748.0, 1646684891.405937,   0.15)); // + 1
    arrivals.push_back(::toArrival("UU", "SPU",  "P", 41.30867,-112.44917, 2086.0, 1646684893.4537792,  0.15));
    arrivals.push_back(::toArrival("UU", "HVU",  "P", 41.77967,-112.775,   1609.0, 1646684902.3164327,  0.06)); // + 1

    ULocator::Origin origin;
    origin.setTime(1646684876.8099976);
    origin.setEpicenter(ULocator::Position::WGS84 {40.5151667,-112.1455});
    origin.setDepth(-2000);
    origin.setArrivals(arrivals);
    origin.setEventType(ULocator::Origin::EventType::QuarryBlast);
    return origin;
}

TEST(ULocator, UtahQuarryBlastTravelTimeDatabase)
{
    ::UtahQuarryBlastTravelTimeDatabase db;

    std::array<ULocator::Optimizers::IOptimizer::Norm, 3> norms
    {
        ULocator::Optimizers::IOptimizer::Norm::LeastSquares,
        ULocator::Optimizers::IOptimizer::Norm::L1,
        ULocator::Optimizers::IOptimizer::Norm::Lp 
    };
    constexpr double p{1.5};
    constexpr double timeWindow{140};
    constexpr bool applyCorrection{true};
    auto origin = createQuarryBlast60484742();
    for (const auto &norm : norms)
    {
        auto [bestQuarryIndex, bestQuarry, bestObjectiveFunction]
            = db.findBestQuarry(origin, norm, p, timeWindow, applyCorrection);
        EXPECT_NEAR(origin.getTime(), bestQuarry.getTime(), 0.5);
        EXPECT_EQ(db[bestQuarryIndex].getName(), "BINGHAM MINE");
        auto arrivals = bestQuarry.getArrivals();
        for (const auto &arrival : arrivals)
        {
            EXPECT_NEAR(arrival.getResidual(), 0, 0.9);
        }
        //std::cout << objectiveFunction << " " << bestQuarry.getName() << std::endl;
    }
    origin = createQuarryBlast60020257();
    for (const auto &norm : norms)
    {
        auto [bestQuarryIndex, bestQuarry, bestObjectiveFunction]
            = db.findBestQuarry(origin, norm, p, timeWindow, applyCorrection);
        EXPECT_NEAR(origin.getTime(), bestQuarry.getTime(), 0.5);
        EXPECT_EQ(db.at(bestQuarryIndex).getName(), "BINGHAM MINE");
        auto arrivals = bestQuarry.getArrivals();
        for (const auto &arrival : arrivals)
        {
            EXPECT_NEAR(arrival.getResidual(), 0, 0.9);
        }
    }
    // First event had 13 station phase pairs, second event had 2 new station phase pairs
    EXPECT_EQ(db.getNumberOfStationPhasePairs(), 15);
}

}
