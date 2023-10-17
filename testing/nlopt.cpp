#include <vector>
#include <string>
#include <umps/logging/standardOut.hpp>
#include "uLocator/nlopt.hpp"
#include "uLocator/nloptOptions.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/firstArrivalRayTracer.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include <gtest/gtest.h>

#define UTM_ZONE 12

using namespace ULocator;

namespace
{

struct TempArrival 
{
    std::string network;
    std::string station;
    std::string phase;
    Position::WGS84 position;
    double elevation;
    double time;
    double standardError;
};

std::vector<Arrival> createArrivals60484742()
{
/*
60484742,UU,CWU,01,EHZ,P,3855859,1646684879.079617,0.75,0,35,8.5,154.6,-0.14,40.44584,-112.10217,1945.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,MID,01,ENZ,P,3855864,1646684879.3381984,1.0,1,35,9.2,271.5,0.05,40.51733,-112.25467,1722.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,NOQ,01,HHZ,P,3855869,1646684880.4258294,1.0,1,35,15.4,7.9,0.13,40.6525,-112.12033,1622.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,WTU,01,EHZ,P,3855874,1646684880.785694,0.5,0,35,17.6,112.4,0.13,40.45483,-111.9535,1552.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,GMU,01,EHZ,P,3855879,1646684883.5544589,0.5,0,35,33.0,78.1,0.2,40.5755,-111.76317,1829.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,SAIU,02,EHZ,P,3855884,1646684883.9708252,1.0,1,35,37.9,355.4,-0.08,40.85483,-112.1815,1384.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,RBU,02,EHZ,P,3855889,1646684885.032984,0.5,0,35,41.0,43.8,0.37,40.78083,-111.80833,1676.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,SNUT,02,EHZ,P,3855894,1646684886.6219738,0.5,0,35,51.3,323.5,0.22,40.88567,-112.509,1652.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,NAIU,01,EHZ,P,3855899,1646684887.0673397,0.5,0,35,56.1,352.9,-0.09,41.01617,-112.228,1472.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,NLU,01,HHZ,P,3855904,1646684888.7377942,0.5,0,35,62.6,174.5,0.35,39.95483,-112.075,2036.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,MOUT,01,HHZ,P,3855909,1646684891.405937,0.5,0,35,79.3,16.3,0.02,41.19894,-111.87953,2748.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,SPU,01,HHZ,P,3855914,1646684893.4537792,0.5,0,35,91.8,344.0,0.09,41.30867,-112.44917,2086.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
60484742,UU,HVU,01,HHZ,P,3855919,1646684902.3164327,0.75,0,26,150.2,339.7,-0.49,41.77967,-112.775,1609.0,40.5151667,-112.1455,-2.0,1646684876.8099976,1.6,l,I,qb
*/
    auto pcwu  = TempArrival
                 { "UU", "CWU",  "P", Position::WGS84{40.44584,-112.10217,  UTM_ZONE}, 1945.0, 1646684879.079617,   0.06};
    auto pmid  = TempArrival
                 { "UU", "MID",  "P", Position::WGS84{40.51733,-112.25467,  UTM_ZONE}, 1722.0, 1646684879.3381984,  0.03};
    auto pnoq  = TempArrival
                 { "UU", "NOQ",  "P", Position::WGS84{40.6525,-112.12033,   UTM_ZONE}, 1622.0, 1646684880.4258294,  0.03};
    auto pwtu  = TempArrival 
                 { "UU", "WTU",  "P", Position::WGS84{40.45483,-111.9535,   UTM_ZONE}, 1552.0, 1646684880.785694,   0.15};
    auto pgmu  = TempArrival
                 { "UU", "GMU",  "P", Position::WGS84{40.5755,-111.76317,   UTM_ZONE}, 1829.0, 1646684883.5544589,  0.15};
    auto psaiu = TempArrival
                 { "UU", "SAIU", "P", Position::WGS84{40.85483,-112.1815,   UTM_ZONE}, 1384.0, 1646684883.9708252,  0.03};
    auto prbu  = TempArrival
                 { "UU", "RBU",  "P", Position::WGS84{40.78083,-111.80833,  UTM_ZONE}, 1676.0, 1646684885.032984,   0.15};
    auto psnut = TempArrival
                 { "UU", "SNUT", "P", Position::WGS84{40.88567,-112.509,    UTM_ZONE}, 1652.0, 1646684886.6219738,  0.15};
    auto pnaiu = TempArrival
                 { "UU", "NAIU", "P", Position::WGS84{41.01617,-112.228,    UTM_ZONE}, 1472.0, 1646684887.0673397,  0.15};
    auto pnlu  = TempArrival
                 {"UU", "NLU",   "P", Position::WGS84{39.95483,-112.075,    UTM_ZONE}, 2036.0, 1646684888.7377942,  0.15};
    auto pmout = TempArrival
                 {"UU", "MOUT",  "P", Position::WGS84{41.19894,-111.87953,  UTM_ZONE}, 2748.0, 1646684891.405937,   0.15};
    auto pspu  = TempArrival
                 {"UU", "SPU",   "P", Position::WGS84{41.30867,-112.44917,  UTM_ZONE}, 2086.0, 1646684893.4537792,  0.15};
    auto phvu  = TempArrival
                 {"UU", "HVU",   "P", Position::WGS84{41.77967,-112.775,    UTM_ZONE}, 1609.0, 1646684902.3164327,  0.06};
    std::vector<TempArrival> tempArrivals{pcwu, pmid, pnoq, pwtu, pgmu, psaiu, prbu, psnut, pnaiu, pnlu, pmout, pspu, phvu};
    std::vector<Arrival> arrivals;
    for (const auto &tempArrival : tempArrivals)
    {
        Station station;
        station.setNetwork(tempArrival.network);
        station.setName(tempArrival.station);
        station.setGeographicPosition(tempArrival.position);
        station.setElevation(tempArrival.elevation);
        Arrival arrival;
        arrival.setStation(station);
        arrival.setPhase(tempArrival.phase);
        arrival.setTime(tempArrival.time);
        arrival.setStandardError(tempArrival.standardError);
        arrivals.push_back(arrival);
    }
    assert(arrivals.size() == 13);
    return arrivals;
}

std::vector<Arrival> createArrivals60484462()
{
/*
60484462,UU,ASU7,01,HHZ,P,3933784,1646479476.6776488,1.0,-1,161,1.0,185.5,-0.04,38.552595,-112.220772,1894.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,LB,MVU,  ,HHZ,P,3854794,1646479477.4650624,1.0,1,102,6.5,174.3,0.0,38.5037,-112.212303,2239.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,ASU7,01,HHZ,S,3933789,1646479477.5781794,0.75,0,161,1.0,185.5,0.16,38.552595,-112.220772,1894.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,TCRU,01,HHZ,P,3854799,1646479479.779858,1.0,-1,91,20.5,285.2,0.04,38.6095,-112.4472,2293.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,TCRU,01,HHZ,S,3933794,1646479483.044294,0.75,0,91,20.5,285.2,0.31,38.6095,-112.4472,2293.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,WCU,02,EHZ,P,3854804,1646479483.8797174,0.75,0,90,46.2,14.0,-0.06,38.96467,-112.09067,2673.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,ECUT,01,HHZ,P,3933804,1646479487.7271104,0.5,0,90,68.3,6.3,-0.17,39.171697,-112.133201,2136.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,DWU,01,EHZ,P,3933809,1646479490.553661,0.5,0,90,84.7,233.5,-0.21,38.10534,-112.9975,2270.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,NMU,01,EHZ,S,3933799,1646479492.619098,0.75,0,90,55.1,265.0,-0.4,38.5165,-112.85,1853.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,CVRU,01,HHZ,P,3933814,1646479492.699132,0.5,0,90,99.1,66.1,0.31,38.9176,-111.1716,1912.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
60484462,UU,SWUT,01,HHZ,P,3933819,1646479496.371336,0.5,0,90,120.0,315.6,-0.26,39.3286,-113.1954,1644.0,38.5615,-112.2196667,2.03,1646479475.880001,1.14,l,F,eq
*/
    auto pasu7 = TempArrival
                 { "UU", "ASU7", "P", Position::WGS84{38.552595,-112.220772, UTM_ZONE}, 1894.0, 1646479476.6776488,   0.03};
    auto pmvu  = TempArrival
                 { "LB", "MVU" , "P", Position::WGS84{38.5037,  -112.212303, UTM_ZONE}, 2239.0, 1646479477.4650624,   0.03};
    auto sasu7 = TempArrival
                 { "UU", "ASU7", "S", Position::WGS84{38.552595,-112.220772, UTM_ZONE}, 1894.0, 1646479477.5781794,   0.75};
    auto ptcru = TempArrival 
                 { "UU", "TCRU", "P", Position::WGS84{38.6095,  -112.4472,   UTM_ZONE}, 2293.0, 1646479479.779858,    0.03};
    auto stcru = TempArrival
                 { "UU", "TCRU", "S", Position::WGS84{38.6095,  -112.4472,   UTM_ZONE}, 2293.0, 1646479483.044294,    0.06};
    auto pwcu  = TempArrival
                 { "UU", "WCU",  "P", Position::WGS84{38.96467, -112.09067,  UTM_ZONE}, 2673.0, 1646479483.8797174,   0.06};
    auto pecut = TempArrival
                 { "UU", "ECUT", "P", Position::WGS84{39.171697,-112.133201, UTM_ZONE}, 2136.0, 1646479487.7271104,   0.15};
    auto pdwu  = TempArrival
                 { "UU", "DWU",  "P", Position::WGS84{38.10534, -112.9975,   UTM_ZONE}, 2270.0, 1646479490.553661,    0.15};
    auto snmu  = TempArrival
                 { "UU", "NMU",  "S", Position::WGS84{38.5165,  -112.85,     UTM_ZONE}, 1853.0, 1646479492.619098,    0.06};
    auto pcvru = TempArrival
                 {"UU", "CVRU", "P", Position::WGS84{38.9176,  -111.1716,   UTM_ZONE},  1912.0, 1646479492.699132,    0.15};
    auto pswut = TempArrival
                 {"UU", "SWUT", "P", Position::WGS84{39.3286,  -113.1954,   UTM_ZONE},  1644.0, 1646479496.371336,    0.15};
    std::vector<TempArrival> tempArrivals{pasu7, pmvu, sasu7, ptcru, stcru, pwcu, pecut, pdwu, snmu, pcvru, pswut};
    std::vector<Arrival> arrivals;
    for (const auto &tempArrival : tempArrivals)
    {
        Station station;
        station.setNetwork(tempArrival.network);
        station.setName(tempArrival.station);
        station.setGeographicPosition(tempArrival.position);
        station.setElevation(tempArrival.elevation);
        Arrival arrival;
        arrival.setStation(station);
        arrival.setPhase(tempArrival.phase);
        arrival.setTime(tempArrival.time);
        arrival.setStandardError(tempArrival.standardError);
        arrivals.push_back(arrival);
    }
    assert(arrivals.size() == 11);
    return arrivals; 
}

Origin createOrigin60484742()
{
    Origin origin;
    origin.setIdentifier(60484742);
    origin.setTime(1646684876.8099976);
    origin.setDepth(-2000.00);
    origin.setEpicenter( Position::WGS84{40.5151667,-112.1455,  UTM_ZONE} );
    origin.setArrivals(createArrivals60484742());
    return origin;
}

Origin createOrigin60484462()
{
    Origin origin;
    origin.setIdentifier(60484462);
    origin.setTime(1646479475.88);
    origin.setDepth(2000.03);
    origin.setEpicenter( Position::WGS84{38.5615, -112.2196667,  UTM_ZONE} );
    origin.setArrivals(createArrivals60484462());
    return origin;
}

std::unique_ptr<ULocator::TravelTimeCalculatorMap>
    createTravelTimeMap(const std::vector<Arrival> &arrivals)
{
    constexpr bool useAlternateModel = false;
    auto map = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    for (const auto &arrival : arrivals)
    {
        auto calculator = std::make_unique<FirstArrivalRayTracer> ();
        if (arrival.getPhase() == "P")
        {
            calculator->initializeUtahP(arrival.getStation(), useAlternateModel);
        }
        else
        {
            calculator->initializeUtahS(arrival.getStation(), useAlternateModel);
        }
        map->insert(arrival.getStation(), arrival.getPhase(),
                    std::move(calculator));
    }
    return map;
}

}

//----------------------------------------------------------------------------//

/*
TEST(NLOpt, LeastSquaresEarthquake)
{
    constexpr bool north{true};
    std::shared_ptr<UMPS::Logging::ILog> logger
    {
        std::make_shared<UMPS::Logging::StandardOut>
            (UMPS::Logging::Level::Debug)
    };

    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::LeastSquares);
    options.setUTMZone(UTM_ZONE, north);
    options.setDefaultElevation(1500);

    NLOpt solver(logger);
    solver.setOptions(options);

    auto origin = createOrigin60484462();
    auto arrivals = origin.getArrivals();
    auto travelTimeMap = createTravelTimeMap(arrivals);

    solver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
    solver.setArrivals(arrivals);

    solver.locateEarthquake();
    auto newOrigin = solver.getOrigin();
    auto newLoss = solver.evaluateLoss(newOrigin);
    auto oldLoss = solver.evaluateLoss(origin);
    std::cout << "new, old loss " << newLoss << " " << oldLoss << std::endl;
}

TEST(NLOpt, L1Earthquake)
{
    constexpr bool north{true};
    std::shared_ptr<UMPS::Logging::ILog> logger
    {   
        std::make_shared<UMPS::Logging::StandardOut>
            (UMPS::Logging::Level::Info)
    };

    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::L1);
    options.setUTMZone(UTM_ZONE, north);
    options.setDefaultElevation(1500);

    NLOpt solver(logger);
    solver.setOptions(options);

    auto origin = createOrigin60484462();
    auto arrivals = origin.getArrivals();
    auto travelTimeMap = createTravelTimeMap(arrivals);
    
    solver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
    solver.setArrivals(arrivals);
    
    solver.locateEarthquake();
    auto newOrigin = solver.getOrigin();
    auto newLoss = solver.evaluateLoss(newOrigin);
    auto oldLoss = solver.evaluateLoss(origin);
    std::cout << "new, old loss " << newLoss << " " << oldLoss << std::endl;
}

TEST(NLOPt, LeastSquaresBlast)
{
    constexpr bool north{true};
    std::shared_ptr<UMPS::Logging::ILog> logger
    {   
        std::make_shared<UMPS::Logging::StandardOut>
            (UMPS::Logging::Level::Info)
    };
    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::LeastSquares);
    options.setUTMZone(UTM_ZONE, north);
    options.setInitialQuarryBlastSearchDepth(-2000);

    NLOpt solver(logger);
    solver.setOptions(options);

    auto origin = createOrigin60484742();
    auto arrivals = origin.getArrivals();
    auto travelTimeMap = createTravelTimeMap(arrivals);

    solver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
    solver.setArrivals(arrivals);
    solver.locateQuarryBlast();
    auto newOrigin = solver.getOrigin();
    auto newLoss = solver.evaluateLoss(newOrigin);
    auto oldLoss = solver.evaluateLoss(origin);
    std::cout << "new, old loss " << newLoss << " " << oldLoss << std::endl;
}

TEST(NLOPt, L1Blast)
{
    constexpr bool north{true};
    std::shared_ptr<UMPS::Logging::ILog> logger
    {   
        std::make_shared<UMPS::Logging::StandardOut>
            (UMPS::Logging::Level::Info)
    };  
    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::L1);
    options.setUTMZone(UTM_ZONE, north);
    options.setInitialQuarryBlastSearchDepth(-2000);

    NLOpt solver(logger);
    solver.setOptions(options);

    auto origin = createOrigin60484742();
    auto arrivals = origin.getArrivals();
    auto travelTimeMap = createTravelTimeMap(arrivals);

    solver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
    solver.setArrivals(arrivals);
    solver.locateQuarryBlast();
    auto newOrigin = solver.getOrigin();
    auto newLoss = solver.evaluateLoss(newOrigin);
    auto oldLoss = solver.evaluateLoss(origin);
    std::cout << "new, old loss " << newLoss << " " << oldLoss << std::endl;
}

TEST(NLOpt, DoubleDifferenceL1Earthquake)
{
    constexpr bool north{true};
    std::shared_ptr<UMPS::Logging::ILog> logger
    {   
        std::make_shared<UMPS::Logging::StandardOut>
            (UMPS::Logging::Level::Debug)
    };

    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setObjectiveFunction(
        NLOptOptions::ObjectiveFunction::DoubleDifferenceL1);
    options.setUTMZone(UTM_ZONE, north);
    options.setDefaultElevation(1500);

    NLOpt solver(logger);
    solver.setOptions(options);

    auto origin = createOrigin60484462();
    auto arrivals = origin.getArrivals();
    auto travelTimeMap = createTravelTimeMap(arrivals);
    
    solver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
    solver.setArrivals(arrivals);
    
    solver.locateEarthquake();
    auto newOrigin = solver.getOrigin();
    auto newLoss = solver.evaluateLoss(newOrigin);
    auto oldLoss = solver.evaluateLoss(origin);
    std::cout << "new, old loss " << newLoss << " " << oldLoss << std::endl;
}
*/
TEST(NLOpt, LpEarthquake)
{
    constexpr bool north{true};
    std::shared_ptr<UMPS::Logging::ILog> logger
    {
        std::make_shared<UMPS::Logging::StandardOut>
            (UMPS::Logging::Level::Debug)
    };

    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::LP);
    options.setUTMZone(UTM_ZONE, north);
    options.setDefaultElevation(1500);

    NLOpt solver(logger);
    solver.setOptions(options);

    auto origin = createOrigin60484462();
    auto arrivals = origin.getArrivals();
    auto travelTimeMap = createTravelTimeMap(arrivals);

    solver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
    solver.setArrivals(arrivals);

    solver.locateEarthquake();
    auto newOrigin = solver.getOrigin();
    auto newLoss = solver.evaluateLoss(newOrigin);
    auto oldLoss = solver.evaluateLoss(origin);
    std::cout << "new, old loss " << newLoss << " " << oldLoss << std::endl;
}
