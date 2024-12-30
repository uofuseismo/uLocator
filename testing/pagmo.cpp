#include <iostream>
#include <iomanip>
#include <string>
#include <uLocator/origin.hpp>
#include <uLocator/position/wgs84.hpp>
#include <uLocator/position/utahRegion.hpp>
#include <uLocator/station.hpp>
#include <uLocator/arrival.hpp>
#include <uLocator/uussRayTracer.hpp>
#include <uLocator/travelTimeCalculatorMap.hpp>
#include <uLocator/topography/constant.hpp>
#include <uLocator/optimizers/pagmo/particleSwarm.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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

/// Problematic FERVO event
std::vector<ULocator::Arrival> createArrivalsFervo()
{   
    std::vector<ULocator::Arrival> arrivals;
    arrivals.push_back(::toArrival("UU", "FORK", "P", 38.501579, -112.886641, 1685.7,  1733256447.945893, 0.03));
    arrivals.push_back(::toArrival("UU", "NMU",  "P", 38.5177,   -112.8509,   1848.91, 1733256448.292053, 0.03));
    arrivals.push_back(::toArrival("UU", "FORU", "P", 38.458801, -112.861487, 1842.14, 1733256448.604687, 0.03));
    arrivals.push_back(::toArrival("UU", "FORK", "S", 38.501579, -112.886641, 1685.7,  1733256461.514645, 0.06));
    arrivals.push_back(::toArrival("UU", "NMU",  "S", 38.5177,   -112.8509,   1848.91, 1733256462.003129, 0.06));
    arrivals.push_back(::toArrival("UU", "FORU", "S", 38.458801, -112.861487, 1842.14, 1733256462.634224, 0.06));
    return arrivals;
}   

}

TEST_CASE("ULocator::Pagmo::ParticleSwarm", "[general]")
{
    ULocator::Origin referenceOrigin;
    referenceOrigin.setTime(1646479475);
    referenceOrigin.setEpicenter(ULocator::Position::WGS84 {38.562, -112.220} );
    referenceOrigin.setDepth(2000);

    auto arrivals = createArrivals();
    auto calculatorMap = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    for (const auto &arrival : arrivals)
    {
        std::unique_ptr<ULocator::ITravelTimeCalculator> calculator;
        if (arrival.getPhase() == "P")
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::P, 
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        else
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::S,
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        calculatorMap->insert(arrival.getStation(), arrival.getPhase(),
                              std::move(calculator));
    }
    auto topography = std::make_unique<ULocator::Topography::Constant> ();
    topography->set(4000);
 
    ULocator::Optimizers::Pagmo::ParticleSwarm pso;
    pso.setGeographicRegion(ULocator::Position::UtahRegion {});
    pso.setTravelTimeCalculatorMap(std::move(calculatorMap));
    pso.setTopography(std::move(topography));
    pso.setArrivals(arrivals);

    ULocator::Origin noGuess;
    pso.setNumberOfParticles(10);
    pso.setNumberOfGenerations(50);
    pso.locate(noGuess,
               ULocator::Optimizers::IOptimizer::LocationProblem::ThreeDimensionsAndTime,
               ULocator::Optimizers::IOptimizer::Norm::Lp);
    auto origin = pso.getOrigin();
    std::cout << std::setprecision(12) << origin.getEpicenter().getLatitude() << " " << origin.getEpicenter().getLongitude() << " " << " " << origin.getDepth() << " " << origin.getTime() << std::endl;

/*
   pso.setArrivals(createArrivalsFervo());
   pso.locate(noGuess,
              ULocator::Optimizers::IOptimizer::LocationProblem::ThreeDimensionsAndTime,
              ULocator::Optimizers::IOptimizer::Norm::Lp); 
    origin = pso.getOrigin();
    std::cout << std::setprecision(12) << origin.getEpicenter().getLatitude() << " " << origin.getEpicenter().getLongitude() << " " << " " << origin.getDepth() << " " << origin.getTime() << std::endl; 
*/
}

/*
TEST_CASE("ULocator::Pagmo::ParticleSwarm", "[fixedDepth]")
{

}

TEST_CASE("ULocator::Pagmo::ParticleSwarm", "[blast"])
{

}
*/
