#include <string>
#include <cmath>
#include <vector>
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/station.hpp"
#include "uLocator/optimizers/pagmo/particleSwarm.hpp"
#include "uLocator/topography/constant.hpp"
#include "uLocator/travelTimeCalculator.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/utah.hpp"
#include "uLocator/position/wgs84.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace ULocator;

class HalfSpaceTravelTime : public ITravelTimeCalculator
{
public:
    ~HalfSpaceTravelTime() override = default;
    HalfSpaceTravelTime(const Station &station,
                        const std::string &phase) :
         mStation(station), 
         mPhase(phase)
    {
        if (phase == "P")
        {
            mVelocity = 6000;
        }
        else if (phase == "S")
        {
            mVelocity = 6000*std::sqrt(0.5); // Poisson solid
        }
        else
        {
            throw std::invalid_argument("Unhandled phase");
        }
        mStationX = mStation.getLocalCoordinates().first;
        mStationY = mStation.getLocalCoordinates().second;
        mStationZ =-mStation.getElevation();
    }
    double evaluate(const double originTime,
                    const double xSource,
                    const double ySource,
                    const double zSource,
                    const bool) const override
    {
        auto dx = mStationX - xSource;
        auto dy = mStationY - ySource;
        auto dz = mStationZ - zSource;
        auto distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        //std::cout << dx*1.e-3 << " " << dy*1e-3 << std::endl;
        return originTime + distance/mVelocity;
    }
    double evaluate(const double originTime,
                    const double xSource,
                    const double ySource,
                    const double zSource,
                    double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
                    const bool) const override
    {
        auto dx = mStationX - xSource;
        auto dy = mStationY - ySource;
        auto dz = mStationZ - zSource;
        auto distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        *dtdt0 = 1;
        *dtdx =-dx/(distance*mVelocity);
        *dtdy =-dy/(distance*mVelocity);
        *dtdz =-dz/(distance*mVelocity);
        return distance/mVelocity;
    }
    Station mStation;
    std::string mPhase;
    double mStationX;
    double mStationY;
    double mStationZ;
    double mVelocity{5000};
};

Station toStation(const std::string &network,
                  const std::string &name,
                  const double latitude,
                  const double longitude,
                  const double elevation,
                  const int utmZone = 12)
{
    Position::Utah utah;
    Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        Position::WGS84 {latitude, longitude, utmZone}, utah);
    station.setElevation(elevation);
    return station;
}
  

std::vector<Station> createStations()
{
    std::vector<Station> stations;
    stations.push_back(::toStation("UU", "HVU",  41.7797, -112.775, 1609));
    stations.push_back(::toStation("UU", "NPI",  42.1473, -112.518, 1640));
    stations.push_back(::toStation("UU", "PTU",  41.9293, -112.325, 2192));
    stations.push_back(::toStation("UU", "MTUT", 41.6709, -112.455, 1373));
    stations.push_back(::toStation("UU", "MLI",  42.0268, -112.126, 1896));
    stations.push_back(::toStation("UU", "LTU",  41.5918, -112.247, 1585));
    stations.push_back(::toStation("UU", "EPU",  41.3915, -112.409, 1436));
    stations.push_back(::toStation("UU", "HONU", 41.615,  -112.05,  1515));
    stations.push_back(::toStation("UU", "WVUT", 41.6102, -111.959, 1828));
    stations.push_back(::toStation("UU", "GZU",  41.4252, -111.975, 2648));
    stations.push_back(::toStation("IE", "PTI",  42.8703, -112.37,  1670));
    stations.push_back(::toStation("US", "HWUT", 41.6069, -111.565, 1830));
    stations.push_back(::toStation("UU", "NAIU", 41.0162, -112.228, 1472));
    stations.push_back(::toStation("UU", "SNUT", 40.8857, -112.509, 1652));
    stations.push_back(::toStation("UU", "MCU",  41.4619, -111.509, 2668));
    stations.push_back(::toStation("UU", "BMUT", 41.9582, -111.234, 2243));
    stations.push_back(::toStation("UU", "SAIU", 40.8548, -112.181, 1384));
    return stations;
}

TEST(ULocator, HalfSpace)
{
    Position::Utah utah;
    const double eventLatitude{41.9576667};
    const double eventLongitude{-112.8155};
    const double eventDepth{5000};
    const double originTime{10};
    auto stations = ::createStations();
    // Create the topography
    auto topography = std::make_unique<ULocator::Topography::Constant> ();
    topography->set(-2500);
    // Create the travel time tables and observations
    std::vector<Arrival> arrivals;
    Position::WGS84 epicenter{eventLatitude, eventLongitude, 12};
    auto [xSource, ySource] = utah.geographicToLocalCoordinates(eventLatitude, eventLongitude);
    auto travelTimeTableMap = std::make_unique<TravelTimeCalculatorMap> ();
    for (const auto &station : stations)
    {
        std::unique_ptr<const ITravelTimeCalculator> pTable{std::make_unique<::HalfSpaceTravelTime> (station, "P")};
        std::unique_ptr<const ITravelTimeCalculator> sTable{std::make_unique<::HalfSpaceTravelTime> (station, "S")};
        auto pTime = pTable->evaluate(originTime, xSource, ySource, eventDepth, false);
        auto sTime = sTable->evaluate(originTime, xSource, ySource, eventDepth, false);
        travelTimeTableMap->insert(station, "P", std::move(pTable));
        travelTimeTableMap->insert(station, "S", std::move(sTable));

        Arrival pArrival;
        pArrival.setTime(pTime);
        pArrival.setPhase("P");
        pArrival.setStandardError(0.1);
        pArrival.setIdentifier(arrivals.size() + 1);
        pArrival.setStation(station);
        arrivals.push_back(pArrival);

        Arrival sArrival;
        sArrival.setTime(sTime);
        sArrival.setPhase("S");
        sArrival.setStandardError(0.2);
        sArrival.setIdentifier(arrivals.size() + 1);
        sArrival.setStation(station);
        if (arrivals.size()%2 == 0){arrivals.push_back(sArrival);}
    }

    ULocator::Optimizers::Pagmo::ParticleSwarm pso;
    pso.setGeographicRegion(utah);
    pso.setTravelTimeCalculatorMap(std::move(travelTimeTableMap));
    pso.setTopography(std::move(topography));
    pso.setArrivals(arrivals);
    pso.locate();
/*
    // Create the locator 
    NLOptOptions options(NLOptOptions::Region::Utah);
    options.setInitialAbsoluteModelTolerance(1.e-4);
    options.setAbsoluteModelTolerance(1.e-7); 
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::LeastSquares);
    NLOpt locator; 
    locator.setOptions(options);
    EXPECT_NO_THROW(locator.setTravelTimeCalculatorMap(std::move(travelTimeTableMap)));
    EXPECT_NO_THROW(locator.setArrivals(arrivals));
    EXPECT_NO_THROW(locator.locateEarthquake());
    auto origin = locator.getOrigin();
    EXPECT_NEAR(eventLatitude,  origin.getEpicenter().getLatitude(),  1.e-3);
    EXPECT_NEAR(eventLongitude, origin.getEpicenter().getLongitude(), 1.e-3);
    EXPECT_NEAR(eventDepth,     origin.getDepth(),                    1);
    EXPECT_NEAR(originTime,     origin.getTime(),                     1.e-2);

    //std::cout << "lat/lon residual: " << eventLatitude - origin.getEpicenter().getLatitude() << "," << eventLongitude - origin.getEpicenter().getLongitude() << std::endl;
    //std::cout << "depth residual: " << eventDepth - origin.getDepth() << std::endl;
    //std::cout << "ot residual: " << originTime - origin.getTime() << std::endl;
*/
}

}
