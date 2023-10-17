#include <string>
#include <cmath>
#include <vector>
//#include "uLocator/hamiltonianMonteCarlo.hpp"
#include "uLocator/nloptOptions.hpp"
#include "uLocator/nlopt.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/station.hpp"
#include "uLocator/travelTimeCalculator.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/wgs84.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace ULocator;

class HomogeneousTravelTime : public ITravelTimeCalculator
{
public:
    ~HomogeneousTravelTime() override = default;
    HomogeneousTravelTime(const Station &station,
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
        mStationEasting  = mStation.getGeographicPosition().getEasting();
        mStationNorthing = mStation.getGeographicPosition().getNorthing(); 
        mElevation = mStation.getElevation();
    }
    double evaluate(const Position::WGS84 &epicenter, const double depth, bool) const override
    {
        auto dx = epicenter.getEasting()  - mStationEasting;
        auto dy = epicenter.getNorthing() - mStationNorthing;
        auto dz = depth - (-mElevation);
        auto distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        //std::cout << dx*1.e-3 << " " << dy*1e-3 << std::endl;
        return distance/mVelocity;
    }
    double evaluate(const Position::WGS84 &epicenter, const double depth,
                    double *dtdx, double *dtdy, double *dtdz,
                    bool) const override
    {
        auto dx = epicenter.getEasting()  - mStationEasting;
        auto dy = epicenter.getNorthing() - mStationNorthing;
        auto dz = depth - (-mElevation);
        auto distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        *dtdx = dx/(distance*mVelocity);
        *dtdy = dy/(distance*mVelocity);
        *dtdz = dz/(distance*mVelocity);
        return distance/mVelocity;
    }
    std::array<double, 3> computeGradient(const Position::WGS84 &epicenter,
                                          const double depth) const
    {
        auto dx = epicenter.getEasting()  - mStationEasting;
        auto dy = epicenter.getNorthing() - mStationNorthing;
        auto dz = depth - (-mElevation);
        auto distance = std::sqrt( dx*dx + dy*dy + dz*dz );
        std::array<double, 3> gradient;
        gradient[0] = dx/(distance*mVelocity);
        gradient[1] = dy/(distance*mVelocity);
        gradient[2] = dz/(distance*mVelocity);
        return gradient;
    }

    Station mStation;
    std::string mPhase;
    double mStationEasting;
    double mStationNorthing;
    double mVelocity{5000};
    double mElevation{0};
};

Station toStation(const std::string &network,
                  const std::string &name,
                  const double latitude, const double longitude,
                  const int utmZone = 12)
{
    Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        Position::WGS84 {latitude, longitude, utmZone});
    station.setElevation(0);
    return station;
}
  

std::vector<Station> createStations()
{
    std::vector<Station> stations;
    stations.push_back(::toStation("UU", "HVU",  41.7797, -112.775));
    stations.push_back(::toStation("UU", "NPI",  42.1473, -112.518));
    stations.push_back(::toStation("UU", "PTU",  41.9293, -112.325));
    stations.push_back(::toStation("UU", "MTUT", 41.6709, -112.455));
    stations.push_back(::toStation("UU", "MLI",  42.0268, -112.126));
    stations.push_back(::toStation("UU", "LTU",  41.5918, -112.247));
    stations.push_back(::toStation("UU", "EPU",  41.3915, -112.409));
    stations.push_back(::toStation("UU", "HONU", 41.615,  -112.05));
    stations.push_back(::toStation("UU", "WVUT", 41.6102, -111.959));
    stations.push_back(::toStation("UU", "GZU",  41.4252, -111.975));
    stations.push_back(::toStation("IE", "PTI",  42.8703, -112.37));
    stations.push_back(::toStation("US", "HWUT", 41.6069, -111.565));
    stations.push_back(::toStation("UU", "NAIU", 41.0162, -112.228));
    stations.push_back(::toStation("UU", "SNUT", 40.8857, -112.509));
    stations.push_back(::toStation("UU", "MCU",  41.4619, -111.509));
    stations.push_back(::toStation("UU", "BMUT", 41.9582, -111.234));
    stations.push_back(::toStation("UU", "SAIU", 40.8548, -112.181));
    return stations;
}

TEST(ULocator, Homogeneous)
{
    const double eventLatitude{41.9576667};
    const double eventLongitude{-112.8155};
    const double eventDepth{5000};//8300};
    const double originTime{10};
    auto stations = ::createStations();
    // Create the travel time tables and observations
    std::vector<Arrival> arrivals;
    Position::WGS84 epicenter{eventLatitude, eventLongitude, 12};
    auto travelTimeTableMap = std::make_unique<TravelTimeCalculatorMap> ();
    int identifier = 100;
    for (const auto &station : stations)
    {
        std::unique_ptr<const ITravelTimeCalculator> pTable{std::make_unique<::HomogeneousTravelTime> (station, "P")};
        std::unique_ptr<const ITravelTimeCalculator> sTable{std::make_unique<::HomogeneousTravelTime> (station, "S")};
        auto pTime = originTime + pTable->evaluate(epicenter, eventDepth);
        auto sTime = originTime + sTable->evaluate(epicenter, eventDepth);
        travelTimeTableMap->insert(station, "P", std::move(pTable));
        travelTimeTableMap->insert(station, "S", std::move(sTable));

        Arrival arrival;
        arrival.setTime(pTime);
        arrival.setPhase("P");
        arrival.setStandardError(0.1);
        arrival.setIdentifier(identifier);
        arrival.setStation(station);
        arrivals.push_back(arrival);
        identifier = identifier + 1;
    }
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
    EXPECT_NEAR(eventDepth,     origin.getDepth(),                    1.e-3);
    EXPECT_NEAR(originTime,     origin.getTime(),                     1.e-2);
    //std::cout << "lat/lon residual: " << eventLatitude - origin.getEpicenter().getLatitude() << "," << eventLongitude - origin.getEpicenter().getLongitude() << std::endl;
    //std::cout << "depth residual: " << eventDepth - origin.getDepth() << std::endl;
    //std::cout << "ot residual: " << originTime - origin.getTime() << std::endl;
}

}
