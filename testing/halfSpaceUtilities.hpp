#ifndef HALF_SPACE_UTILITIES_HPP
#define HALF_SPACE_UTILITIES_HPP
#include <map>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/utah.hpp"
#include "uLocator/travelTimeCalculator.hpp"
#include "uLocator/topography/constant.hpp"

namespace
{

#define P_VELOCITY 5000
#define S_VELOCITY 3535.533905932738

class QuadraticUtahTopography : public ULocator::Topography::ITopography
{
public:
    [[nodiscard]] double lagrange1d(int node,
                                    const double xi,
                                    const double xi1,
                                    const double xi2,
                                    const double xi3,
                                    double *dfdxi) const
    {
        if (node == 1)
        {
            double denominator = (xi1 - xi2)*(xi1 - xi3);
            if (dfdxi != nullptr){*dfdxi = (2*xi - xi2 - xi3)/denominator;}
            return (xi - xi2)*(xi - xi3)/denominator;
        }
        else if (node == 2)
        { 
            double denominator = (xi2 - xi1)*(xi2 - xi3);
            if (dfdxi != nullptr){*dfdxi = (2*xi - xi1  - xi3)/denominator;}
            return (xi - xi1)*(xi - xi3)/denominator;
        }
        else if (node == 3)
        {
            double denominator = (xi3 - xi1)*(xi3 - xi2);
            if (dfdxi != nullptr){*dfdxi = (2*xi - xi1 - xi2)/denominator;}
            return (xi - xi1)*(xi - xi2)/denominator;
        }
#ifndef NDEBUG
        assert(false);
#endif
    }
    [[nodiscard]] double evaluate(const double x, const double y,
                                  double *dEdx, double *dEdy) const override
    {
        double elevation{0};
        double dedx{0};
        double dedy{0};
        for (int i = 0; i < 9; ++i)
        {
            double dldx, dldy; // No Jacobian necessary
            double fxi  = lagrange1d(mNodes[i].first,  x,
                                     mXi[0],  mXi[1],  mXi[2], &dldx);
            double feta = lagrange1d(mNodes[i].second, y,
                                     mEta[0], mEta[1], mEta[2], &dldy);
            elevation = elevation + mElevations[i]*(fxi*feta);
            // dEdx = dE/dxi dxi/dx + dE/deta deta/dx = dE/dxi
            dedx = dedx + mElevations[i]*feta*dldx;
            dedy = dedy + mElevations[i]*fxi*dldy;
        }
        if (dEdx != nullptr){*dEdx = dedx;}
        if (dEdy != nullptr){*dEdy = dedy;}
        return elevation;
    }
    bool haveTopography() const noexcept override 
    {
        return true;
    }
    std::pair<double, double> getMinimumAndMaximumElevation() const override
    {
        return {0, 4000};
    }
    virtual ~QuadraticUtahTopography() = default;
    double evaluate(const double x, const double y) const override
    {
         return evaluate(x, y, nullptr, nullptr); 
    }
    ULocator::Position::Utah mRegion;
    std::array<double, 9> mElevations{0, 0, 0, 0, 0, 0, 0, 0, 4000};
    std::array<std::pair<int, int>, 9> mNodes{ std::pair{1, 1},
                                               std::pair{3, 1},
                                               std::pair{3, 3},
                                               std::pair{1, 3},
                                               std::pair{2, 1},
                                               std::pair{3, 2},
                                               std::pair{2, 3},
                                               std::pair{1, 2},
                                               std::pair{2, 2} }; 
    std::array<double, 3> mXi{mRegion.getExtentInX().first,
                              0.5*(mRegion.getExtentInX().first + mRegion.getExtentInX().second),
                              mRegion.getExtentInX().second};
    std::array<double, 3> mEta{mRegion.getExtentInY().first,
                               0.5*(mRegion.getExtentInY().first + mRegion.getExtentInY().second),
                               mRegion.getExtentInY().second};

};

class HalfSpaceTravelTime : public ULocator::ITravelTimeCalculator
{
public:
    ~HalfSpaceTravelTime() override = default;
    HalfSpaceTravelTime(const ULocator::Station &station,
                        const std::string &phase) :
         mStation(station), 
         mPhase(phase)
    {   
        if (phase == "P")
        {
            mVelocity = P_VELOCITY;
        }
        else if (phase == "S")
        {
            mVelocity = S_VELOCITY;
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
        if (dtdt0 != nullptr){*dtdt0 = 1;}
        if (dtdx != nullptr){*dtdx =-dx/(distance*mVelocity);}
        if (dtdy != nullptr){*dtdy =-dy/(distance*mVelocity);}
        if (dtdz != nullptr){*dtdz =-dz/(distance*mVelocity);}
        return originTime + distance/mVelocity;
    }   
    ULocator::Station mStation;
    std::string mPhase;
    double mStationX;
    double mStationY;
    double mStationZ;
    double mVelocity{5000};
};

ULocator::Station toStation(const std::string &network,
                            const std::string &name,
                            const double latitude,
                            const double longitude,
                            const double elevation,
                            const int utmZone = 12) 
{
    ULocator::Position::Utah utah;
    ULocator::Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        ULocator::Position::WGS84 {latitude, longitude, utmZone}, utah);
    station.setElevation(elevation);
    return station;
}


std::vector<ULocator::Station> createStations()
{
    std::vector<ULocator::Station> stations;
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

template<typename T>
void fillObjectiveFunction(T &objectiveFunction)
{

/*
{
::QuadraticUtahTopography topo;
int nx = 10;
int ny = 20;
//std::ofstream outfl("elev.txt");
for (int j = 0; j < ny; ++j)
{
    for (int i = 0; i < nx; ++i)
    {
        ULocator::Position::Utah u;
        auto x0 = u.getExtentInX().first; 
        auto x1 = u.getExtentInX().second;
        auto dx = (x1 - x0)/(nx - 1);
        auto y0 = u.getExtentInY().first;
        auto y1 = u.getExtentInY().second;
        auto dy = (y1 - y0)/(ny - 1);
        auto x = x0 + i*dx;
        auto y = y0 + j*dy;
        double dedx, dedy;
        auto elevation = topo.evaluate(x, y, &dedx, &dedy);
        //std::cout << dedy << " " << (topo.evaluate(x, y + 0.01) - elevation)/0.01 << std::endl;
        //outfl << x << " " << y << " "  << topo.evaluate(x, y) << std::endl; 
    }
    //outfl << std::endl;
}
//outfl.close();
}
*/
    ULocator::Position::Utah utah;
    const double eventLatitude{41.9576667};
    const double eventLongitude{-112.8155};
    const double eventDepth{5000};
    const double originTime{10};
    const double penalty{2};
    auto stations = ::createStations();
    // Create the topography
    //constexpr double elevation{2500};
    //auto topography = std::make_unique<ULocator::Topography::Constant> (); 
    //topography->set(elevation);
    auto topography = std::make_unique<::QuadraticUtahTopography> ();
    // Create the travel time tables and observations
    std::vector<ULocator::Arrival> arrivals;
    ULocator::Position::WGS84 epicenter{eventLatitude, eventLongitude, 12};
    auto [xSource, ySource] = utah.geographicToLocalCoordinates(eventLatitude, eventLongitude);
    auto travelTimeTableMap = std::make_unique<ULocator::TravelTimeCalculatorMap> (); 
    std::vector<std::pair<ULocator::Station, std::string>> stationPhases;
    std::vector<double> weights;
    std::vector<double> observations;
    for (const auto &station : stations)
    {   
        std::unique_ptr<const ULocator::ITravelTimeCalculator>
            pTable{std::make_unique<::HalfSpaceTravelTime> (station, "P")};
        std::unique_ptr<const ULocator::ITravelTimeCalculator>
            sTable{std::make_unique<::HalfSpaceTravelTime> (station, "S")};
        auto pTime = pTable->evaluate(originTime, xSource, ySource, eventDepth, false);
        auto sTime = sTable->evaluate(originTime, xSource, ySource, eventDepth, false);
        travelTimeTableMap->insert(station, "P", std::move(pTable));
        travelTimeTableMap->insert(station, "S", std::move(sTable));

        ULocator::Arrival pArrival;
        pArrival.setTime(pTime);
        pArrival.setPhase("P");
        pArrival.setStandardError(0.1);
        pArrival.setIdentifier(arrivals.size() + 1); 
        pArrival.setStation(station);
        arrivals.push_back(pArrival);

        stationPhases.push_back(std::pair {station, "P"});
        observations.push_back(pTime);
        weights.push_back(1./0.1);

        ULocator::Arrival sArrival;
        sArrival.setTime(sTime);
        sArrival.setPhase("S");
        sArrival.setStandardError(0.2);
        sArrival.setIdentifier(arrivals.size() + 1); 
        sArrival.setStation(station);
        if (arrivals.size()%2 == 0)
        {
            arrivals.push_back(sArrival);
            stationPhases.push_back(std::pair {station, "S"});
            weights.push_back(1./0.2); 
            observations.push_back(sTime);
        }
    }
    double minArrivalTime = *std::min_element(observations.begin(),
                                              observations.end());
    objectiveFunction.mTravelTimeCalculatorMap = travelTimeTableMap.release();
    double maxElevation = topography->getMinimumAndMaximumElevation().second;
    //objectiveFunction.mTopography = topography.release();
    objectiveFunction.mStationPhases = stationPhases;
    objectiveFunction.mObservations = observations;
    objectiveFunction.mWeights = weights;
    std::vector<double> lowerBounds{minArrivalTime - 200,
                                    utah.getExtentInX().first,
                                    utah.getExtentInY().first,
                                    -maxElevation};
    std::vector<double> upperBounds{minArrivalTime,
                                    utah.getExtentInX().second,
                                    utah.getExtentInY().second,
                                    65000}; 
    objectiveFunction.mLowerBounds = lowerBounds; 
    objectiveFunction.mUpperBounds = upperBounds;
    objectiveFunction.mPenaltyCoefficient = penalty;
}

}
#endif
