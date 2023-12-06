#include <string>
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/ynp.hpp"
#include <gtest/gtest.h>

using namespace ULocator;

namespace
{

TEST(ULocator, Arrival)
{
    const std::string network{"WY"};
    const std::string stationName{"YML"};
    constexpr double stationLatitude{44.60533};
    constexpr double stationLongitude{-110.64317};
    constexpr double stationElevation{2653.0};
    Position::YNP ynp;
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
    EXPECT_NEAR(copy.getTime(),     time, 1.e-14);
    EXPECT_NEAR(copy.getStandardError(), standardError, 1.e-14);
    EXPECT_NEAR(copy.getResidual(), residual, 1.e-14);
    EXPECT_NEAR(copy.getDistance(), distance, 1.e-12);
    EXPECT_NEAR(copy.getAzimuth(), azimuth, 1.e-14);
    EXPECT_NEAR(copy.getBackAzimuth(), backAzimuth, 1.e-14); 
    EXPECT_EQ(copy.getIdentifier(), identifier);
    EXPECT_EQ(copy.getPhase(), "S");
    EXPECT_EQ(copy.getStationReference().getNetwork(), network);
    EXPECT_TRUE(copy.haveStation());
    EXPECT_NEAR(copy.getStationReference().getElevation(),
                stationElevation, 1.e-8);
    EXPECT_NEAR(copy.getStationReference().getLocalCoordinates().first,
                8480.827504103971, 1.e-4);
    EXPECT_NEAR(copy.getStationReference().getLocalCoordinates().second,
                11710.12450692187, 1.e-4);
    EXPECT_EQ(copy.getStation().getName(), stationName);

    arrival.clear();
    EXPECT_FALSE(arrival.haveTime());
}

}
