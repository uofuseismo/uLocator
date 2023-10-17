#include <string>
#include "uLocator/arrival.hpp"
#include <gtest/gtest.h>

using namespace ULocator;

namespace
{

TEST(ULocator, Arrival)
{
    constexpr double time{24};
    constexpr double residual{-1};
    constexpr double standardError{0.12};
    constexpr double azimuth{230};
    constexpr double backAzimuth{50};
    constexpr double distance{3000};
    constexpr int64_t identifier{3838};
    const Arrival::PhaseType phaseType{Arrival::PhaseType::S};
    Arrival arrival;

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
}

}
