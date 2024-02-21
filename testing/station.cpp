#include <string>
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/position/utahRegion.hpp"
#include <gtest/gtest.h>

using namespace ULocator;

namespace
{

TEST(ULocator, Station)
{
    Position::UtahRegion utah;
    const std::string network{"UU"};
    const std::string name{"FPU"};
    constexpr double latitude{41.02633};
    constexpr double longitude{-111.83683};
    constexpr double elevation{2816};
    Position::WGS84 position{latitude, longitude};
    Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(position, utah);
    station.setElevation(elevation);

    Station sCopy(station);
    EXPECT_EQ(sCopy.getNetwork(), network);
    EXPECT_EQ(sCopy.getName(), name);
    EXPECT_NEAR(sCopy.getElevation(), elevation, 1.e-10);
    auto positionBack = sCopy.getGeographicPosition();
    EXPECT_NEAR(positionBack.getLatitude(),  latitude, 1.e-8);
    EXPECT_NEAR(positionBack.getLongitude(), longitude, 1.e-8);
    auto localCoordinates = sCopy.getLocalCoordinates();
    EXPECT_NEAR(localCoordinates.first,  -28327.81065000718, 1.e-3);
    EXPECT_NEAR(localCoordinates.second,  155642.5444867647, 1.e-3);
}

}
