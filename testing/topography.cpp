#include <vector>
#include "uLocator/topography/constant.hpp"
#include "uLocator/topography/gridded.hpp"
#include <gtest/gtest.h>

using namespace ULocator::Topography;

namespace
{

TEST(ULocatorTopography, Constant)
{
    constexpr double elevation{2500};
    Constant topography;
    topography.set(elevation);
    EXPECT_NEAR(topography(41, -111), elevation, 1.e-10);
    topography.clear();
    EXPECT_FALSE(topography.haveTopography());
}

TEST(ULocatorTopography, Gridded)
{
    std::vector<double> latitudes{2, 3, 4};
    std::vector<double> longitudes{5, 6, 7, 8};
    std::vector<double> elevations(latitudes.size()*longitudes.size());
    auto nLatitudes = static_cast<int> (latitudes.size());
    auto nLongitudes = static_cast<int> (longitudes.size());
    double dElevDlat = 2;
    double dElevDlon =-3;
    for (int iLat = 0; iLat < nLatitudes; ++iLat)
    {
        for (int iLon = 0; iLon < nLongitudes; ++iLon)
        {
            auto indx = iLat*nLongitudes + iLon;
            auto elev = dElevDlat*latitudes[iLat] + dElevDlon*longitudes[iLon];
            elevations.at(indx) = elev;
        }
    }
    Gridded topography;
    EXPECT_NO_THROW(topography.set(latitudes.size(), latitudes.data(),
                                   longitudes.size(), longitudes.data(),
                                   elevations.size(), elevations.data()));
    EXPECT_TRUE(topography.haveTopography());
    constexpr double tol{1.e-7};
    // Corners
    // 8  9  10 11
    // 4  5   6  7
    // 0  1   2  3
    EXPECT_NEAR(elevations.at(0),  topography.evaluate(2, 5), tol);
    EXPECT_NEAR(elevations.at(11), topography.evaluate(4, 8), tol);
    EXPECT_NEAR(elevations.at(3),  topography.evaluate(2, 8), tol);
    EXPECT_NEAR(elevations.at(8),  topography.evaluate(4, 5), tol);
    // Out of bounds 
    EXPECT_NEAR(elevations.at(0),  topography.evaluate(2 - 1, 5 - 1), tol);
    EXPECT_NEAR(elevations.at(11), topography.evaluate(4 + 1, 8 + 1), tol);
    EXPECT_NEAR(elevations.at(3),  topography.evaluate(2 - 1, 8 + 1), tol);
    EXPECT_NEAR(elevations.at(8),  topography.evaluate(4 + 1, 5 - 1), tol);
    // Along bottom, top, left, right
    EXPECT_NEAR(elevations.at(1),  topography.evaluate(2 - 1, 6), tol);
    EXPECT_NEAR(elevations.at(9),  topography.evaluate(4 + 1, 6), tol);
    EXPECT_NEAR(elevations.at(4),  topography.evaluate(3, 5 - 1), tol);
    EXPECT_NEAR(elevations.at(7),  topography.evaluate(3, 8 + 1), tol);

    // Inside model
    EXPECT_NEAR(elevations.at(5),  topography.evaluate(3, 6), tol);
    EXPECT_NEAR(dElevDlat*3.3 + dElevDlon*6.8,
                topography.evaluate(3.3, 6.8), tol);
    EXPECT_NEAR(dElevDlat*2.1 + dElevDlon*6.4,
                topography.evaluate(2.1, 6.4), tol);
    EXPECT_NEAR(dElevDlat*2 + dElevDlon*6.4,
                topography.evaluate(2, 6.4), tol);
    EXPECT_NEAR(dElevDlat*3.5 + dElevDlon*7.5,
                topography.evaluate(3.5, 7.5), tol);
    EXPECT_NEAR(dElevDlat*3.5 + dElevDlon*5.1,
                topography.evaluate(3.5, 5.1), tol);
}

}
