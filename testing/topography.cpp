#include <vector>
#include "uLocator/topography/constant.hpp"
#include "uLocator/topography/gridded.hpp"
#include "uLocator/position/utah.hpp"
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
    EXPECT_NEAR(topography.getMinimumAndMaximumElevation().first,
                elevation, 1.e-10);
    EXPECT_NEAR(topography.getMinimumAndMaximumElevation().second,
                elevation, 1.e-10);
    double dElevationDx, dElevationDy;
    EXPECT_NEAR(topography(41, -111, &dElevationDx, &dElevationDy),
                elevation, 1.e-10);
    EXPECT_NEAR(dElevationDx, 0, 1.e-10);
    EXPECT_NEAR(dElevationDy, 0, 1.e-10); 
    topography.clear();
    EXPECT_FALSE(topography.haveTopography());
}

TEST(ULocatorTopography, Gridded)
{
    std::vector<double> xs{5, 6, 7, 8};
    std::vector<double> ys{2, 3, 4};
    std::vector<double> elevations(xs.size()*ys.size());
    auto nx = static_cast<int> (xs.size());
    auto ny = static_cast<int> (ys.size());
    double dElevDx = 2;
    double dElevDy =-3;
    for (int iy = 0; iy < ny; ++iy)
    {
        for (int ix = 0; ix < nx; ++ix)
        {
            auto indx = iy*nx + ix;
            auto elevation = dElevDx*xs[ix] + dElevDy*ys[iy];
            elevations.at(indx) = elevation;
        }
    }
    Gridded topography;
    EXPECT_NO_THROW(topography.set(xs.size(), xs.data(),
                                   ys.size(), ys.data(),
                                   elevations.size(), elevations.data()));
    EXPECT_NEAR(topography.getMinimumAndMaximumElevation().first,
                *std::min_element(elevations.begin(), elevations.end()),
                1.e-10);
    EXPECT_NEAR(topography.getMinimumAndMaximumElevation().second,
                *std::max_element(elevations.begin(), elevations.end()),
                1.e-10);
    EXPECT_TRUE(topography.haveTopography());
    constexpr double tol{1.e-7};
    // Corners
    // 8  9  10 11
    // 4  5   6  7
    // 0  1   2  3
    EXPECT_NEAR(elevations.at(0),  topography.evaluate(5, 2), tol);
    EXPECT_NEAR(elevations.at(11), topography.evaluate(8, 4), tol);
    EXPECT_NEAR(elevations.at(3),  topography.evaluate(8, 2), tol);
    EXPECT_NEAR(elevations.at(8),  topography.evaluate(5, 4), tol);
    // Out of bounds 
    EXPECT_NEAR(elevations.at(0),  topography.evaluate(5 - 1, 2 - 1), tol);
    EXPECT_NEAR(elevations.at(11), topography.evaluate(8 + 1, 4 + 1), tol);
    EXPECT_NEAR(elevations.at(3),  topography.evaluate(8 + 1, 2 - 1), tol);
    EXPECT_NEAR(elevations.at(8),  topography.evaluate(5 - 1, 4 + 1), tol);
    // Along bottom, top, left, right
    EXPECT_NEAR(elevations.at(1),  topography.evaluate(6, 2 - 1), tol);
    EXPECT_NEAR(elevations.at(9),  topography.evaluate(6, 4 + 1), tol);
    EXPECT_NEAR(elevations.at(4),  topography.evaluate(5 - 1, 3), tol);
    EXPECT_NEAR(elevations.at(7),  topography.evaluate(8 + 1, 3), tol);

    // Inside model
    EXPECT_NEAR(elevations.at(5),  topography.evaluate(6, 3), tol);
    EXPECT_NEAR(dElevDx*6.8 + dElevDy*3.3,
                topography.evaluate(6.8, 3.3), tol);
    EXPECT_NEAR(dElevDx*6.4 + dElevDy*2.1,
                topography.evaluate(6.4, 2.1), tol);
    EXPECT_NEAR(dElevDx*6.4 + dElevDy*2,
                topography.evaluate(6.4, 2), tol);
    EXPECT_NEAR(dElevDx*7.5 + dElevDy*3.5,
                topography.evaluate(7.5, 3.5), tol);
    EXPECT_NEAR(dElevDx*5.1 + dElevDy*3.5,
                topography.evaluate(5.1, 3.5), tol);
    double dEdx, dEdy;
    EXPECT_NEAR(dElevDx*6.4 + dElevDy*2.1,
                topography.evaluate(6.4, 2.1, &dEdx, &dEdy), tol);
    EXPECT_NEAR(dElevDx, dEdx, tol);
    EXPECT_NEAR(dElevDy, dEdy, tol);
}

/*
TEST(ULocatorTopography, LoadHDF5)
{
    ULocator::Position::Utah utah;
    Gridded topography;
    topography.load("utahTopo.h5", utah);
}
*/

}
