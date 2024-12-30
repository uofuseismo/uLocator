#include <vector>
#include "uLocator/topography/constant.hpp"
#include "uLocator/topography/gridded.hpp"
#include "uLocator/position/utahRegion.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace ULocator::Topography;


TEST_CASE("ULocator::Topography", "[constant]")
{
    constexpr double elevation{2500};
    Constant topography;
    topography.set(elevation);
    REQUIRE_THAT(topography(41, -111),
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    REQUIRE_THAT(topography.getMinimumAndMaximumElevation().first,
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    REQUIRE_THAT(topography.getMinimumAndMaximumElevation().second,
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    double dElevationDx, dElevationDy;
    REQUIRE_THAT(topography(41, -111, &dElevationDx, &dElevationDy),
                 Catch::Matchers::WithinAbs(elevation, 1.e-10));
    REQUIRE_THAT(dElevationDx,
                 Catch::Matchers::WithinAbs(0, 1.e-10));
    REQUIRE_THAT(dElevationDy,
                 Catch::Matchers::WithinAbs(0, 1.e-10));
    SECTION("clear")
    {
    topography.clear();
    REQUIRE_FALSE(topography.haveTopography());
    }
}

TEST_CASE("ULocator::Topography", "[gridded]")
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
    REQUIRE_NOTHROW(topography.set(xs.size(), xs.data(),
                                   ys.size(), ys.data(),
                                   elevations.size(), elevations.data()));
    REQUIRE_THAT(topography.getMinimumAndMaximumElevation().first,
                 Catch::Matchers::WithinAbs(
                  *std::min_element(elevations.begin(), elevations.end()),
                  1.e-10));
    REQUIRE_THAT(topography.getMinimumAndMaximumElevation().second,
                 Catch::Matchers::WithinAbs(
                  *std::max_element(elevations.begin(), elevations.end()),
                  1.e-10));
    REQUIRE(topography.haveTopography());
    constexpr double tol{1.e-7};
    // Corners
    // 8  9  10 11
    // 4  5   6  7
    // 0  1   2  3
    REQUIRE_THAT(elevations.at(0),
                 Catch::Matchers::WithinAbs(topography.evaluate(5, 2), tol));
    REQUIRE_THAT(elevations.at(11),
                 Catch::Matchers::WithinAbs(topography.evaluate(8, 4), tol));
    REQUIRE_THAT(elevations.at(3),
                 Catch::Matchers::WithinAbs(topography.evaluate(8, 2), tol));
    REQUIRE_THAT(elevations.at(8),
                 Catch::Matchers::WithinAbs(topography.evaluate(5, 4), tol));
    // Out of bounds 
    REQUIRE_THAT(elevations.at(0),
                 Catch::Matchers::WithinAbs(topography.evaluate(5 - 1, 2 - 1), tol));
    REQUIRE_THAT(elevations.at(11),
                 Catch::Matchers::WithinAbs(topography.evaluate(8 + 1, 4 + 1), tol));
    REQUIRE_THAT(elevations.at(3),
                 Catch::Matchers::WithinAbs(topography.evaluate(8 + 1, 2 - 1), tol));
    REQUIRE_THAT(elevations.at(8),
                 Catch::Matchers::WithinAbs(topography.evaluate(5 - 1, 4 + 1), tol));
    // Along bottom, top, left, right
    REQUIRE_THAT(elevations.at(1),
                 Catch::Matchers::WithinAbs(topography.evaluate(6, 2 - 1), tol));
    REQUIRE_THAT(elevations.at(9),
                 Catch::Matchers::WithinAbs(topography.evaluate(6, 4 + 1), tol));
    REQUIRE_THAT(elevations.at(4),
                 Catch::Matchers::WithinAbs(topography.evaluate(5 - 1, 3), tol));
    REQUIRE_THAT(elevations.at(7),
                 Catch::Matchers::WithinAbs(topography.evaluate(8 + 1, 3), tol));

    // Inside model
    REQUIRE_THAT(elevations.at(5),
                 Catch::Matchers::WithinAbs(topography.evaluate(6, 3), tol));
    REQUIRE_THAT(dElevDx*6.8 + dElevDy*3.3,
                 Catch::Matchers::WithinAbs(topography.evaluate(6.8, 3.3), tol));
    REQUIRE_THAT(dElevDx*6.4 + dElevDy*2.1,
                 Catch::Matchers::WithinAbs(topography.evaluate(6.4, 2.1), tol));
    REQUIRE_THAT(dElevDx*6.4 + dElevDy*2,
                 Catch::Matchers::WithinAbs(topography.evaluate(6.4, 2), tol));
    REQUIRE_THAT(dElevDx*7.5 + dElevDy*3.5,
                 Catch::Matchers::WithinAbs(topography.evaluate(7.5, 3.5), tol));
    REQUIRE_THAT(dElevDx*5.1 + dElevDy*3.5,
                 Catch::Matchers::WithinAbs(topography.evaluate(5.1, 3.5), tol));
    double dEdx, dEdy;
    REQUIRE_THAT(dElevDx*6.4 + dElevDy*2.1,
                 Catch::Matchers::WithinAbs(topography.evaluate(6.4, 2.1, &dEdx, &dEdy), tol));
    REQUIRE_THAT(dElevDx,
                 Catch::Matchers::WithinAbs(dEdx, tol));
    REQUIRE_THAT(dElevDy,
                 Catch::Matchers::WithinAbs(dEdy, tol));
}

/*
TEST_CASE("ULocator::Topography", "[loadHDF5]")
{
    ULocator::Position::Utah utah;
    Gridded topography;
    topography.load("utahTopo.h5", utah);
}
*/
