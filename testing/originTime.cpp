#include <iostream>
#include <vector>
#include "uLocator/optimizers/originTime.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("ULocator::Optimizers::OriginTime", "[leastSquares]")
{
    ULocator::Optimizers::OriginTime originTime;
    originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::LeastSquares, 2);
    std::vector<double> arrivals{101, 102, 103, 104, 105};
    std::vector<double> weights{0.1, 0.2, 0.3, 0.2, 0.1};
    std::vector<double> travelTimes{1, 2, 3, 4, 5};
    originTime.setArrivalTimes(arrivals, weights);
    originTime.setTravelTimes(travelTimes);
    originTime.optimize();
    REQUIRE_THAT(originTime.getTime(),
                 Catch::Matchers::WithinAbs(100, 1.e-10));
    REQUIRE_THAT(originTime.computeObjectiveFunction(),
                 Catch::Matchers::WithinAbs(0, 1.e-10));
}

TEST_CASE("ULocator::Optimizers::OriginTime", "[l1]")
{
    ULocator::Optimizers::OriginTime originTime;
    originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::L1, 1); 
    std::vector<double> arrivals{101, 102, 104, 106, 107};
    std::vector<double> weights{0.1, 0.2, 0.3, 0.2, 0.1};
    std::vector<double> travelTimes{1, 2, 4, 6, 7}; 
    originTime.setArrivalTimes(arrivals, weights);
    originTime.setTravelTimes(travelTimes);
    originTime.optimize();
    REQUIRE_THAT(originTime.getTime(),
                 Catch::Matchers::WithinAbs(100, 1.e-10));
    REQUIRE_THAT(originTime.computeObjectiveFunction(),
                 Catch::Matchers::WithinAbs(0, 1.e-10));
}

TEST_CASE("ULocator::Optimizers::OriginTime", "[lp]")
{
    ULocator::Optimizers::OriginTime originTime;
    originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::Lp, 1.5);
    std::vector<double> arrivals{101, 103, 107, 111, 113};
    std::vector<double> weights{0.3, 0.2, 0.1, 0.2, 0.3};
    std::vector<double> travelTimes{1, 3, 7, 11, 13};
    originTime.setTolerance(1.e-12);
    originTime.setTimeWindow(100);
    originTime.setNumberOfIterations(1000);
    REQUIRE_THAT(originTime.getTolerance(),
                 Catch::Matchers::WithinAbs(1.e-12, 1.e-14));
    REQUIRE_THAT(originTime.getTimeWindow(),
                 Catch::Matchers::WithinAbs(100, 1.e-10));
    REQUIRE(originTime.getNumberOfIterations() == 1000);
    REQUIRE_NOTHROW(originTime.setArrivalTimes(arrivals, weights));
    REQUIRE_NOTHROW(originTime.setTravelTimes(travelTimes));
    originTime.optimize();
    REQUIRE_THAT(originTime.getTime(),
                 Catch::Matchers::WithinAbs(100, 1.e-10));
    REQUIRE_THAT(originTime.computeObjectiveFunction(),
                 Catch::Matchers::WithinAbs(0, 1.e-10)); 
}

