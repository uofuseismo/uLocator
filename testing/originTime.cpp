#include <vector>
#include "uLocator/optimizers/originTime.hpp"
#include <gtest/gtest.h>

TEST(ULocatorOptimizersOriginTime, LeastSquares)
{
    ULocator::Optimizers::OriginTime originTime;
    originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::LeastSquares, 2);
    std::vector<double> arrivals{101, 102, 103, 104, 105};
    std::vector<double> weights{0.1, 0.2, 0.3, 0.2, 0.1};
    std::vector<double> travelTimes{1, 2, 3, 4, 5};
    originTime.setArrivalTimes(arrivals, weights);
    originTime.setTravelTimes(travelTimes);
    originTime.optimize();
    EXPECT_NEAR(originTime.getTime(), 100, 1.e-10);
}

TEST(ULocatorOptimizersOriginTime, L1)
{
    ULocator::Optimizers::OriginTime originTime;
    originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::L1, 1); 
    std::vector<double> arrivals{101, 102, 104, 106, 107};
    std::vector<double> weights{0.1, 0.2, 0.3, 0.2, 0.1};
    std::vector<double> travelTimes{1, 2, 4, 6, 7}; 
    originTime.setArrivalTimes(arrivals, weights);
    originTime.setTravelTimes(travelTimes);
    originTime.optimize();
    EXPECT_NEAR(originTime.getTime(), 100, 1.e-10);
}

TEST(ULocatorOptimizersOriginTime, Lp) 
{
    ULocator::Optimizers::OriginTime originTime;
    originTime.setNorm(ULocator::Optimizers::IOptimizer::Norm::Lp, 1.5);
    std::vector<double> arrivals{101, 103, 107, 111, 113};
    std::vector<double> weights{0.3, 0.2, 0.1, 0.2, 0.3};
    std::vector<double> travelTimes{1, 3, 7, 11, 13};
    originTime.setTolerance(1.e-12);
    originTime.setTimeWindow(100);
    originTime.setNumberOfIterations(1000);
    EXPECT_NEAR(originTime.getTolerance(), 1.e-12, 1.e-14);
    EXPECT_NEAR(originTime.getTimeWindow(), 100, 1.e-10);
    EXPECT_EQ(originTime.getNumberOfIterations(), 1000);
    originTime.setArrivalTimes(arrivals, weights);
    originTime.setTravelTimes(travelTimes);
    originTime.optimize();
    EXPECT_NEAR(originTime.getTime(), 100, 1.e-10);
}

