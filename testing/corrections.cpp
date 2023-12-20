#include <iostream>
#include <vector>
#include <string>
#include "uLocator/corrections/static.hpp"
#include "weightedMean.hpp"
#include <gtest/gtest.h>

using namespace ULocator::Corrections;

namespace
{

TEST(ULocatorCorrections, WeightedMean)
{
    std::vector<double> residuals{14.424, 14.421, 14.417, 14.413, 14.41};
    std::vector<double> weights{3058.0, 8826.0, 56705.0, 30657.0, 12984.0};
    EXPECT_NEAR(::weightedMean(residuals, weights), 14.415602815646439, 1.e-12);
}

TEST(ULocatorCorrections, Static)
{
    std::string network{"UU"};
    std::string station{"FORU"};
    std::string phase{"S"};
    const double staticCorrection{1};
    const double estimateTime{30};
    Static correction;
    EXPECT_NO_THROW(correction.setStationNameAndPhase(network, station, phase));
    EXPECT_EQ(correction.getNetwork(), network);
    EXPECT_EQ(correction.getStation(), station);
    EXPECT_EQ(correction.getPhase(), phase);
    correction.setCorrection(staticCorrection);
    EXPECT_NEAR(correction.getCorrection(), staticCorrection, 1.e-14); 
    EXPECT_NEAR(correction(0), staticCorrection, 1.e-14);
    EXPECT_NEAR(correction(estimateTime),
                staticCorrection + estimateTime, 1.e-14);
    double dtdt0, dtdx, dtdy, dtdz;
    EXPECT_NEAR(correction.evaluate(estimateTime, &dtdt0, &dtdx, &dtdy, &dtdz), 
                staticCorrection + estimateTime, 1.e-14);
    EXPECT_NEAR(dtdt0, 0, 1.e-14);
    EXPECT_NEAR(dtdx,  0, 1.e-14);
    EXPECT_NEAR(dtdy,  0, 1.e-14);
    EXPECT_NEAR(dtdz,  0, 1.e-14);
    const std::vector<double> observed{3.4, 4.8, 5.1, 4.9};
    const std::vector<double> estimated{3.2, 5.0, 5.3, 4.3};
    const std::vector<double> weights{1./3., 1./3., 4./15., 1./15};
    correction.train(observed, estimated, Static::Method::Mean);
    EXPECT_NEAR(correction.getCorrection(), 0.1, 1.e-10);
    correction.train(observed, estimated, weights, Static::Method::Mean);
    EXPECT_NEAR(correction.getCorrection(), -0.01333333333333333, 1.e-10);

}

}
