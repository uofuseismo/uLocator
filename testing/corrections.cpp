#include <iostream>
#include <cmath>
#include <array>
#include <vector>
#include <string>
#include <random>
#include "uLocator/corrections/static.hpp"
#include "uLocator/corrections/sourceSpecific.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "weightedMean.hpp"
#include "weightedMedian.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>


using namespace ULocator::Corrections;

namespace
{

std::vector<int>
    findNeighborsIndices(const double x, const double y, const double z,
                         const int nNeighbors,
                         const std::vector<double> &xSources,
                         const std::vector<double> &ySources,
                         const std::vector<double> &zSources)
{
    std::vector<std::pair<double, int>> distances(xSources.size());
    std::vector<int> result(nNeighbors, -1);
    for (int i = 0; i < static_cast<int> (xSources.size()); ++i)
    {
        auto distance = std::sqrt( std::pow(x - xSources[i], 2)
                                 + std::pow(y - ySources[i], 2)
                                 + std::pow(z - zSources[i], 2) );
        distances[i] = std::pair {distance, i}; 
    }
    std::sort(distances.begin(), distances.end(),
              [](const auto &lhs, const auto &rhs)
              {
                  return lhs.first < rhs.first;
              });
    for (int i = 0; i < nNeighbors; ++i)
    {
        result.at(i) = distances.at(i).second;
    }
    return result;
}

}

TEST_CASE("ULocator::Corrections", "[weightedMean]")
{
    std::vector<double> residuals{14.424, 14.421, 14.417, 14.413, 14.41};
    std::vector<double> weights{3058.0, 8826.0, 56705.0, 30657.0, 12984.0};
    REQUIRE_THAT(::weightedMean(residuals, weights),
                 Catch::Matchers::WithinAbs(14.415602815646439, 1.e-12));
}

TEST_CASE("ULocator::Corrections", "[weightedMedian]")
{
    std::vector<std::pair<double, int>> workSpace;
    REQUIRE_THAT(::weightedMedian(std::vector<double> {-2, -1, 1, 3},
                                  std::vector<double> {1, 1, 1, 1},
                                  workSpace),
                                  Catch::Matchers::WithinAbs(0, 1.e-14));
    REQUIRE_THAT(::weightedMedian(std::vector<double> {-2, 1, -1, 2, 4},
                                  std::vector<double> { 1, 1,  1, 1, 1},
                                  workSpace),
                                  Catch::Matchers::WithinAbs(1, 1.e-14));
    REQUIRE_THAT(::weightedMedian(std::vector<double> {1, 2, 3, 4},
                                  std::vector<double> {0.1, 0.4, 0.4, 0.1},
                                  workSpace),
                                  Catch::Matchers::WithinAbs(2.5, 1.e-14));
    REQUIRE_THAT(::weightedMedian(std::vector<double> {1, 4, 5, 6, 11},
                                  std::vector<double> {0.1, 0.2, 0.3, 0.2, 0.1},
                                  workSpace),
                                  Catch::Matchers::WithinAbs(5, 1.e-14));
}

TEST_CASE("ULocator::Corrections", "[static]")
{
    std::string network{"UU"};
    std::string station{"FORU"};
    std::string phase{"S"};
    const double staticCorrection{1};
    const double estimateTime{30};
    Static correction;
    REQUIRE_NOTHROW(correction.setStationNameAndPhase(network, station, phase));
    REQUIRE(correction.getNetwork() == network);
    REQUIRE(correction.getStation() == station);
    REQUIRE(correction.getPhase() == phase);
    correction.setCorrection(staticCorrection);
    REQUIRE_THAT(correction.getCorrection(),
                 Catch::Matchers::WithinAbs(staticCorrection, 1.e-14));
    REQUIRE_THAT(correction(0),
                 Catch::Matchers::WithinAbs( staticCorrection, 1.e-14));
    REQUIRE_THAT(correction(estimateTime),
                 Catch::Matchers::WithinAbs(staticCorrection + estimateTime, 1.e-14));
    double dtdt0, dtdx, dtdy, dtdz;
    REQUIRE_THAT(correction.evaluate(&dtdt0, &dtdx, &dtdy, &dtdz, estimateTime), 
                 Catch::Matchers::WithinAbs(staticCorrection + estimateTime, 1.e-14));
    REQUIRE_THAT(dtdt0,
                 Catch::Matchers::WithinAbs(0, 1.e-14));
    REQUIRE_THAT(dtdx,
                 Catch::Matchers::WithinAbs(0, 1.e-14));
    REQUIRE_THAT(dtdy,
                 Catch::Matchers::WithinAbs(0, 1.e-14));
    REQUIRE_THAT(dtdz,
                 Catch::Matchers::WithinAbs(0, 1.e-14));
    const std::vector<double> observed{3.4, 4.8, 5.1, 4.9};
    const std::vector<double> estimated{3.2, 5.0, 5.3, 4.3};
    const std::vector<double> weights{1./3., 1./3., 4./15., 1./15};
    correction.train(observed, estimated, Static::Method::Mean);
    REQUIRE_THAT(correction.getCorrection(),
                 Catch::Matchers::WithinAbs(0.1, 1.e-10));
    correction.train(observed, estimated, weights, Static::Method::Mean);
    REQUIRE_THAT(correction.getCorrection(),
                 Catch::Matchers::WithinAbs(-0.01333333333333333, 1.e-10));
}

TEST_CASE("ULocator::Corrections", "[sourceSpecific]")
{
    std::mt19937 rng(86754309);
    std::uniform_real_distribution<double> distribution(-0.5, 0.5);
    std::uniform_int_distribution<int> intDistribution(-500, 500);
    ULocator::Position::UtahRegion utah;
    SourceSpecific correction; 
  
    std::string network{"UU"};
    std::string station{"CTU"};
    std::string phase{"P"};
    constexpr int nPhi{5};
    constexpr int nTheta{7};
    constexpr int nNeighbors{nPhi*nTheta};
    constexpr double radius{5000}; // Radius of 5 km
    constexpr double maximumDistance{radius}; // Works b/c random pert of radius
    std::array<double, 4> center1{40.0, -111.95, 5000.0, 10};
    std::array<double, 4> center2{41.0, -111.9,  6000.0, 11};
    std::array<double, 4> center3{42.0, -111.0,  7000.0, 12};
    std::array<std::array<double, 4>, 3>
        centroids( {center1, center2, center3} );
    std::vector<double> xSources, ySources, zSources,
                        estimates, observations, weights, residuals;
    for (const auto &centroid : centroids)
    {
        auto [x, y]
            = utah.geographicToLocalCoordinates(centroid.at(0), centroid.at(1));
        auto z = centroid.at(2);
        auto bias = centroid.at(3);
        for (int iPhi = 0; iPhi < nPhi; ++iPhi)
        {
            for (int iTheta = 0; iTheta < nTheta; ++iTheta)
            {
                double phi = iPhi*(2*M_PI - 0)/(nPhi); // Don't wrap around
                double theta = iTheta*(M_PI - 0)/(nTheta - 1);
                double randomRadius = radius + intDistribution(rng);
                double xi = x + randomRadius*std::sin(theta)*std::cos(phi);
                double yi = y + randomRadius*std::sin(theta)*std::sin(phi);
                double zi = z + randomRadius*std::cos(theta);
                xSources.push_back(xi);
                ySources.push_back(yi);
                zSources.push_back(zi);
                observations.push_back(bias);
                estimates.push_back(bias + distribution(rng));
                auto residual = observations.back() - estimates.back();
                residuals.push_back(residual);
                weights.push_back(0.5 + distribution(rng));
            }
        }
    }

    REQUIRE_NOTHROW(correction.setStationNameAndPhase(network, station, phase));
    REQUIRE_NOTHROW(correction.setNumberOfNeighbors(nNeighbors));
    REQUIRE_NOTHROW(correction.setEvaluationMethod(SourceSpecific::EvaluationMethod::InverseDistanceWeighted));
    REQUIRE_NOTHROW(correction.setMaximumDistance(maximumDistance));
    
    REQUIRE(correction.getNetwork() == network);
    REQUIRE(correction.getStation() == station);
    REQUIRE(correction.getPhase() ==   phase);
    REQUIRE_THAT(correction.getMaximumDistance(),
                 Catch::Matchers::WithinAbs( maximumDistance, 1.e-8));
   
    REQUIRE_NOTHROW(correction.train(xSources, ySources, zSources,
                                     observations, estimates, weights));
    REQUIRE(correction.haveModel());
    std::array<SourceSpecific::EvaluationMethod, 3> evaluationMethods
    {
        SourceSpecific::EvaluationMethod::WeightedAverage,
        SourceSpecific::EvaluationMethod::MaximumDistanceWeighted,
        SourceSpecific::EvaluationMethod::InverseDistanceWeighted
    };
    // Simplest test - importantly; this tests the KNN lookup
    correction.setEvaluationMethod(
        SourceSpecific::EvaluationMethod::WeightedAverage);
    for (const auto &centroid : centroids)
    {
        auto [x, y]
            = utah.geographicToLocalCoordinates(centroid.at(0), centroid.at(1));
        auto z = centroid.at(2);
        auto indices = ::findNeighborsIndices(x, y, z,
                                              nNeighbors,
                                              xSources, ySources, zSources);
        double numerator{0};
        double denominator{0};
        for (const auto &index : indices)
        {
            auto residual = residuals.at(index);
            auto weight = weights.at(index);
            numerator = numerator + residual*weight; 
            denominator = denominator + weight;
        }
        double referenceCorrection = numerator/denominator;
        constexpr double travelTime{10};
        REQUIRE_THAT(correction.evaluate(x, y, z, travelTime), 
                     Catch::Matchers::WithinAbs(travelTime + referenceCorrection, 1.e-10));
        double dcdt0, dcdx, dcdy, dcdz;
        REQUIRE_THAT(correction.evaluate(x, y, z,
                                        &dcdt0, &dcdx, &dcdy, &dcdz,
                                        travelTime),
                     Catch::Matchers::WithinAbs(travelTime + referenceCorrection, 1.e-10)); 
        REQUIRE_THAT(dcdt0, Catch::Matchers::WithinAbs(0, 2.e-16));
        REQUIRE_THAT(dcdx,  Catch::Matchers::WithinAbs(0, 2.e-16));
        REQUIRE_THAT(dcdy,  Catch::Matchers::WithinAbs(0, 2.e-16));
        REQUIRE_THAT(dcdz,  Catch::Matchers::WithinAbs(0, 2.e-16));
    }
    // Have to be a bit more intelligent for weighted distance
    correction.setEvaluationMethod(
        SourceSpecific::EvaluationMethod::MaximumDistanceWeighted);
    for (const auto &centroid : centroids)
    {   
        auto [x, y]
            = utah.geographicToLocalCoordinates(centroid.at(0), centroid.at(1));
        auto z = centroid.at(2);
        auto indices = ::findNeighborsIndices(x, y, z,
                                              nNeighbors,
                                              xSources, ySources, zSources);
        double numerator{0};
        double denominator{0};
        for (const auto &index : indices)
        {
            auto dx = x - xSources.at(index);
            auto dy = y - ySources.at(index);
            auto dz = z - zSources.at(index);
            auto residual = residuals.at(index);
            auto weight = weights.at(index);
            auto distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            if (distance <= maximumDistance)
            {
                numerator = numerator + residual*weight; 
                denominator = denominator + weight;
            }
        }
        double referenceCorrection = numerator/denominator;
        constexpr double travelTime{11};
        REQUIRE_THAT(correction.evaluate(x, y, z, travelTime), 
                     Catch::Matchers::WithinAbs(travelTime + referenceCorrection, 1.e-10));
        double dcdt0, dcdx, dcdy, dcdz;
        REQUIRE_THAT(correction.evaluate(x, y, z,
                                        &dcdt0, &dcdx, &dcdy, &dcdz,
                                        travelTime),
                     Catch::Matchers::WithinAbs(travelTime + referenceCorrection, 1.e-10));
        REQUIRE_THAT(dcdt0, Catch::Matchers::WithinAbs(0, 2.e-16)); 
        REQUIRE_THAT(dcdx,  Catch::Matchers::WithinAbs(0, 2.e-16)); 
        REQUIRE_THAT(dcdy,  Catch::Matchers::WithinAbs(0, 2.e-16)); 
        REQUIRE_THAT(dcdz,  Catch::Matchers::WithinAbs(0, 2.e-16)); 
        // Let's go far away
        REQUIRE_THAT(correction.evaluate(1.e10, 1.e10, 1.e10, 0),
                     Catch::Matchers::WithinAbs(0, 2.e-16)); 
    }
    // This is the hard one...
    correction.setEvaluationMethod(
        SourceSpecific::EvaluationMethod::InverseDistanceWeighted);
    for (const auto &centroid : centroids)
    {   
        auto [x, y]
            = utah.geographicToLocalCoordinates(centroid.at(0), centroid.at(1));
        auto z = centroid.at(2);
        auto indices = ::findNeighborsIndices(x, y, z,
                                              nNeighbors,
                                              xSources, ySources, zSources);
        double numerator{0};
        double denominator{0};
        for (const auto &index : indices)
        {
            auto dx = x - xSources.at(index);
            auto dy = y - ySources.at(index);
            auto dz = z - zSources.at(index);
            auto residual = residuals.at(index);
            auto weight = weights.at(index);
            auto distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            double inverseDistanceWeight{0};
            if (distance > 0){inverseDistanceWeight = 1./distance;}
            if (distance <= maximumDistance)
            {
                numerator = numerator + residual*(weight*inverseDistanceWeight);
                denominator = denominator + (weight*inverseDistanceWeight);
            }
        }
        double referenceCorrection = numerator/denominator;
        constexpr double travelTime{12};
        REQUIRE_THAT(correction.evaluate(x, y, z, travelTime), 
                     Catch::Matchers::WithinAbs(travelTime + referenceCorrection, 1.e-10));
        double dcdt0, dcdx, dcdy, dcdz;
        REQUIRE_THAT(correction.evaluate(x, y, z,
                                         &dcdt0, &dcdx, &dcdy, &dcdz,
                                         travelTime),
                     Catch::Matchers::WithinAbs(travelTime + referenceCorrection, 1.e-10));

        REQUIRE_THAT(dcdt0,
                     Catch::Matchers::WithinAbs(0, 2.e-16));
        // A gnarly thing can happen at the centroid.  Basically,
        // we switch class membership.  So let's move away.
        double xp = x + 520;
        double yp = y + 530;
        double zp = z + 490;
        double c = correction.evaluate(xp, yp, zp,
                                       &dcdt0, &dcdx, &dcdy, &dcdz,
                                       travelTime);
        constexpr double h{1.e-5};
        double cx = correction.evaluate(xp + h, yp, zp, travelTime);
        double cy = correction.evaluate(xp, yp + h, zp, travelTime);
        double cz = correction.evaluate(xp, yp, zp + h, travelTime);
        double dcdxfd = (cx - c)/h;
        double dcdyfd = (cy - c)/h;
        double dcdzfd = (cz - c)/h;
        REQUIRE_THAT(dcdx,
                     Catch::Matchers::WithinAbs(dcdxfd, 1.e-4));
        REQUIRE_THAT(dcdy,
                     Catch::Matchers::WithinAbs(dcdyfd, 1.e-4));
        REQUIRE_THAT(dcdz, 
                     Catch::Matchers::WithinAbs(dcdzfd, 1.e-4));
    }
}

