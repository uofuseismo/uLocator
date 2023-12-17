#include <algorithm>
#include <numeric>
#include "optimizers/objectiveFunctions.hpp"
#include "halfSpaceUtilities.hpp"
#include <gtest/gtest.h>

namespace
{

std::vector<double> ATx(const int nRows, const int nColumns,
                        const std::vector<double> &A,
                        const std::vector<double> &x,
                        const double alpha = 1)
{
    std::vector<double> y(nColumns, 0);
    for (int j = 0; j < nColumns; ++j)
    {
        y.at(j) = 0;
        for (int i = 0; i < nRows; ++i)
        {
            auto indx = i*nColumns + j;
            y[j] = y[j] + alpha*(A.at(indx)*x.at(i));
        }
    }
    return y;
}

TEST(ULocatorOptimizersObjectiveFunctions, L1Standard)
{
    std::vector<double> weights{1, 0.5, 3};
    std::vector<double> observations{4, 3, 6};
    std::vector<double> estimates{5, 1, 8};
    EXPECT_NEAR(::l1(weights, observations, estimates,
                     Measurement::Standard), 8, 1.e-10);
}

TEST(ULocatorOptimizersObjectiveFunctions, L1DoubleDifference)
{
    std::vector<double> weights{1, 0.5, 3}; 
    std::vector<double> observations{4, 3, 6}; 
    std::vector<double> estimates{5, 1, 8}; 
    double reference = std::sqrt(1*0.5)*std::abs( (4 - 3) - (5 - 1) )
                     + std::sqrt(1*3)*std::abs( (4 - 6) - (5 - 8) )
                     + std::sqrt(0.5*3)*std::abs( (3 - 6) - (1 - 8) );
    EXPECT_NEAR(::l1(weights, observations, estimates,
                     Measurement::DoubleDifference), reference, 1.e-10);
}

TEST(ULocatorOptimizersObjectiveFunctions, LeastSquaresStandard)
{
    std::vector<double> weights{0.1, 0.2, 0.3};
    std::vector<double> observations{7, 1, 4};
    std::vector<double> estimates{3, -1, 5};
    EXPECT_NEAR(::leastSquares(weights, observations, estimates,
                               Measurement::Standard), 0.41, 1.e-10);
}

TEST(ULocatorOptimizersObjectiveFunctions, LeastSquaresDoubleDifference)
{
    std::vector<double> weights{0.1, 0.2, 0.3};
    std::vector<double> observations{7, 1, 4};
    std::vector<double> estimates{3, -1, 5};
    double reference = std::pow(std::sqrt(0.1*0.2)*( (7 - 1) - ( 3 - -1) ), 2)
                     + std::pow(std::sqrt(0.1*0.3)*( (7 - 4) - ( 3 -  5) ), 2)
                     + std::pow(std::sqrt(0.2*0.3)*( (1 - 4) - (-1 -  5) ), 2);
    EXPECT_NEAR(::leastSquares(weights, observations, estimates,
                               Measurement::DoubleDifference),
                reference, 1.e-10);
}

TEST(ULocatorOptimizersObjectiveFunctions, LpStandard)
{
    constexpr double p{1.5};
    std::vector<double> weights{0.1, 0.3, 0.5};
    std::vector<double> observations{4, 9, 2};
    std::vector<double> estimates{3, 11, 0.5};
    EXPECT_NEAR(::lp(weights, observations, estimates,
                p, Measurement::Standard), 1.0950428241353785, 1.e-10);
}

TEST(ULocatorOptimizersObjectiveFunctions, LpDoubleDifference)
{
    constexpr double p{1.5};
    std::vector<double> weights{0.1, 0.3, 0.5};
    std::vector<double> observations{4, 9, 2}; 
    std::vector<double> estimates{3, 11, 0.5};
    double reference
        = std::pow(
             std::pow(std::sqrt(0.1*0.3)*std::abs((4 - 9) - (3 - 11)),   p)
           + std::pow(std::sqrt(0.1*0.5)*std::abs((4 - 2) - (3 - 0.5)),  p)
           + std::pow(std::sqrt(0.3*0.5)*std::abs((9 - 2) - (11 - 0.5)), p),
           1/p);
    EXPECT_NEAR(::lp(weights, observations, estimates,
                     p, Measurement::DoubleDifference),
                reference, 1.e-10);
}

TEST(ULocatorOptimizersObjectiveFunctions, ThreeDimensionsAndTime)
{
    constexpr int nParameters{4};
    constexpr double p{1.5};
    ::LeastSquaresProblem3DAndTime obj;
    ::fillObjectiveFunction(obj);
    // This is an in-region test
    std::vector<double> x{1, 2, 3, 4};
    std::vector<double> gradient(nParameters, 0);
    std::vector<double> estimates;
    double referenceL1{0};
    double referenceL2{0};
    double referenceLp{0};
    std::vector<double> jacobian(4*obj.mObservations.size());
    std::vector<double> l1Residuals(obj.mObservations.size());
    std::vector<double> l2Residuals(obj.mObservations.size());
    int i = 0;
    for (const auto &stationPhase : obj.mStationPhases)
    {
        constexpr bool applyCorrection{true};
        auto station = stationPhase.first;
        auto phase = stationPhase.second;
        double dtdt0, dtdx, dtdy, dtdz;
        auto estimate
           = obj.mTravelTimeCalculatorMap->evaluate(
                station, phase,
                x.at(0), x.at(1), x.at(2), x.at(3),
                &dtdt0, &dtdx, &dtdy, &dtdz,
                applyCorrection);
        auto xStation = station.getLocalCoordinates().first;
        auto yStation = station.getLocalCoordinates().second;
        auto zStation =-station.getElevation();
        double deltaX = xStation - x.at(1);
        double deltaY = yStation - x.at(2);
        double deltaZ = zStation - x.at(3);
        double distance = std::sqrt( std::pow(deltaX, 2)
                                   + std::pow(deltaY, 2)
                                   + std::pow(deltaZ, 2) );
        double velocity = P_VELOCITY;
        if (phase == "S"){velocity = S_VELOCITY;}
        jacobian.at(4*i + 0) = dtdt0;
        jacobian.at(4*i + 1) =-deltaX/(distance*velocity);
        jacobian.at(4*i + 2) =-deltaY/(distance*velocity);
        jacobian.at(4*i + 3) =-deltaZ/(distance*velocity);
        EXPECT_NEAR(dtdt0, 1, 1.e-14);
        estimates.push_back(estimate); 
        double residual = obj.mObservations.at(i) - estimate;
        double weightedResidual = obj.mWeights.at(i)*residual;
        double l1AbsoluteWeightedResidual = std::abs(weightedResidual);
        double l2WeightedResidual = std::pow(weightedResidual, 2);
        referenceL1 = referenceL1 + l1AbsoluteWeightedResidual;
        referenceL2 = referenceL2 + l2WeightedResidual;
        referenceLp = referenceLp + std::pow(l1AbsoluteWeightedResidual, p);
        l1Residuals.at(i) =-1*obj.mWeights[i]*::sgn(weightedResidual);
        l2Residuals.at(i) =-2*std::pow(obj.mWeights[i], 2)*residual;
        i = i + 1;
    }
    referenceLp = std::pow(referenceLp, 1./p);
    // g = -J^T r 
    auto gradientL1 = ::ATx(l1Residuals.size(), nParameters, jacobian, l1Residuals, -1); 
    auto gradientL2 = ::ATx(l2Residuals.size(), nParameters, jacobian, l2Residuals, -1);
    obj.mNorm = ::Norm::LeastSquares;
    auto logLikelihood = obj.logLikelihood(x.size(), x.data());
    EXPECT_NEAR(logLikelihood, -referenceL2, 1.e-4);
    // Compute gradient
    logLikelihood
        = obj.logLikelihoodAndGradient(x.size(), x.data(), gradient.data());
    EXPECT_NEAR(logLikelihood, -referenceL2, 1.e-4);
    std::vector<double> gradientFiniteDifference(4);
    obj.finiteDifference(x.size(), x.data(), gradientFiniteDifference.data());
    EXPECT_NEAR(logLikelihood, -referenceL2, 1.e-4); 
    for (int i = 0; i < nParameters; ++i)
    {
        EXPECT_NEAR(gradientL2[i], gradient[i], 1.e-8);
        EXPECT_NEAR(gradientL2[i], gradientFiniteDifference[i], 1);
    }
    //------------------------------------l1----------------------------------//
    obj.mNorm = ::Norm::L1;
    logLikelihood = obj.logLikelihood(x.size(), x.data());
    EXPECT_NEAR(logLikelihood, -referenceL1, 1.e-4);
    // Compute the gradient
    logLikelihood
        = obj.logLikelihoodAndGradient(x.size(), x.data(), gradient.data());
    EXPECT_NEAR(logLikelihood, -referenceL1, 1.e-4);
    obj.finiteDifference(x.size(), x.data(), gradientFiniteDifference.data());
    for (int i = 0; i < nParameters; ++i)
    {
        EXPECT_NEAR(gradientL1[i], gradient[i], 1.e-8);
        EXPECT_NEAR(gradientL1[i], gradientFiniteDifference[i], 1.e-3);
    }
    //------------------------------------lp----------------------------------//
    obj.mNorm = ::Norm::Lp;
    obj.mPNorm = p;
    logLikelihood = obj.logLikelihood(x.size(), x.data());
    EXPECT_NEAR(logLikelihood, -referenceLp, 1.e-4);
    // Compute gradient.  This is a PITA.  By now we've confirmed the
    // finite differencing works.  So let's just compare with that.
    logLikelihood
        = obj.logLikelihoodAndGradient(x.size(), x.data(), gradient.data());
    EXPECT_NEAR(logLikelihood, -referenceLp, 1.e-4);
    obj.finiteDifference(x.size(), x.data(), gradientFiniteDifference.data());
    for (int i = 0; i < nParameters; ++i)
    {
        if (i == 0)
        {
            EXPECT_NEAR(gradient[i], gradientFiniteDifference[i], 1.e-3);
        }
        else
        {
            EXPECT_NEAR(gradient[i], gradientFiniteDifference[i], 1.e-6);
        }
    }
    //------------------------------------------------------------------------//
    //                            Test the Penalty                            //
    //------------------------------------------------------------------------//
    double outOfBoundsSource{-4400};
    x[1] = 4000;
    x[2] = 5000;
    double dedx, dedy;
    auto elevation =-obj.mTopography->evaluate(x[1], x[2], &dedx, &dedy);
    dedx =-dedx; // Flip orientation
    dedy =-dedy; // Flip orientation
    x.at(3) = outOfBoundsSource; // Blow this out
    referenceL1 = 0;
    referenceL2 = 0;
    referenceLp = 0;
    i = 0;
    for (const auto &stationPhase : obj.mStationPhases)
    {
        constexpr bool applyCorrection{true};
        auto station = stationPhase.first;
        auto phase = stationPhase.second;
        double dtdt0, dtdx, dtdy, dtdz;
        auto estimate
           = obj.mTravelTimeCalculatorMap->evaluate(
                station, phase,
                x.at(0), x.at(1), x.at(2), x.at(3),
                &dtdt0, &dtdx, &dtdy, &dtdz,
                applyCorrection);
        auto xStation = station.getLocalCoordinates().first;
        auto yStation = station.getLocalCoordinates().second;
        auto zStation =-station.getElevation();
        double deltaX = xStation - x.at(1);
        double deltaY = yStation - x.at(2);
        double deltaZ = zStation - x.at(3);
        double distance = std::sqrt( std::pow(deltaX, 2)
                                   + std::pow(deltaY, 2)
                                   + std::pow(deltaZ, 2) );
        double velocity = P_VELOCITY;
        if (phase == "S"){velocity = S_VELOCITY;}
        jacobian.at(4*i + 0) = dtdt0;
        jacobian.at(4*i + 1) =-deltaX/(distance*velocity);
        jacobian.at(4*i + 2) =-deltaY/(distance*velocity);
        jacobian.at(4*i + 3) =-deltaZ/(distance*velocity);
        EXPECT_NEAR(dtdt0, 1, 1.e-14);
        estimates.push_back(estimate);
        double residual = obj.mObservations.at(i) - estimate;
        double weightedResidual = obj.mWeights.at(i)*residual;
        double l1AbsoluteWeightedResidual = std::abs(weightedResidual);
        double l2WeightedResidual = std::pow(weightedResidual, 2);
        referenceL1 = referenceL1 + l1AbsoluteWeightedResidual;
        referenceL2 = referenceL2 + l2WeightedResidual;
        referenceLp = referenceLp + std::pow(l1AbsoluteWeightedResidual, p);
        l1Residuals.at(i) =-1*obj.mWeights[i]*::sgn(weightedResidual);
        l2Residuals.at(i) =-2*std::pow(obj.mWeights[i], 2)*residual;
        i = i + 1;
    }
    referenceLp = std::pow(referenceLp, 1./p);
    // g = -J^T r 
    gradientL1 = ::ATx(l1Residuals.size(), nParameters, jacobian, l1Residuals, -1);
    gradientL2 = ::ATx(l2Residuals.size(), nParameters, jacobian, l2Residuals, -1);
    obj.mNorm = ::Norm::LeastSquares;
    logLikelihood = obj.logLikelihood(x.size(), x.data());
    auto l1Penalty = obj.mPenaltyCoefficient*(elevation - outOfBoundsSource);
    auto l2Penalty = obj.mPenaltyCoefficient*(elevation - outOfBoundsSource);
    auto lpPenalty = obj.mPenaltyCoefficient*(elevation - outOfBoundsSource);
    gradientL2[1] = gradientL2[1] - obj.mPenaltyCoefficient*dedx;
    gradientL2[2] = gradientL2[2] - obj.mPenaltyCoefficient*dedy;
    gradientL2[3] = gradientL2[3] + obj.mPenaltyCoefficient;
    gradientL1[1] = gradientL1[1] - obj.mPenaltyCoefficient*dedx;
    gradientL1[2] = gradientL1[2] - obj.mPenaltyCoefficient*dedy;
    gradientL1[3] = gradientL1[3] + obj.mPenaltyCoefficient;
    EXPECT_NEAR(logLikelihood, -(referenceL2 + l2Penalty), 1.e-4);
    logLikelihood
        = obj.logLikelihoodAndGradient(x.size(), x.data(), gradient.data());
    EXPECT_NEAR(logLikelihood, -(referenceL2 + l2Penalty), 1.e-4);
    obj.finiteDifference(x.size(), x.data(), gradientFiniteDifference.data());
    for (int i = 0; i < nParameters; ++i)
    {   
        EXPECT_NEAR(gradientL2[i], gradient[i], 1.e-8);
        EXPECT_NEAR(gradientL2[i], gradientFiniteDifference[i], 5.e-1); 
    } 
    //---------------------------------------l1-------------------------------//
    obj.mNorm = ::Norm::L1;
    logLikelihood = obj.logLikelihood(x.size(), x.data());
    EXPECT_NEAR(logLikelihood, -(referenceL1 + l1Penalty), 1.e-4);
    // Compute the gradient
    logLikelihood
        = obj.logLikelihoodAndGradient(x.size(), x.data(), gradient.data());
    obj.finiteDifference(x.size(), x.data(), gradientFiniteDifference.data());
    for (int i = 0; i < nParameters; ++i)
    {
        EXPECT_NEAR(gradientL1[i], gradient[i], 1.e-8);
        EXPECT_NEAR(gradientL1[i], gradientFiniteDifference[i], 5.e-1);
    }
    //---------------------------------------lp-------------------------------//
    obj.mNorm = ::Norm::Lp;
    logLikelihood = obj.logLikelihood(x.size(), x.data());
    EXPECT_NEAR(logLikelihood, -(referenceLp + lpPenalty), 1.e-4);
    logLikelihood
        = obj.logLikelihoodAndGradient(x.size(), x.data(), gradient.data());
    obj.finiteDifference(x.size(), x.data(), gradientFiniteDifference.data());
    for (int i = 0; i < nParameters; ++i)
    {
        if (i == 0)
        {   
            EXPECT_NEAR(gradient[i], gradientFiniteDifference[i], 1.e-3);
        }
        else
        {   
            EXPECT_NEAR(gradient[i], gradientFiniteDifference[i], 1.e-6);
        }
    }

}

}
