#include <vector>
#include <cmath>
#define MCMC_ENABLE_EIGEN_WRAPPERS
#include <Eigen/Dense>
//#include <mcmc.hpp>

namespace
{
struct Position
{
    Eigen::VectorXd mSigma;
    std::vector<std::tuple<double, double, double>> mStations;
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    std::vector<double> mStaticCorrections;
    // Tabulates the sum-squared of the residual objective function
    [[nodiscard]]
    double objectiveFunction(const double xs, const double ys,
                             const double zs, const double t0)
    {
        double l2Function{0};
        for (size_t i = 0; i < mStations.size(); ++i)
        {
            auto xr = std::get<0> (mStations[i]);
            auto yr = std::get<1> (mStations[i]); 
            auto zr = std::get<2> (mStations[i]);
            auto tEstimate = t0
                           + estimateTravelTime(xr, yr, zr, xs, ys, zs)
                           + mStaticCorrections[i];
            auto weightedResidual = mWeights[i]*(mObservations[i] - tEstimate);
            l2Function = l2Function + weightedResidual*weightedResidual;               
        }
        return l2Function;
    }
    // Tabulates elements of the gradient
    void computeGradient(const double xs, const double ys,
                         const double zs, const double t0,
                         Eigen::VectorXd *gradientOut)
    {
        double dtdx{0};
        double dtdy{0};
        double dtdz{0};
        double dtdt0{0};
        for (size_t i = 0; i < mStations.size(); ++i)
        {
            auto xr = std::get<0> (mStations[i]);
            auto yr = std::get<1> (mStations[i]); 
            auto zr = std::get<2> (mStations[i]);
            auto dx = xs - xr;
            auto dy = ys - yr;
            auto dz = zs - zr;
            auto tEstimate = t0
                           + estimateTravelTime(xr, yr, zr, xs, ys, zs)
                           + mStaticCorrections[i];
            auto weightedResidual = mWeights[i]*(mObservations[i] - tEstimate);

            auto distance = std::sqrt(dx*dx + dy*dy + dz*dz);
            distance = std::max(std::numeric_limits<double>::epsilon()*100,
                                distance);
            auto dWeightedResidualDx = dx/distance*mSlowness; 
            auto dWeightedResidualDy = dy/distance*mSlowness;
            auto dWeightedResidualDz = dz/distance*mSlowness;
            const double dWeightedResidualDT0{1};
            dtdx  = dtdx - 2*weightedResidual*dWeightedResidualDx;
            dtdy  = dtdy - 2*weightedResidual*dWeightedResidualDy;
            dtdz  = dtdz - 2*weightedResidual*dWeightedResidualDz;
            dtdt0 = dtdt0 - 2*weightedResidual*dWeightedResidualDT0;
        }
        (*gradientOut)(0, 0) = dtdx;
        (*gradientOut)(1, 0) = dtdy;
        (*gradientOut)(2, 0) = dtdz;
        (*gradientOut)(3, 0) = dtdt0;
    }
    double estimateTravelTime(
        const double xr, const double yr, const double zr,
        const double xs, const double ys, const double zs)
    {
        auto dx = xs - xr;
        auto dy = ys - yr;
        auto dz = zs - zr;
        auto distance = std::sqrt(dx*dx + dy*dy + dz*dz);
        auto travelTime = distance*mSlowness;
        return travelTime;
    }
    double mSlowness{1./5000};
};
 
/// Likelihood 
double logLikelihoodDensity(const Eigen::VectorXd &values,
                            Eigen::VectorXd *gradientOut,
                            void *logLikelihoodData)
{
    auto position = reinterpret_cast<::Position *> (logLikelihoodData);
    auto xs = values[0];
    auto ys = values[1];
    auto zs = values[2];
    auto t0 = values[3];
    // log \{ exp(-1/2 \sum [ ((t_i - t_est_i)/sigma_i)^2 ] \}
    auto logLikelihood = -0.5*position->objectiveFunction(xs, ys, zs, t0);
    // Calculate the gradient
    if (gradientOut)
    {
        gradientOut->resize(4, 1);
        position->computeGradient(xs, ys, zs, t0, gradientOut);
    }
    return logLikelihood;
}

/// Prior is uninformative \log ( exp(-0) )
double logPriorDensity(const Eigen::VectorXd &values, void *logLikelihoodData)
{
    return 0;
}

/// log p(m | d) \propto log ( p(d | m) p(m) ) = log p(d | m) + log p(m)
double logPosterior(const Eigen::VectorXd &inputValues,
                    Eigen::VectorXd *outputGradient, void *logLikelihoodData)
{
    return logLikelihoodDensity(inputValues, outputGradient, logLikelihoodData)
         + logPriorDensity(inputValues, logLikelihoodData);
}

}
