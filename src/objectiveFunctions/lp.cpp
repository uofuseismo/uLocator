#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <Eigen/Dense>
#include <boost/math/tools/minima.hpp>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "../logging/standardOut.hpp"
#endif
#include "uLocator/objectiveFunctions/lp.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"
#include "optimizeSourceTime.hpp"
#include "sgn.hpp"

using namespace ULocator::ObjectiveFunctions;

namespace
{
struct SourceTimeObjectiveFunction
{
    double operator()(const double &x)
    {
        double pNorm{0};
        for (int i = 0; i < mResiduals.size(); ++i)
        {
            // 1/sigma | t_obs_i - (t_est_i + t_0) }
            // weights | residual - t_0 |
            double residual = mWeights[i]*std::abs(mResiduals[i] - x);
            //std::cout << std::setprecision(16) << x << " " << residual << std::endl;
            pNorm = pNorm + std::pow(residual, mP); 
        }
//std::cout << "------" << pNorm << "-----" << std::endl;
        return std::pow(pNorm, 1./mP);
    }
    std::vector<double> mResiduals;
    std::vector<double> mWeights;
    double mP{1.5};
};
}

class LP::LPImpl
{
public:
    LPImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr, double p = 1.5) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
        if (p < 1)
        {
            mLogger->warn("p norm must be >= 1; using defualt = " + std::to_string(mP));
        }
        else
        {
            mP = p;
        }
    }
    [[nodiscard]] double evaluateEvidence(const Eigen::VectorXd &estimates)
    {
        if (!mUseDiagonalDataPrecisionMatrix)
        {
            throw std::runtime_error("No idea what to do here");
        }
        double sumOfResiduals{0};
        for (int i = 0; i < mData.rows(); ++i)
        {
            auto weightedAbsResidual
                = mDiagonalDataPrecisionMatrix[i]
                 *std::abs(mData[i] - estimates[i]);
            mResiduals[i] = std::pow(weightedAbsResidual, mP);
            sumOfResiduals = sumOfResiduals + mResiduals[i];
        }
        if (mUseDiagonalDataPrecisionMatrix)
        {
            return std::pow(sumOfResiduals, 1./mP);
        }
#ifndef NDEBUG
        assert(false);
#endif
    }
    [[nodiscard]] double evaluatePrior(
        const Eigen::VectorXd &model)
    {
        return 0;
        //auto residuals = model - mInitialModel;
        //return (residuals.transpose()*mModelPrecisionMatrix*residuals)(0, 0);
    }
    [[nodiscard]] double defaultComputeSourceTime(
        const Eigen::VectorXd &)
    {
        return 0;
    }
    [[nodiscard]] double computeSourceTime(const Eigen::VectorXd &estimates)
    {
        // Create some `ballpark' estimate
#ifndef NDEBUG
        assert(estimates.rows() == mData.rows());
#endif
        double reductionTime{mData[0]};
        for (int i = 0; i < static_cast<int> (estimates.size()); ++i)
        {
            mResidualsSourceTime[i] = mData[i] - reductionTime - estimates[i]; 
        }
        auto sumOfWeights = std::accumulate(mWeightsSourceTime.begin(), mWeightsSourceTime.end(), 0.0);
        double t0{0};
        if (std::abs(mP - 2) < 1.e-10)
        {
            t0 = ::optimizeSourceTimeLeastSquares(mResidualsSourceTime,
                                                  mWeightsSourceTime,
                                                  sumOfWeights);
        }
        else if (std::abs(mP - 1) < 1.e-10)
        {
            t0 = ::optimizeSourceTimeL1(mResidualsSourceTime,
                                        mWeightsSourceTime);
        }
        else
        {
            int nBits = std::numeric_limits<double>::digits;
            std::uintmax_t maxIterations{5000};
            struct SourceTimeObjectiveFunction
                 objectiveFunction{mResidualsSourceTime, mWeightsSourceTime, mP};
            std::pair<double, double> result
                = boost::math::tools::brent_find_minima(
                      objectiveFunction, t0 - 100, t0 + 100, nBits, maxIterations);
            t0 = result.first; 
            //std::cout << std::setprecision(16) << t0 + reductionTime << " " << result.first + reductionTime << " " << result.second << std::endl;
        }
        return t0 + reductionTime;
    }
    void computeJacobian(const std::map<std::string, std::vector<double>> &sensitivities,
                         Eigen::MatrixXd &jacobian)
    {
        const std::vector<double> &dtdx = sensitivities.at("dtdx");
        const std::vector<double> &dtdy = sensitivities.at("dtdy");
        auto nObservations = static_cast<int> (mData.size());
        if (mInversionStrategy == InversionStrategy::LatitudeLongitudeDepth)
        {
            const std::vector<double> &dtdz = sensitivities.at("dtdz");
            if (jacobian.rows() != dtdx.size() ||
                jacobian.cols() != 3)
            {
                jacobian.setZero(dtdx.size(), 3);
            }
            for (int i = 0 ; i < nObservations; ++i)
            {
                jacobian(i, 0) =-dtdx[i];
                jacobian(i, 1) =-dtdy[i];
                jacobian(i, 2) =-dtdz[i];
            }
        }
        else if (mInversionStrategy ==
                 InversionStrategy::LatitudeLongitudeFixedDepth ||
                 mInversionStrategy == 
                 InversionStrategy::LatitudeLongitudeDepthAtFreeSurface)
        {
            if (jacobian.rows() != dtdx.size() ||
                jacobian.cols() != 2)
            {
                jacobian.setZero(dtdx.size(), 2);
            }
            for (int i = 0 ; i < nObservations; ++i)
            {
                jacobian(i, 0) =-dtdx[i];
                jacobian(i, 1) =-dtdy[i];
            }
        }
        else
        {
            const std::vector<double> &dtdz = sensitivities.at("dtdz");
            if (jacobian.rows() != dtdx.size() ||
                jacobian.cols() != 4)
            {
                jacobian.setZero(dtdx.size(), 4);
            }
            for (int i = 0; i < nObservations; ++i)
            {
                jacobian(i, 0) =-dtdx[i];
                jacobian(i, 1) =-dtdy[i];
                jacobian(i, 2) =-dtdz[i];
                jacobian(i, 3) =-1;
            }
        }
    }
    void computeGradient(const Eigen::MatrixXd &jacobian,
                         const Eigen::VectorXd &estimates,
                         Eigen::VectorXd &gradient)
    {
        if (gradient.rows() != jacobian.cols())
        {
            gradient.setZero(jacobian.cols());
        }
        for (int i = 0; i < mData.rows(); ++i)
        {
            auto weightedAbsResidual
                = mDiagonalDataPrecisionMatrix[i]
                 *std::abs(mData[i] - estimates[i]);
            mSgnResiduals[i] = (mP*mDiagonalDataPrecisionMatrix[i])
                             *(::sgn(weightedAbsResidual)
                               *std::pow(weightedAbsResidual, mP - 1));
        }
        if (mUseDiagonalDataPrecisionMatrix)
        {
            gradient =-jacobian.transpose()*mSgnResiduals;
        }
        else
        {
            gradient =-jacobian.transpose()*(mDataPrecisionMatrix*mSgnResiduals);
        }
    }

//private:
    std::function<double (const Eigen::VectorXd &estimates) > mEvaluateEvidence
    {
        std::bind(&LPImpl::evaluateEvidence,
                  this,
                  std::placeholders::_1)
    };
/*
    std::function<double (const Eigen::VectorXd &model,
                          const Eigen::VectorXd &initialModel,
                          const Eigen::MatrixXd &modelCovarianceMatrix) > mEvaluatePrior
    {
        std::bind(&LPImpl::evaluatePrior,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3);
    };
*/
    std::function<double (const Eigen::VectorXd &estimates) >
    mDefaultComputeSourceTime
    {
        std::bind(&LPImpl::defaultComputeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<double (const Eigen::VectorXd &estimates) >
    mComputeSourceTime
    {
        std::bind(&LPImpl::computeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<void (const Eigen::MatrixXd &jacobian,
                        const Eigen::VectorXd &estimates,
                        Eigen::VectorXd &gradient) >
    mComputeGradient
    {
        std::bind(&LPImpl::computeGradient,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3)
    };
    std::function<void (const std::map<std::string, std::vector<double>> &sensitivities,
                        Eigen::MatrixXd &jacobian)>
    mComputeJacobian
    {
        std::bind(&LPImpl::computeJacobian,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2)
    };
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    Eigen::MatrixXd mDataPrecisionMatrix;
    Eigen::MatrixXd mModelPrecisionMatrix;
    Eigen::VectorXd mDiagonalDataPrecisionMatrix;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> mDiagonalModelPrecisionMatrix;
    Eigen::VectorXd mInitialModel;
    Eigen::VectorXd mData;
    Eigen::VectorXd mResiduals;
    Eigen::VectorXd mSgnResiduals;
    std::vector<double> mResidualsSourceTime;
    std::vector<double> mWeightsSourceTime;
    std::vector<std::pair<double, int>> mSourceTimeWork;
    double mP{1.5};
    InversionStrategy mInversionStrategy;
    bool mUseDiagonalDataPrecisionMatrix{true};
    bool mUseDiagonalModelPrecisionMatrix{true};
};

LP::LP(double p) :
    IObjectiveFunction(),
    pImpl(std::make_unique<LPImpl> (nullptr, p))
{
}

LP::LP(std::shared_ptr<UMPS::Logging::ILog> &logger, double p) :
    IObjectiveFunction(logger),
    pImpl(std::make_unique<LPImpl> (logger, p))
{
    IObjectiveFunction::setLikelihoodCallbackFunction(pImpl->mEvaluateEvidence);
    IObjectiveFunction::setComputeSourceTimeCallbackFunction(
       pImpl->mComputeSourceTime);
    IObjectiveFunction::setJacobianCallbackFunction(
       pImpl->mComputeJacobian);
    IObjectiveFunction::setGradientCallbackFunction(
       pImpl->mComputeGradient);
    IObjectiveFunction::toggleComputeSourceTimeInFAndG(true);
}

void LP::setArrivals(const std::vector<Arrival> &arrivals)
{
    pImpl->mLogger->debug("Setting arrivals in base class...");
    IObjectiveFunction::setArrivals(arrivals);
    const auto &inputArrivals = IObjectiveFunction::getArrivalsReference();
    pImpl->mLogger->debug("Setting data...");
    auto nArrivals = static_cast<int> (inputArrivals.size());
    pImpl->mResiduals.setZero(nArrivals);
    pImpl->mSgnResiduals.setZero(nArrivals);
    pImpl->mData.resize(nArrivals);
    for (int i = 0; i < nArrivals; ++i)
    {
        pImpl->mData[i] = inputArrivals[i].getTime();
    }
    pImpl->mLogger->debug("Setting diagonal precision matrix...");
    pImpl->mDiagonalDataPrecisionMatrix.resize(nArrivals);
    pImpl->mWeightsSourceTime.resize(nArrivals);
    for (int i = 0; i < nArrivals; ++i)
    {
        auto standardError = inputArrivals[i].getStandardError();
#ifndef NDEBUG
        assert(standardError > 0);
#endif
        pImpl->mDiagonalDataPrecisionMatrix(i) = 1./standardError;
        pImpl->mWeightsSourceTime[i] = 1./standardError;
    }
    pImpl->mUseDiagonalDataPrecisionMatrix = true;
    // Terms for origin time optimization 
    pImpl->mResidualsSourceTime.resize(nArrivals);
    pImpl->mSourceTimeWork.resize(nArrivals);
}

/*
void LP::initialize()
{
auto job = InversionStrategy::LatitudeLongitudeDepth;
    // Set the forward problem
    if (job == InversionStrategy::LatitudeLongitudeDepth)
    {
    }
    //setLikelihoodFunction(forwardProblem);
    // Set the way to tabulate residuals

}
*/

LP::~LP() = default;

double LP::operator()(const unsigned int n,
                      const double *x,
                      double *gradient)
{
    if (gradient == nullptr)
    {
        return f(n, x);
    }
    else
    {
        return g(n, x, gradient);
    }
}

void LP::setInversionStrategy(
    const InversionStrategy strategy) noexcept
{
    pImpl->mInversionStrategy = strategy;
    if (strategy == InversionStrategy::LatitudeLongitudeDepth ||
        strategy == InversionStrategy::LatitudeLongitudeFixedDepth ||
        strategy == InversionStrategy::LatitudeLongitudeDepthAtFreeSurface)
    {
        pImpl->mLogger->debug(
            "Will independently optimize source time during inversion");
        IObjectiveFunction::setComputeSourceTimeCallbackFunction(
            pImpl->mComputeSourceTime);
    }
    else
    {
        pImpl->mLogger->debug(
            "Will jointly optimize source time during inversion");
        IObjectiveFunction::setComputeSourceTimeCallbackFunction(
            pImpl->mDefaultComputeSourceTime);
    }
    IObjectiveFunction::setInversionStrategy(strategy);
    // Prior model
    pImpl->mInitialModel.setZero(getNumberOfModelParameters());
    pImpl->mModelPrecisionMatrix.setZero( getNumberOfModelParameters(),
                                          getNumberOfModelParameters() );
    pImpl->mDiagonalModelPrecisionMatrix.setZero(getNumberOfModelParameters());
}
