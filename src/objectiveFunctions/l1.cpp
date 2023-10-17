#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <Eigen/Dense>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "../logging/standardOut.hpp"
#endif
#include "uLocator/objectiveFunctions/l1.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"
#include "optimizeSourceTime.hpp"
#include "sgn.hpp"

using namespace ULocator::ObjectiveFunctions;

class L1::L1Impl
{
public:
    L1Impl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    [[nodiscard]] double evaluateEvidence(const Eigen::VectorXd &estimates)
    {
        for (int i = 0; i < mData.rows(); ++i)
        {
            mResiduals[i] = std::abs(mData[i] - estimates[i]);
        }
        if (mUseDiagonalDataPrecisionMatrix)
        {
            return mDiagonalDataPrecisionMatrix.dot(mResiduals);
        }
        return (mDataPrecisionMatrix*mResiduals).sum();
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
    [[nodiscard]] double computeSourceTime(
        const Eigen::VectorXd &estimates)
    {
#ifndef NDEBUG
        assert(estimates.rows() == mData.rows());
#endif
        for (int i = 0; i < static_cast<int> (estimates.size()); ++i)
        {
            mResidualsSourceTime[i] = mData[i] - estimates[i]; 
        }
        return ::weightedMedian(mResidualsSourceTime, mWeightsSourceTime, mSourceTimeWork);
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
            mSgnResiduals[i] = ::sgn(mData[i] - estimates[i]);
        }
        if (mUseDiagonalDataPrecisionMatrix)
        {
            gradient =-jacobian.transpose()
                      *(mDiagonalDataPrecisionMatrix.cwiseProduct(mSgnResiduals));
        }
        else
        {
            gradient =-jacobian.transpose()*(mDataPrecisionMatrix*mSgnResiduals);
        }
    }

//private:
    std::function<double (const Eigen::VectorXd &estimates) > mEvaluateEvidence
    {
        std::bind(&L1Impl::evaluateEvidence,
                  this,
                  std::placeholders::_1)
    };
/*
    std::function<double (const Eigen::VectorXd &model,
                          const Eigen::VectorXd &initialModel,
                          const Eigen::MatrixXd &modelCovarianceMatrix) > mEvaluatePrior
    {
        std::bind(&L1Impl::evaluatePrior,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3);
    };
*/
    std::function<double (const Eigen::VectorXd &estimates) >
    mDefaultComputeSourceTime
    {
        std::bind(&L1Impl::defaultComputeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<double (const Eigen::VectorXd &estimates) >
    mComputeSourceTime
    {
        std::bind(&L1Impl::computeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<void (const Eigen::MatrixXd &jacobian,
                        const Eigen::VectorXd &estimates,
                        Eigen::VectorXd &gradient) >
    mComputeGradient
    {
        std::bind(&L1Impl::computeGradient,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3)
    };
    std::function<void (const std::map<std::string, std::vector<double>> &sensitivities,
                        Eigen::MatrixXd &jacobian)>
    mComputeJacobian
    {
        std::bind(&L1Impl::computeJacobian,
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
    InversionStrategy mInversionStrategy;
    bool mUseDiagonalDataPrecisionMatrix{true};
    bool mUseDiagonalModelPrecisionMatrix{true};
};

L1::L1() :
    IObjectiveFunction(),
    pImpl(std::make_unique<L1Impl> (nullptr))
{
}

L1::L1(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    IObjectiveFunction(logger),
    pImpl(std::make_unique<L1Impl> (logger))
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

void L1::setArrivals(const std::vector<Arrival> &arrivals)
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
void L1::initialize()
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

L1::~L1() = default;

double L1::operator()(const unsigned int n,
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

void L1::setInversionStrategy(
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
/*
double L1::f(const int n, const double x[]) const
{
}

double L1::g(const int n, const double x[], double g[]) const
{
}
*/

/*
int L1::getNumberOfModelParameters() const
{
    getInversionStrategy();
}
*/

