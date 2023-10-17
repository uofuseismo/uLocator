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
#include "uLocator/objectiveFunctions/doubleDifferenceL1.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"
#include "optimizeSourceTime.hpp"
#include "sgn.hpp"

using namespace ULocator::ObjectiveFunctions;

class DoubleDifferenceL1::DoubleDifferenceL1Impl
{
public:
    DoubleDifferenceL1Impl(
        std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    [[nodiscard]] double evaluateEvidence(const Eigen::VectorXd &estimates)
    {
        for (int i = 0; i < mNumberOfObservations; ++i)
        {
            for (int j = i + 1; j < mNumberOfObservations; ++j)
            {
                auto estimate = estimates[i] - estimates[j];
                auto indx = ( i*(2*mNumberOfObservations - i - 3) )/2 + j - 1;
                mDoubleDifferenceResiduals[indx]
                    = std::abs(mDoubleDifferenceData[i] - estimate);
            }
        }
        if (mUseDiagonalDataPrecisionMatrix)
        {
            return mDiagonalDoubleDifferenceDataPrecisionMatrix.dot(
                       mDoubleDifferenceResiduals);
        }
        return (mDataPrecisionMatrix*mDoubleDifferenceResiduals).sum();
    }
    [[nodiscard]] double evaluatePrior(
        const Eigen::VectorXd &model)
    {
        return 0;
    }
/*
    [[nodiscard]] double defaultComputeSourceTime(
        const Eigen::VectorXd &)
    {
        return 0;
    }
*/
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
        return ::weightedMedian(mResidualsSourceTime,
                                mWeightsSourceTime, mSourceTimeWork);
    }
    void computeJacobian(const std::map<std::string, std::vector<double>> &sensitivities,
                         Eigen::MatrixXd &jacobian)
    {
        const std::vector<double> &dtdx = sensitivities.at("dtdx");
        const std::vector<double> &dtdy = sensitivities.at("dtdy");
        if (mInversionStrategy == InversionStrategy::LatitudeLongitudeDepth)
        {
            const std::vector<double> &dtdz = sensitivities.at("dtdz");
            if (jacobian.rows() != mNumberOfDoubleDifferenceObservations ||
                jacobian.cols() != 3)
            {
                jacobian.setZero(mNumberOfDoubleDifferenceObservations, 3);
            }
            for (int i = 0 ; i < mNumberOfObservations; ++i)
            {
                for (int j = i + 1; j < mNumberOfObservations; ++j)
                {
                    auto indx
                        = ( i*(2*mNumberOfObservations - i - 3) )/2 + j - 1;
                    jacobian(indx, 0) =-(dtdx[i] - dtdx[j]);
                    jacobian(indx, 1) =-(dtdy[i] - dtdy[j]);
                    jacobian(indx, 2) =-(dtdz[i] - dtdz[j]);
                }
            }
        }
        else if (mInversionStrategy ==
                 InversionStrategy::LatitudeLongitudeFixedDepth ||
                 mInversionStrategy == 
                 InversionStrategy::LatitudeLongitudeDepthAtFreeSurface)
        {
            if (jacobian.rows() != mNumberOfDoubleDifferenceObservations ||
                jacobian.cols() != 2)
            {
                jacobian.setZero(mNumberOfDoubleDifferenceObservations, 2);
            }
            for (int i = 0 ; i < mNumberOfObservations; ++i)
            {
                for (int j = i + 1; j < mNumberOfObservations; ++j)
                {
                    auto indx
                        = ( i*(2*mNumberOfObservations - i - 3) )/2 + j - 1;
                    jacobian(indx, 0) =-(dtdx[i] - dtdx[j]);
                    jacobian(indx, 1) =-(dtdy[i] - dtdy[j]);
                }
            }
        }
        else
        {
            throw std::invalid_argument("Cannot invert for source time");
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
        for (int i = 0; i < mNumberOfObservations; ++i)
        {
            for (int j = i + 1; j < mNumberOfObservations; ++j)
            {
                auto indx = ( i*(2*mNumberOfObservations - i - 3) )/2 + j - 1;
                auto estimate = estimates[i] - estimates[j];
                mDoubleDifferenceSgnResiduals[indx]
                    = ::sgn(mDoubleDifferenceData[i] - estimate);
            }
        }
        if (mUseDiagonalDataPrecisionMatrix)
        {
            gradient =-jacobian.transpose()
                      *(mDiagonalDoubleDifferenceDataPrecisionMatrix
                          .cwiseProduct(
                              mDoubleDifferenceSgnResiduals));
        }
        else
        {
            gradient =-jacobian.transpose()
                      *(mDataPrecisionMatrix*mDoubleDifferenceSgnResiduals);
        }
    }

//private:
    std::function<double (const Eigen::VectorXd &estimates) > mEvaluateEvidence
    {
        std::bind(&DoubleDifferenceL1Impl::evaluateEvidence,
                  this,
                  std::placeholders::_1)
    };
/*
    std::function<double (const Eigen::VectorXd &model,
                          const Eigen::VectorXd &initialModel,
                          const Eigen::MatrixXd &modelCovarianceMatrix) > mEvaluatePrior
    {
        std::bind(&DoubleDifferenceL1Impl::evaluatePrior,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3);
    };
    std::function<double (const Eigen::VectorXd &estimates) >
    mDumbyComputeSourceTime
    {
        std::bind(&DoubleDifferenceL1Impl::defaultComputeSourceTime,
                  this,
                  std::placeholders::_1)
    };
*/
    std::function<double (const Eigen::VectorXd &estimates) >
    mComputeSourceTime
    {
        std::bind(&DoubleDifferenceL1Impl::computeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<void (const Eigen::MatrixXd &jacobian,
                        const Eigen::VectorXd &estimates,
                        Eigen::VectorXd &gradient) >
    mComputeGradient
    {
        std::bind(&DoubleDifferenceL1Impl::computeGradient,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3)
    };
    std::function<void (const std::map<std::string, std::vector<double>> &sensitivities,
                        Eigen::MatrixXd &jacobian)>
    mComputeJacobian
    {
        std::bind(&DoubleDifferenceL1Impl::computeJacobian,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2)
    };
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    Eigen::MatrixXd mDataPrecisionMatrix;
    Eigen::MatrixXd mModelPrecisionMatrix;
    Eigen::VectorXd mDiagonalDoubleDifferenceDataPrecisionMatrix;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic> mDiagonalModelPrecisionMatrix;
    Eigen::VectorXd mInitialModel;
    Eigen::VectorXd mDoubleDifferenceData;
    Eigen::VectorXd mData;
    Eigen::VectorXd mDoubleDifferenceResiduals;
    Eigen::VectorXd mDoubleDifferenceSgnResiduals;
    std::vector<double> mResidualsSourceTime;
    std::vector<double> mWeightsSourceTime;
    std::vector<std::pair<double, int>> mSourceTimeWork;
    InversionStrategy mInversionStrategy;
    int mNumberOfObservations{0};
    int mNumberOfDoubleDifferenceObservations{0};
    bool mUseDiagonalDataPrecisionMatrix{true};
    bool mUseDiagonalModelPrecisionMatrix{true};
};

DoubleDifferenceL1::DoubleDifferenceL1() :
    IObjectiveFunction(),
    pImpl(std::make_unique<DoubleDifferenceL1Impl> (nullptr))
{
}

DoubleDifferenceL1::DoubleDifferenceL1(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    IObjectiveFunction(logger),
    pImpl(std::make_unique<DoubleDifferenceL1Impl> (logger))
{
    IObjectiveFunction::setLikelihoodCallbackFunction(pImpl->mEvaluateEvidence);
    IObjectiveFunction::setComputeSourceTimeCallbackFunction(
       pImpl->mComputeSourceTime);
    IObjectiveFunction::setJacobianCallbackFunction(
       pImpl->mComputeJacobian);
    IObjectiveFunction::setGradientCallbackFunction(
       pImpl->mComputeGradient);
    IObjectiveFunction::toggleComputeSourceTimeInFAndG(false);
}

void DoubleDifferenceL1::setArrivals(const std::vector<Arrival> &arrivals)
{
    if (arrivals.empty()){throw std::invalid_argument("No arrivals");}
    pImpl->mLogger->debug("Setting arrivals in base class...");
    IObjectiveFunction::setArrivals(arrivals);
    const auto &inputArrivals = IObjectiveFunction::getArrivalsReference();
    if (inputArrivals.size() < 2)
    {
        clearArrivals();
        throw std::invalid_argument("At least 2 arrivals required");
    }
    pImpl->mLogger->debug("Setting data...");
    auto nArrivals = static_cast<int> (inputArrivals.size());
    pImpl->mNumberOfObservations = nArrivals;
    pImpl->mNumberOfDoubleDifferenceObservations
        = static_cast<int> (((nArrivals - 1)*nArrivals)/2);
    pImpl->mDoubleDifferenceResiduals.setZero(
        pImpl->mNumberOfDoubleDifferenceObservations);
    pImpl->mDoubleDifferenceSgnResiduals.setZero(
        pImpl->mNumberOfDoubleDifferenceObservations);
    pImpl->mDoubleDifferenceData.resize(
        pImpl->mNumberOfDoubleDifferenceObservations);
    pImpl->mData.resize(pImpl->mNumberOfObservations);
    for (int i = 0; i < nArrivals; ++i)
    {
        for (int j = i + 1; j < nArrivals; ++j)
        {
            auto indx = ( i*(2*nArrivals - i - 3) )/2 + j - 1;  
            pImpl->mDoubleDifferenceData(indx)
               = inputArrivals[i].getTime() - inputArrivals[j].getTime();
        }
        pImpl->mData[i] = inputArrivals[i].getTime();
    }
    pImpl->mLogger->debug("Setting diagonal precision matrix...");
    pImpl->mDiagonalDoubleDifferenceDataPrecisionMatrix.resize(
        pImpl->mNumberOfDoubleDifferenceObservations);
    pImpl->mWeightsSourceTime.resize(nArrivals);
    for (int i = 0; i < nArrivals; ++i)
    {
        auto standardError_i = inputArrivals[i].getStandardError();
#ifndef NDEBUG
        assert(standardError_i > 0);
#endif
        for (int j = i + 1; j < nArrivals; ++j)
        {
            auto standardError_j = inputArrivals[j].getStandardError(); 
            auto indx = ( i*(2*nArrivals - i - 3) )/2 + j - 1;
            pImpl->mDiagonalDoubleDifferenceDataPrecisionMatrix(indx)
                = 1./(standardError_i + standardError_j);
                //= 1./std::sqrt(standardError_i*standardError_j);
         }
         pImpl->mWeightsSourceTime[i] = 1./standardError_i;
    }
    pImpl->mUseDiagonalDataPrecisionMatrix = true;
    // Terms for origin time optimization 
    pImpl->mResidualsSourceTime.resize(nArrivals);
    pImpl->mSourceTimeWork.resize(nArrivals);
}

DoubleDifferenceL1::~DoubleDifferenceL1() = default;

double DoubleDifferenceL1::operator()(const unsigned int n,
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

void DoubleDifferenceL1::setInversionStrategy(
    const InversionStrategy strategy)
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
        throw std::invalid_argument("Invalid inversion strategy");
    }
    IObjectiveFunction::setInversionStrategy(strategy);
    // Prior model
    pImpl->mInitialModel.setZero(getNumberOfModelParameters());
    pImpl->mModelPrecisionMatrix.setZero( getNumberOfModelParameters(),
                                          getNumberOfModelParameters() );
    pImpl->mDiagonalModelPrecisionMatrix.setZero(getNumberOfModelParameters());
}
