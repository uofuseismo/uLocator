#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <Eigen/Dense>
#include <umps/logging/standardOut.hpp>
#include "uLocator/objectiveFunctions/leastSquares.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"

using namespace ULocator::ObjectiveFunctions;

class LeastSquares::LeastSquaresImpl
{
public:
    LeastSquaresImpl(std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    [[nodiscard]] double evaluateEvidence(const Eigen::VectorXd &estimates)
    {
        mResiduals = mData - estimates;
        if (mUseDiagonalDataPrecisionMatrix)
        {
            return mResiduals.dot(
                     mDiagonalDataPrecisionMatrix.cwiseProduct(mResiduals));
        }
        return (mResiduals.transpose()*mDataPrecisionMatrix*mResiduals)(0, 0);
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
        auto numerator = mDiagonalDataPrecisionMatrix.dot(mData - estimates);
        return numerator/mOriginTimeNormalization;
    }
    void computeJacobian(const std::map<std::string, std::vector<double>> &sensitivities,
                         Eigen::MatrixXd &jacobian)
    {
        //mResiduals = mData - estimates;
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
                jacobian(i, 0) =-2*dtdx[i];
                jacobian(i, 1) =-2*dtdy[i];
                jacobian(i, 2) =-2*dtdz[i];
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
                jacobian(i, 0) =-2*dtdx[i];
                jacobian(i, 1) =-2*dtdy[i];
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
                jacobian(i, 0) =-2*dtdx[i];
                jacobian(i, 1) =-2*dtdy[i];
                jacobian(i, 2) =-2*dtdz[i]; 
                jacobian(i, 3) =-2;
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
        mResiduals = mData - estimates;
        if (mUseDiagonalDataPrecisionMatrix)
        {
            gradient =-jacobian.transpose()
                      *(mDiagonalDataPrecisionMatrix.cwiseProduct(mResiduals));
        }
        else
        {
            gradient =-jacobian.transpose()*(mDataPrecisionMatrix*mResiduals);
        }
    }

//private:
    std::function<double (const Eigen::VectorXd &estimates) > mEvaluateEvidence
    {
        std::bind(&LeastSquaresImpl::evaluateEvidence,
                  this,
                  std::placeholders::_1)
    };
/*
    std::function<double (const Eigen::VectorXd &model,
                          const Eigen::VectorXd &initialModel,
                          const Eigen::MatrixXd &modelCovarianceMatrix) > mEvaluatePrior
    {
        std::bind(&LeastSquaresImpl::evaluatePrior,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3);
    };
*/
    std::function<double (const Eigen::VectorXd &estimates) >
    mDefaultComputeSourceTime
    {
        std::bind(&LeastSquaresImpl::defaultComputeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<double (const Eigen::VectorXd &estimates) >
    mComputeSourceTime
    {
        std::bind(&LeastSquaresImpl::computeSourceTime,
                  this,
                  std::placeholders::_1)
    };
    std::function<void (const Eigen::MatrixXd &jacobian,
                        const Eigen::VectorXd &estimates,
                        Eigen::VectorXd &gradient) >
    mComputeGradient
    {
        std::bind(&LeastSquaresImpl::computeGradient,
                  this,
                  std::placeholders::_1,
                  std::placeholders::_2,
                  std::placeholders::_3)
    };
    std::function<void (const std::map<std::string, std::vector<double>> &sensitivities,
                        Eigen::MatrixXd &jacobian)>
    mComputeJacobian
    {
        std::bind(&LeastSquaresImpl::computeJacobian,
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
    double mOriginTimeNormalization{1};
    InversionStrategy mInversionStrategy;
    bool mUseDiagonalDataPrecisionMatrix{true};
    bool mUseDiagonalModelPrecisionMatrix{true};
};

LeastSquares::LeastSquares() :
    IObjectiveFunction(),
    pImpl(std::make_unique<LeastSquaresImpl> (nullptr))
{
}

LeastSquares::LeastSquares(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    IObjectiveFunction(logger),
    pImpl(std::make_unique<LeastSquaresImpl> (logger))
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

void LeastSquares::setArrivals(const std::vector<Arrival> &arrivals)
{
    pImpl->mLogger->debug("Setting arrivals in base class...");
    IObjectiveFunction::setArrivals(arrivals);
    const auto &inputArrivals = IObjectiveFunction::getArrivalsReference();
    pImpl->mLogger->debug("Setting data...");
    auto nArrivals = static_cast<int> (inputArrivals.size());
    pImpl->mResiduals.setZero(nArrivals);
    pImpl->mData.resize(nArrivals);
    for (int i = 0; i < nArrivals; ++i)
    {
        pImpl->mData[i] = inputArrivals[i].getTime();
    }
    pImpl->mLogger->debug("Setting diagonal precision matrix...");
    pImpl->mDiagonalDataPrecisionMatrix.resize(nArrivals);
    for (int i = 0; i < nArrivals; ++i)
    {
        auto standardError = inputArrivals[i].getStandardError();
#ifndef NDEBUG
        assert(standardError > 0);
#endif
        pImpl->mDiagonalDataPrecisionMatrix(i)
           = 1./(standardError*standardError);
    }
    pImpl->mUseDiagonalDataPrecisionMatrix = true;
    // Terms for origin time optimization 
    pImpl->mOriginTimeNormalization = pImpl->mDiagonalDataPrecisionMatrix.sum();
}

/*
void LeastSquares::initialize()
{
auto job = InversionStrategy::LatitudeLongitudeDepth;
    // Set the forward problem
    if (job == InversionStrategy::LatitudeLongitudeDepth)
    {
        auto forwardProblem
            = std::bind(
                 &LeastSquaresImpl::forwardProblemLatitudeLongitudeDepth,
                 this, 
                 std::placeholders::_1,
                 std::placeholders::_2);
    }
    //setLikelihoodFunction(forwardProblem);
    // Set the way to tabulate residuals

}
*/

LeastSquares::~LeastSquares() = default;

double LeastSquares::operator()(const unsigned int n,
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

void LeastSquares::setInversionStrategy(
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
double LeastSquares::f(const int n, const double x[]) const
{
}

double LeastSquares::g(const int n, const double x[], double g[]) const
{
}
*/

/*
int LeastSquares::getNumberOfModelParameters() const
{
    getInversionStrategy();
}
*/

