#include <iostream>
#include <iomanip>
#include <functional>
#include <memory>
#include <vector>
#include <numeric>
#ifndef NDEBUG
#include <cassert>
#endif
#include <Eigen/Dense>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "../logging/standardOut.hpp"
#endif
#include "uLocator/objectiveFunctions/objectiveFunction.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/travelTimeCalculator.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/topography.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"

using namespace ULocator::ObjectiveFunctions;

class IObjectiveFunction::IObjectiveFunctionImpl
{
public:
    IObjectiveFunctionImpl(std::shared_ptr<UMPS::Logging::ILog> logger) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    void setForwardProblem()
    {
        mHaveForwardProblem = false;
        mModelParameters =-1;
        if (mInversionStrategy == InversionStrategy::LatitudeLongitudeDepth)
        {
            mForwardProblem = std::bind(
                 &IObjectiveFunctionImpl::forwardProblemLatitudeLongitudeDepth,
                 this, 
                 std::placeholders::_1,
                 std::placeholders::_2,
                 std::placeholders::_3,
                 std::placeholders::_4,
                 std::placeholders::_5,
                 std::placeholders::_6);
            mModelParameters = 3;
            mHaveForwardProblem = true;
        } 
        else if (mInversionStrategy ==
                 InversionStrategy::LatitudeLongitudeFixedDepth)
        {
            mForwardProblem = std::bind(
                 &IObjectiveFunctionImpl::
                     forwardProblemLatitudeLongitudeFixedDepth,
                 this, 
                 std::placeholders::_1,
                 std::placeholders::_2,
                 std::placeholders::_3,
                 std::placeholders::_4,
                 std::placeholders::_5,
                 std::placeholders::_6);
            mModelParameters = 2;
            mHaveForwardProblem = true;
        }
        else if (mInversionStrategy ==
                 InversionStrategy::LatitudeLongitudeDepthAtFreeSurface)
        {
            if (!mHaveTopographyCallbackFunction)
            {
                throw std::runtime_error("Topography callback not set!");
            }
            mForwardProblem = std::bind(
                 &IObjectiveFunctionImpl::
                     forwardProblemLatitudeLongitudeDepthAtFreeSurface,
                 this, 
                 std::placeholders::_1,
                 std::placeholders::_2,
                 std::placeholders::_3,
                 std::placeholders::_4,
                 std::placeholders::_5,
                 std::placeholders::_6);
            mModelParameters = 2;
            mHaveForwardProblem = true;
        }
        else
        {
            throw std::runtime_error("Unhandled inversion strategy");
        }
    }
    void forwardProblemLatitudeLongitudeDepth(
        const std::vector<double> &x,
        std::vector<double> *estimates,
        std::vector<double> *dtdx,
        std::vector<double> *dtdy,
        std::vector<double> *dtdz,
        const bool applyCorrection)
    {
#ifndef NDEBUG
        assert(x.size() == 3);
        assert(mTravelTimeCalculatorMap != nullptr);
#endif
        if (estimates->size() != x.size())
        {
            estimates->resize(x.size());
        }
        auto utmX = x[0];
        auto utmY = x[1];
        ULocator::Position::WGS84 epicenter{mUTMZone, mNorth, utmX, utmY};
        auto sourceDepth = x[2];
 
/*
        auto longitude = x[0];
        auto latitude = x[1];
        auto sourceDepth = x[2];
        ULocator::Position::WGS84 epicenter{latitude, longitude, mUTMZone};
*/
        try
        {
            if (dtdx == nullptr && dtdy == nullptr && dtdz == nullptr)
            {
                mTravelTimeCalculatorMap->evaluate(mStationPhasePairs, 
                                                   epicenter,
                                                   sourceDepth, estimates,
                                                   applyCorrection);
            }
            else
            {
                mTravelTimeCalculatorMap->evaluate(mStationPhasePairs,
                                                   epicenter,
                                                   sourceDepth, estimates,
                                                   dtdx, dtdy, dtdz,
                                                   applyCorrection);
            }
        }
        catch (const std::exception &e)
        {
            mLogger->error(e.what());
            throw std::runtime_error("lat/lon/depth forward problem failed: "
                                   + std::string(e.what()));
        }
    }
    void forwardProblemLatitudeLongitudeFixedDepth(
        const std::vector<double> &x, 
        std::vector<double> *estimates,
        std::vector<double> *dtdx,
        std::vector<double> *dtdy,
        std::vector<double> *dtdz,
        const bool applyCorrection)
    {
#ifndef NDEBUG
        assert(x.size() == 2); 
        assert(mTravelTimeCalculatorMap != nullptr);
        assert(mHaveFixedDepthCallbackFunction);
#endif
        if (estimates->size() != x.size())
        {
            estimates->resize(x.size());
        }
        auto utmX = x[0];
        auto utmY = x[1];
        ULocator::Position::WGS84 epicenter{mUTMZone, mNorth, utmX, utmY};
        auto sourceDepth = mFixedDepthCallbackFunction();
/*
        auto longitude = x[0];
        auto latitude = x[1];
        auto sourceDepth = mFixedDepth;
        ULocator::Position::WGS84 epicenter{latitude, longitude, mUTMZone};
*/
        try
        {
            if (dtdx == nullptr && dtdy == nullptr && dtdz == nullptr)
            {
                mTravelTimeCalculatorMap->evaluate(mStationPhasePairs,
                                                   epicenter,
                                                   sourceDepth, estimates,
                                                   applyCorrection);
            }
            else
            {
                mTravelTimeCalculatorMap->evaluate(mStationPhasePairs,
                                                   epicenter,
                                                   sourceDepth, estimates,
                                                   dtdx, dtdy, dtdz,
                                                   applyCorrection);
            }
        }
        catch (const std::exception &e)
        {
            mLogger->error(e.what());
            throw std::runtime_error(
                "lat/lon/fixed depth forward problem failed: "
              + std::string(e.what()));
        }
    }
    void forwardProblemLatitudeLongitudeDepthAtFreeSurface(
        const std::vector<double> &x,
        std::vector<double> *estimates,
        std::vector<double> *dtdx,
        std::vector<double> *dtdy,
        std::vector<double> *dtdz,
        const bool applyCorrection)
    {
#ifndef NDEBUG
        assert(x.size() == 2); 
        assert(mTravelTimeCalculatorMap != nullptr);
#endif
        if (estimates->size() != x.size())
        {
            estimates->resize(x.size());
        }
        auto utmX = x[0];
        auto utmY = x[1];
        ULocator::Position::WGS84 epicenter{mUTMZone, mNorth, utmX, utmY};
        auto latitude = epicenter.getLatitude();
        auto longitude = epicenter.getLongitude();
        auto sourceDepth =-mTopographyCallbackFunction(latitude, longitude);
/*
        auto longitude = x[0];
        auto latitude = x[1];
        auto sourceDepth =-mTopographyCallbackFunction(latitude, longitude);
        ULocator::Position::WGS84 epicenter{latitude, longitude, mUTMZone};
*/
        try
        {
            if (dtdx == nullptr && dtdy == nullptr && dtdz == nullptr)
            {
                mTravelTimeCalculatorMap->evaluate(mStationPhasePairs,
                                                   epicenter,
                                                   sourceDepth, estimates,
                                                   applyCorrection);
            }
            else
            {
                mTravelTimeCalculatorMap->evaluate(mStationPhasePairs,
                                                   epicenter,
                                                   sourceDepth, estimates,
                                                   dtdx, dtdy, dtdz,
                                                   applyCorrection);
            }
        }
        catch (const std::exception &e)
        {
            mLogger->error(e.what());
            throw std::runtime_error(
                "lat/lon/depth free surface forward problem failed: "
              + std::string(e.what()));
        }
    }
/*
    double objectiveFunctionOnly(const int n, const double x[])
    {
        constexpr bool applyCorrection{true};
        if (mModel.size() != n){mModel.resize(n);}
        std::copy(x, x + n, mModel.begin());
        mForwardProblem(mModel, &mEstimates, applyCorrection);
        mTabulateResiduals(mObservations, mEstimates, mWeights, mSumOfWeights,
                           &mResiduals);
        auto objectiveFunction = mApplyNorm(mResiduals);
        return objectiveFunction;
    }
*/
    // Need to know how to compute predictions
    std::function
    <
        void (const std::vector<double> &x,
              std::vector<double> *y,
              std::vector<double> *dtdx,
              std::vector<double> *dtdy,
              std::vector<double> *dtdz,
              const bool applyCorrection)
    > mForwardProblem;
/*
    std::function<void (const std::vector<double> &observations,
                        const std::vector<double> &estimates,
                        const std::vector<double> &weights,
                        const double sumOfWeights,
                        std::vector<double> *residuals)>
        mTabulateResiduals;
*/
    std::function<double () > mFixedDepthCallbackFunction;
    std::function<double (const Eigen::VectorXd &estimates) > mLikelihoodFunction;
    std::function
    <
        double (const Eigen::VectorXd &estimates)
    > mComputeSourceTimeCallbackFunction;
    std::function<
        double (const double latitude, const double longitude)
    > mTopographyCallbackFunction;
    std::function
    <
        void (const Eigen::MatrixXd &jacobian,
              const Eigen::VectorXd &estimates,
              Eigen::VectorXd &gradient)
    > mGradientCallbackFunction;
    std::function
    <
         void (const std::map<std::string, std::vector<double>> &sensitivities,
               Eigen::MatrixXd &jacobian)
    > mJacobianCallbackFunction;
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::unique_ptr<const ULocator::Topography> mTopography{nullptr};
    std::unique_ptr<const ULocator::TravelTimeCalculatorMap>
        mTravelTimeCalculatorMap{nullptr};
    std::vector<std::pair<ULocator::Station, std::string>> mStationPhasePairs;
    std::vector<Arrival> mInputArrivals;
    std::vector<double> mResiduals;
    Eigen::MatrixXd mJacobian;
    Eigen::VectorXd mEstimates;
    Eigen::VectorXd mGradient;
    std::vector<double> mWeights;
    mutable std::vector<double> mModel;
    double mFiniteDifferenceStep{0.5}; // meters
    int mModelParameters{-1};
    int mObjectiveFunctionCounter{0};
    int mGradientCounter{0};
    int mUTMZone{12};
    InversionStrategy mInversionStrategy;
    bool mUseFiniteDifference{false};
    bool mNorth{true};
    bool mApplyCorrection{true};
    bool mOptimizeSourceTime{true};
    bool mHaveInversionStrategy{false};
    bool mHaveTopographyCallbackFunction{false};
    bool mHaveLikelihoodCallbackFunction{false};
    bool mHaveComputeSourceTimeCallbackFunction{false};
    bool mHaveGradientCallbackFunction{false};
    bool mHaveJacobianCallbackFunction{false};
    bool mHaveFixedDepthCallbackFunction{false};
    bool mHaveForwardProblem{false};
    bool mComputeSourceTimeInFAndG{true};
};

/// Constructor
IObjectiveFunction::IObjectiveFunction() :
    pImpl(std::make_unique<IObjectiveFunctionImpl> (nullptr))
{
}

/// Constructor
IObjectiveFunction::IObjectiveFunction(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<IObjectiveFunctionImpl> (logger))
{
}

/// Destructor
IObjectiveFunction::~IObjectiveFunction() = default;

/// Sets the observed arrivals


/// Sets the observations
/*
void IObjectiveFunction::setObservations(
    const std::vector<double> &observations)
{
    std::vector<double> weights(observations.size(), 1.0); 
    setObservations(observations, weights);
}

*/

/*
*/

void IObjectiveFunction::resetCounters() noexcept
{
    pImpl->mObjectiveFunctionCounter = 0;
    pImpl->mGradientCounter = 0;
}
    

bool IObjectiveFunction::haveInversionStrategy() const noexcept
{
    return pImpl->mHaveInversionStrategy;
}

void IObjectiveFunction::setInversionStrategy(
    const InversionStrategy strategy)
{
    pImpl->mInversionStrategy = strategy;
    pImpl->setForwardProblem();
    pImpl->mHaveInversionStrategy = true;
}

void IObjectiveFunction::clearArrivals() noexcept
{
    pImpl->mInputArrivals.clear();
}

void IObjectiveFunction::setArrivals(const std::vector<Arrival> &arrivals)
{
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::invalid_argument("Travel time calculator map not set");
    }
    clearArrivals();
    std::vector<bool> keepArrival(arrivals.size());
    std::vector<std::string> arrivalNames;
    int nArrivals = 0;
    for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
    {
        keepArrival[i] = false;
        if (!arrivals[i].haveStation())
        {
            pImpl->mLogger->error("Arrival " + std::to_string(i)
                               + "'s station not set.  Skipping...");
            continue;
        }
        const auto &station = arrivals[i].getStationReference();
        auto name = station.getNetwork() + "." + station.getName(); 
        if (!arrivals[i].havePhase())
        {
            pImpl->mLogger->warn("Arrival " + name
                               + " does not have phase.  Skipping...");
            continue;
        }
        auto phase = arrivals[i].getPhase();
        name = name + "." + phase;
        if (!arrivals[i].haveTime())
        {
            pImpl->mLogger->warn("Arrival " + name
                               + " does not have time.  Skipping...");
            continue;
        }
        if (!arrivals[i].haveStandardError())
        {
            pImpl->mLogger->warn("Arrival " + name
                               + " does not have standard error.  Skipping...");
            continue;
        }
        if (!pImpl->mTravelTimeCalculatorMap->contains(station, phase))
        {
            pImpl->mLogger->warn("Arrival " + name
                               + " is not in the travel time calculator map."
                               + "  Skipping...");
            continue;
        }
        bool arrivalExists{false};
        for (const auto &arrivalName : arrivalNames)
        {
            if (name == arrivalName)
            {
                pImpl->mLogger->warn("Arrival " + name
                                   + " exists.  Skipping...");
                arrivalExists = true;
                break;
            }
        }
        if (arrivalExists){continue;}
        // Mark it as a keeper
        arrivalNames.push_back(name);
        nArrivals = nArrivals + 1;
        keepArrival[i] = true;
    }
    // Create the station phase pairs and copy the arrivals
    pImpl->mStationPhasePairs.clear();
    pImpl->mInputArrivals.clear();
    pImpl->mStationPhasePairs.reserve(nArrivals);
    for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
    {
        if (keepArrival[i])
        {
            auto station = arrivals[i].getStationReference();
            auto phase = arrivals[i].getPhase();
            if (pImpl->mLogger->getLevel() == UMPS::Logging::Level::Debug)
            {
                auto name = station.getNetwork()
                          + "." + station.getName() + "." + phase;
                auto time = arrivals[i].getTime();
                auto standardError = arrivals[i].getStandardError();
                pImpl->mLogger->debug("Adding: " + name + "; onset time = "
                                    + std::to_string(time) + " +/- "
                                    + std::to_string(standardError) + " (s)");
            }
            pImpl->mStationPhasePairs.push_back(std::pair {station, phase});
            pImpl->mInputArrivals.push_back(arrivals[i]);
        }
    }
    // Reserve auxiliary space
    pImpl->mEstimates.resize(pImpl->mInputArrivals.size());
    pImpl->mLogger->debug("Set " + std::to_string(nArrivals) + " arrivals");
}
 
const std::vector<ULocator::Arrival> &IObjectiveFunction::getArrivalsReference() const
{
    return *&pImpl->mInputArrivals;
}

std::vector<ULocator::Arrival> IObjectiveFunction::getArrivals() const
{
    return pImpl->mInputArrivals;
}

int IObjectiveFunction::getNumberOfArrivals() const
{
    return static_cast<int> (pImpl->mInputArrivals.size());
}

/// Travel time calcuators
void IObjectiveFunction::setTravelTimeCalculatorMap(
    std::unique_ptr<const ULocator::TravelTimeCalculatorMap> &&calculators)
{
    if (calculators == nullptr)
    {    
        throw std::invalid_argument("Travel time calculator is NULL");
    }    
    if (calculators->size() == 0)
    {    
        throw std::invalid_argument("No travel time calculalors");
    }
    clearArrivals();
    pImpl->mTravelTimeCalculatorMap = std::move(calculators);
    calculators = nullptr;
}

/// Release the travel time calculator
std::unique_ptr<const ULocator::TravelTimeCalculatorMap>
IObjectiveFunction::releaseTravelTimeCalculatorMap()
{
    auto result = std::move(pImpl->mTravelTimeCalculatorMap);
    pImpl->mTravelTimeCalculatorMap = nullptr;
    return result;
}

/// Have the calculator?
bool IObjectiveFunction::haveTravelTimeCalculatorMap() const noexcept
{
    return (pImpl->mTravelTimeCalculatorMap != nullptr);
}


/*
void IObjectiveFunction::setModel(const std::vector<double> &x) const
{
    if (x.empty())
    {
        throw std::runtime_error("x is empty");
    }
    setModel(x.size(), x.data());
}
*/

/*
void IObjectiveFunction::evaluateForwardProblem(const int n, const double x[], double y[])
{
    //pImpl->mForwardProblem(n, x, y);
}

void IObjectiveFunction::computeResiduals(const int n, const double x[], const double y[], double residuals[])
{
    pImpl->mComputeResiduals(n, x, y, residuals);
}
*/



double IObjectiveFunction::f(const int n, const double x[]) const
{
    pImpl->mObjectiveFunctionCounter = pImpl->mObjectiveFunctionCounter + 1;
#ifndef NDEBUG
    assert(haveLikelihoodCallbackFunction());
    assert(x != nullptr);
#endif
    if (pImpl->mModel.size() != n)
    {
        pImpl->mModel.resize(n);
    }
    std::copy(x, x + n, pImpl->mModel.data());
    std::vector<double> estimates;
    pImpl->mForwardProblem(pImpl->mModel, &estimates,
                           nullptr, nullptr, nullptr,
                           applyCorrection());
    auto nObservations = static_cast<int> (pImpl->mEstimates.size());
    if (pImpl->mEstimates.size() != estimates.size())
    {
        pImpl->mEstimates.resize(estimates.size());
    }
    std::copy(estimates.begin(), estimates.end(), pImpl->mEstimates.data());
    if (pImpl->mComputeSourceTimeInFAndG)
    {
        double originTime
            = pImpl->mComputeSourceTimeCallbackFunction(pImpl->mEstimates);
        for (int i = 0; i < nObservations; ++i)
        {
            pImpl->mEstimates[i] += originTime;
        }
    }
    // Compute the evidence
    double evidence = pImpl->mLikelihoodFunction(pImpl->mEstimates);
    double prior = 0;
    return evidence + prior;
}

double IObjectiveFunction::g(const int n, const double x[], double g[]) const
{
    pImpl->mGradientCounter = pImpl->mGradientCounter + 1;
#ifndef NDEBUG
    assert(haveLikelihoodCallbackFunction());
    assert(haveGradientCallbackFunction());
    assert(haveJacobianCallbackFunction());
    assert(g != nullptr);
#endif
    if (pImpl->mModel.size() != n)
    {
        pImpl->mModel.resize(n);
    }
    std::copy(x, x + n, pImpl->mModel.data());
    std::vector<double> estimates;
    std::vector<double> dtdx;
    std::vector<double> dtdy;
    std::vector<double> dtdz;
    pImpl->mForwardProblem(pImpl->mModel, &estimates,
                           &dtdx,
                           &dtdy,
                           &dtdz,
                           applyCorrection());
    auto nObservations = static_cast<int> (pImpl->mEstimates.size());
    if (pImpl->mEstimates.size() != estimates.size())
    {
        pImpl->mEstimates.resize(estimates.size());
    }
    std::copy(estimates.begin(), estimates.end(), pImpl->mEstimates.data());
    if (pImpl->mComputeSourceTimeInFAndG)
    {   
        double originTime
            = pImpl->mComputeSourceTimeCallbackFunction(pImpl->mEstimates);
        for (int i = 0; i < nObservations; ++i)
        {   
            pImpl->mEstimates[i] += originTime;
        }   
    }
    // Compute the Jacobian and gradient
    std::map<std::string, std::vector<double>> sensitivities
    {
        std::pair {"dtdx", std::move(dtdx)},
        std::pair {"dtdy", std::move(dtdy)},
        std::pair {"dtdz", std::move(dtdz)}
    };
    try
    {
        pImpl->mJacobianCallbackFunction(sensitivities, pImpl->mJacobian);
    }
    catch (const std::exception &e)
    {
        auto errorMessage = "Failed to compute Jacobian.  Failed with: " 
                          + std::string{e.what()};
        pImpl->mLogger->error(errorMessage);
        throw std::runtime_error(errorMessage);
    }
    try
    {
        pImpl->mGradientCallbackFunction(pImpl->mJacobian,
                                         pImpl->mEstimates,
                                         pImpl->mGradient);
    }
    catch (const std::exception &e)
    {
        auto errorMessage = "Failed to compute gradient.  Failed with: " 
                          + std::string{e.what()};
        pImpl->mLogger->error(errorMessage);
        throw std::runtime_error(errorMessage);
    }
    for (int i = 0; i < n; ++i)
    {
        g[i] = pImpl->mGradient[i];
    }
    // Compute the evidence
    double evidence = pImpl->mLikelihoodFunction(pImpl->mEstimates);
    double prior{0};
    double f0 = evidence + prior;
    if (pImpl->mUseFiniteDifference)
    {
        auto h = pImpl->mFiniteDifferenceStep;
        std::vector<double> xWork(n);
        for (int i = 0; i < n; ++i)
        {
            std::copy(x, x + n, xWork.data());
            xWork[i] = xWork[i] + h;
std::cout << g[i] << " ";
            g[i] = (f(xWork.size(), xWork.data()) - f0)/h;
std::cout << g[i] << std::endl;
        }
getchar();
    }
    return f0;
}

int IObjectiveFunction::getNumberOfModelParameters() const
{
    if (!haveInversionStrategy())
    {
        throw std::invalid_argument("Inversion strategy not set");
    }
    return pImpl->mModelParameters;
}

/// UTM zone
void IObjectiveFunction::setUTMZone(const int zone, const bool north)
{
    if (zone < 0 || zone > 60)
    {
        if (zone !=-1)
        {
            throw std::invalid_argument(
                "UTM zone must be -1 or in range [1,60]");
        }
    }
    pImpl->mUTMZone = zone; 
    pImpl->mNorth = north;
}

int IObjectiveFunction::getUTMZone() const noexcept
{
    return pImpl->mUTMZone;
}

/// Toggles applying the static correction
void IObjectiveFunction::toggleApplyCorrection(
    const bool applyCorrection) noexcept
{
    pImpl->mApplyCorrection = applyCorrection;
}

bool IObjectiveFunction::applyCorrection() const noexcept
{
    return pImpl->mApplyCorrection;
}

/// Topography callback function
void IObjectiveFunction::setTopographyCallbackFunction(
    const std::function<double(const double latitude, const double longitude)> &callback)
{
    pImpl->mTopographyCallbackFunction = callback;
    pImpl->mHaveTopographyCallbackFunction = true;
}

bool IObjectiveFunction::haveTopographyCallbackFunction() const noexcept
{
    return pImpl->mHaveTopographyCallbackFunction;
}

/// Evidence callback
void IObjectiveFunction::setLikelihoodCallbackFunction(
    const std::function<double(const Eigen::VectorXd &estimates)> &callback)
{
    pImpl->mLikelihoodFunction = callback;
    pImpl->mHaveLikelihoodCallbackFunction = true;
}

bool IObjectiveFunction::haveLikelihoodCallbackFunction() const noexcept
{
    return pImpl->mHaveLikelihoodCallbackFunction;
}

/// Source time callback
void IObjectiveFunction::setComputeSourceTimeCallbackFunction(
    const std::function<double(const Eigen::VectorXd &estimates)> &callback)
{
    pImpl->mComputeSourceTimeCallbackFunction = callback;
    pImpl->mHaveComputeSourceTimeCallbackFunction = true;
}

bool IObjectiveFunction::haveComputeSourceTimeCallbackFunction() const noexcept
{
    return pImpl->mHaveComputeSourceTimeCallbackFunction;
}

/// Gradient callback
void IObjectiveFunction::setGradientCallbackFunction(
    const std::function<void(const Eigen::MatrixXd &jacobian,
                             const Eigen::VectorXd &estimates,
                             Eigen::VectorXd &gradient)> &callback)
{
    pImpl->mGradientCallbackFunction = callback;
    pImpl->mHaveGradientCallbackFunction = true;
}

bool IObjectiveFunction::haveGradientCallbackFunction() const noexcept
{
    return pImpl->mHaveGradientCallbackFunction;
}

/// Jacobian callback
void IObjectiveFunction::setJacobianCallbackFunction(
    const std::function<void(
        const std::map<std::string, std::vector<double>> &sensitivities,
        Eigen::MatrixXd &jacobian)> &callback)
{
    pImpl->mJacobianCallbackFunction = callback;
    pImpl->mHaveJacobianCallbackFunction = true;
}

bool IObjectiveFunction::haveJacobianCallbackFunction() const noexcept
{
    return pImpl->mHaveJacobianCallbackFunction;
}

/// Fixed depth callback
void IObjectiveFunction::setFixedDepthCallbackFunction(
   const std::function<double ()> &callback)
{
    pImpl->mFixedDepthCallbackFunction = callback;
    pImpl->mHaveFixedDepthCallbackFunction = true;
}

bool IObjectiveFunction::haveFixedDepthCallbackFunction() const noexcept
{
    return pImpl->mHaveFixedDepthCallbackFunction;
}

double IObjectiveFunction::evaluateLoss(const ULocator::Origin &origin) const
{
    if (!origin.haveEpicenter())
    {
        throw std::invalid_argument("Epicenter not set");
    }
    if (!origin.haveDepth())
    {
        throw std::invalid_argument("Depth not set");
    }
    if (getNumberOfArrivals() < 1)
    {
        throw std::runtime_error("No arrivals set");
    }
    pImpl->mLogger->debug("Predicting travel times for given origin");
    // Compute travel times 
    std::vector<double> travelTimes;
    pImpl->mTravelTimeCalculatorMap->evaluate(
        pImpl->mStationPhasePairs,
        origin.getEpicenter(),
        origin.getDepth(), &travelTimes,
        applyCorrection());
    // Deal with origin time
    double originTime{0};
    if (origin.haveTime())
    {
        originTime = origin.getTime();
    }
    else
    {
        pImpl->mLogger->debug("Origin time not specified - estimating...");
        Eigen::VectorXd estimatesVector(travelTimes.size());
        std::copy(travelTimes.begin(), travelTimes.end(),
                  estimatesVector.data());
        originTime
            = pImpl->mComputeSourceTimeCallbackFunction(estimatesVector);
    }
    // Compute estimates
    Eigen::VectorXd estimates;
    estimates.setZero(travelTimes.size());
    std::copy(travelTimes.begin(), travelTimes.end(), estimates.data());
    for (int i = 0; i < static_cast<int> (travelTimes.size()); ++i)
    {
        estimates[i] += originTime;
    }
    // Compute the evidence
    double evidence = pImpl->mLikelihoodFunction(estimates);
    double prior = 0;
    return evidence + prior;

}

ULocator::Origin
IObjectiveFunction::predict(const ULocator::Origin &origin) const
{
    if (!haveTravelTimeCalculatorMap())
    {
        throw std::runtime_error("Travel time calculator map not set");
    }
    if (!origin.haveEpicenter())
    {
        throw std::invalid_argument("Epicenter not set");
    }
    if (!origin.haveDepth())
    {
        throw std::invalid_argument("Depth not set");
    }
    auto arrivals = getArrivals();
    if (arrivals.empty())
    {
        throw std::runtime_error("No arrivals set");
    }
    auto nArrivals = static_cast<int> (arrivals.size());
    pImpl->mLogger->debug("Predicting travel times for given origin");
    // Copy some stuff to keep origin safe
    auto result = origin;
    auto epicenter = result.getEpicenter();
    auto depth = result.getDepth();
    // Compute travel times 
    std::vector<double> travelTimes;
    pImpl->mTravelTimeCalculatorMap->evaluate(
        pImpl->mStationPhasePairs, epicenter,
        depth, &travelTimes,
        applyCorrection());
    // Deal with origin time
    double originTime{0};
    if (result.haveTime())
    {
        originTime = result.getTime();
    }
    else
    {
        pImpl->mLogger->debug("Origin time not specified - estimating...");
        Eigen::VectorXd estimates(travelTimes.size());
        std::copy(travelTimes.begin(), travelTimes.end(), estimates.data());
        originTime
            = pImpl->mComputeSourceTimeCallbackFunction(estimates);
        result.setTime(originTime);
    }
    // Update arrivals
    for (int i = 0; i < nArrivals; ++i) 
    {
        auto residual = arrivals.at(i).getTime() - (travelTimes[i] + originTime);
        arrivals[i].setResidual(residual);
        auto stationPosition
            = arrivals[i].getStation().getGeographicPosition(); 
        double greatCircleDistance, distance, azimuth, backAzimuth;
        Position::computeDistanceAzimuth(epicenter, stationPosition,
                                         &greatCircleDistance,
                                         &distance,
                                         &azimuth,
                                         &backAzimuth);
        arrivals[i].setDistance(distance);
        arrivals[i].setAzimuth(azimuth);
        arrivals[i].setBackAzimuth(backAzimuth);
    }
    // Set them
    result.setArrivals(std::move(arrivals));
    return result;
}

int IObjectiveFunction::getNumberOfObjectiveFunctionEvaluations() const noexcept
{
    return pImpl->mObjectiveFunctionCounter;
}

int IObjectiveFunction::getNumberOfGradientEvaluations() const noexcept
{
    return pImpl->mGradientCounter;
}

void IObjectiveFunction::toggleComputeSourceTimeInFAndG(
    const bool computeSourceTime) noexcept
{
    pImpl->mComputeSourceTimeInFAndG = computeSourceTime;
}

void IObjectiveFunction::toggleUseFiniteDifference(
    const bool useFiniteDifference, const double perturbation)
{
    if (perturbation == 0)
    {
        throw std::invalid_argument("Perturbation must be non-negative");
    }
    pImpl->mFiniteDifferenceStep = perturbation;
    pImpl->mUseFiniteDifference = useFiniteDifference;
}
