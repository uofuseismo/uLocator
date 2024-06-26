#include <cmath>
#include <functional>
#include <algorithm>
#include <vector>
#include "uLocator/optimizers/originTime.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "weightedMean.hpp"
#include "weightedMedian.hpp"
#include "objectiveFunctions.hpp"

using namespace ULocator::Optimizers;

namespace
{

template<typename T>
[[nodiscard]]
std::vector<T> computeResiduals(const std::vector<T> &observations,
                                const std::vector<T> &estimates)
{
#ifndef NDEBUG
    assert(observations.size() == estimates.size());
#endif
    std::vector<T> residuals(observations.size());
    std::transform(observations.begin(), observations.end(),
                   estimates.begin(), residuals.begin(), std::minus<T>());
    return residuals;
}

/// @brief Golden section search from Algorithms for Optimization (Alg 3.3)
///        with an additional break for when the bracket becomes very
///        small.
/// @param[in] f   The 1D function to whose minima will be bracketed i.e.,
///                we find a <= x* <= b such that f(a) <= f(x*) <= f(b).
/// @param[in] leftBracket   The start of the initial interval.
/// @param[in] rightBracket  The end of the initial interval.
/// @param[in] nIterations   The number of bracket refinement iterations.
/// @param[in] tolerance     If the width of the interval is less than this
///                          value then we terminate early.
/// @result result.first is the left bracket and result.second is the right
///         bracket.
[[nodiscard]]
std::pair<double, double> 
    goldenSectionSearch(const std::function<double (double x)> &f,
                        const double leftBracket,
                        const double rightBracket,
                        int nIterations = 1500,
                        const double tolerance = 1.e-5)
{
    constexpr double phi{1.618033988749895}; // (1 + sqrt(5))/2
    double a{leftBracket};
    double b{rightBracket};
    double rho = phi - 1;
    double d = rho*b + (1 - rho)*a;
    double yd = f(d);
    for (int i = 0; i < nIterations; ++i)
    {
        double c = rho*a + (1 - rho)*b;
        double yc = f(c);
        if (yc < yd)
        {
            b = d;
            d = c;
            yd = yc;
        }
        else
        {
            a = b;
            b = c;
        } 
        // Is the bracket tight enough?
        if (std::abs(a - b) < tolerance){break;}
    }
    if (a < b)
    {
        return std::pair {a, b};
    }
    else
    {
        return std::pair {b, a};
    } 
} 

/// @brief Optimizes for the origin time in an least-squares optimization which
///        the weighted mean.
/// @param[in] observations   The times, UTC, of the arrivals in seconds since
///                           the epoch.
/// @param[in] estimates      The times, UTC, of the estimates in seconds since
///                           the epoch.
/// @param[in] weights        The corresponding weights.
/// @result The origin time, UTC, in seconds since the epoch to add to the 
///         estimates.
double optimizeOriginTimeLeastSquares(const std::vector<double> &observations,
                                      const std::vector<double> &estimates,
                                      const std::vector<double> &weights)
{
    constexpr double zero{0}; 
    double sumOfWeights = std::accumulate(weights.begin(), weights.end(), zero);
#ifndef NDEBUG
    assert(sumOfWeights > 0); 
#endif
    auto residuals = ::computeResiduals(observations, estimates);
    return ::weightedMean(residuals, weights);
}

/// @brief Optimizes for the origin time in an Lp norm optimization.
/// @param[in] observations   The times, UTC, of the arrivals in seconds since
///                           the epoch.
/// @param[in] estimates      The times, UTC, of the estimates in seconds since
///                           the epoch.
/// @param[in] weights        The corresponding weights.
/// @param[in] p              The norm.
/// @param[in] timeWindow     The time window to search in seconds before the
///                           the first arrival - i.e., we refined the search
///                           bracket [firstArrival - timeWindow, firstArrival)
///                           using the golden-section search.
/// @param[in] maxIterations  The maximum number of iterations to perform.
/// @param[in] tolerance      If the window width is less than this size in 
///                           seconds then the golden-section search is
///                           terminated.
/// @result The origin time, UTC, in seconds since the epoch to add to the 
///         estimates.
[[nodiscard]] [[maybe_unused]]
double optimizeOriginTimeLp(const std::vector<double> &observations,
                            const std::vector<double> &estimates,
                            const std::vector<double> &weights,
                            const double p = 1.5,
                            const double timeWindow = 200,
                            const int maxIterations = 1500,
                            const double tolerance = 1.e-5)
{
    /// Initialize origin structure
    struct OriginTimeObjectiveFunction
    {
        OriginTimeObjectiveFunction(const std::vector<double> &observationsIn,
                                    const std::vector<double> &estimatesIn,
                                    const std::vector<double> &weightsIn,
                                    const double pIn) :
            observations(observationsIn),
            estimates(estimatesIn),
            weights(weightsIn),
            p(pIn)
        {
            estimatesWork.resize(observations.size());
        }
        /// Evaluate objective function
        double operator()(double const &x) 
        {
            std::transform(estimates.begin(), estimates.end(),
                           estimatesWork.begin(),
                           [&](auto &value)
                           {
                               return value + x;
                           });
            return ::lp(weights, observations, estimatesWork, p,
                        ::Measurement::Standard);
        }
        std::function<double (double)> f
        {
            std::bind(&OriginTimeObjectiveFunction::operator(),
                      this,
                      std::placeholders::_1)
        };
        const std::vector<double> &observations;
        const std::vector<double> &estimates;
        const std::vector<double> &weights;
        std::vector<double> estimatesWork;
        double p{1.5};
    };
    OriginTimeObjectiveFunction
        objectiveFunction{observations, estimates, weights, p}; 
    double tMax = *std::min_element(observations.begin(),
                                    observations.end());
    double tMin = tMax - timeWindow;
    auto bracket
        = ::goldenSectionSearch(objectiveFunction.f, tMin, tMax,
                                maxIterations, tolerance);
    return 0.5*(bracket.first + bracket.second);
}

/// @brief Optimizes for the origin time in an L1 optimization which is the
///        weighted median.
/// @param[in] observations   The times, UTC, of the arrivals in seconds since
///                           the epoch.
/// @param[in] estimates      The times, UTC, of the estimates in seconds since
///                           the epoch.
/// @param[in] weights        The corresponding weights.
/// @result The origin time, UTC, in seconds since the epoch to add to the 
///         estimates.
[[nodiscard]]
double optimizeOriginTimeL1(const std::vector<double> &observations,
                            const std::vector<double> &estimates,
                            const std::vector<double> &weights)
{
    auto residuals = ::computeResiduals(observations, estimates);
    std::vector<std::pair<double, int>> workSpace;
    return ::weightedMedian(residuals, weights, workSpace);
}

}

class OriginTime::OriginTimeImpl
{
public:
    std::vector<double> mObservations;
    std::vector<double> mWeights;
    std::vector<double> mTravelTimes;
    double mReductionTime{0};
    double mOriginTime{0};
    double mTimeWindow{200};
    double mP{1.5};
    double mTolerance{1.e-5}; // Seconds
    ULocator::Optimizers::IOptimizer::Norm
        mNorm{ULocator::Optimizers::IOptimizer::Norm::Lp};
    int mMaxIterations{100};
    bool mHaveOriginTime{false};
    bool mReduceTimes{true};
    bool mHaveTime{false};
};

/// Constructor
OriginTime::OriginTime() :
    pImpl(std::make_unique<OriginTimeImpl> ())
{
}

/// Copy constructor
OriginTime::OriginTime(const OriginTime &originTime)
{
    *this = originTime;
}

/// Move constructor
OriginTime::OriginTime(OriginTime &&originTime) noexcept
{
    *this = std::move(originTime);
}

/// Copy assignment
OriginTime& OriginTime::operator=(const OriginTime &originTime)
{
    if (&originTime == this){return *this;}
    pImpl = std::make_unique<OriginTimeImpl> (*originTime.pImpl);
    return *this;
}

/// Move assignment
OriginTime& OriginTime::operator=(OriginTime &&originTime) noexcept
{
    if (&originTime == this){return *this;}
    pImpl = std::move(originTime.pImpl);
    return *this;
}

/// Iterations
void OriginTime::setNumberOfIterations(const int iterations)
{
    if (iterations < 1)
    {
        throw std::invalid_argument("Number of iterations must be positive");
    }
    pImpl->mMaxIterations = iterations;
}

int OriginTime::getNumberOfIterations() const noexcept
{
    return pImpl->mMaxIterations;
}

/// Tolerance
void OriginTime::setTolerance(const double tolerance)
{
    pImpl->mTolerance = tolerance;
}

double OriginTime::getTolerance() const noexcept
{
    return pImpl->mTolerance;
}

/// Sets the time window
void OriginTime::setTimeWindow(const double window)
{
    if (window <= 0)
    {
        throw std::invalid_argument("The search window must be positive");
    }
    pImpl->mTimeWindow = window;
}

double OriginTime::getTimeWindow() const noexcept
{
    return pImpl->mTimeWindow;
}

/// Reset the class
void OriginTime::clear() noexcept
{
    pImpl = std::make_unique<OriginTimeImpl> ();
}

/// Destructor
OriginTime::~OriginTime() = default;

/// Norm
void OriginTime::setNorm(const IOptimizer::Norm norm, const double p)
{
    if (norm == IOptimizer::Norm::Lp)
    {
        if (p < 1){throw std::invalid_argument("p is less than 1");}
    }
    pImpl->mHaveOriginTime = false;
    if (norm == IOptimizer::Norm::LeastSquares ||
        norm == IOptimizer::Norm::L1 ||
        norm == IOptimizer::Norm::Lp)
    {
        pImpl->mNorm = norm;
    }
    else
    {
        throw std::invalid_argument("Unhandled norm");
    }
    if (norm == IOptimizer::Norm::Lp && std::abs(p - 1) < 1.e-5)
    {
        pImpl->mNorm = IOptimizer::Norm::L1;
    }
    else if (norm == IOptimizer::Norm::Lp && std::abs(p - 2) < 1.e-5)
    {
        pImpl->mNorm = IOptimizer::Norm::LeastSquares;
    }
}

IOptimizer::Norm OriginTime::getNorm() const noexcept
{
    return pImpl->mNorm;
}

void OriginTime::setArrivalTimes(const std::vector<double> &observations,
                                 const std::vector<double> &weights)
{
    if (observations.empty())
    {
        throw std::invalid_argument("At least 1 observation is required");
    }
    if (observations.size() != weights.size())
    {
        throw std::invalid_argument("observations size != weights size"); 
    }
    if (*std::min_element(weights.begin(), weights.end()) <= 0)
    {
        throw std::invalid_argument("All weights must be positive");
    }
    pImpl->mHaveOriginTime = false;
    pImpl->mReductionTime = 0;
    if (reduceTimes())
    {
        ::reduceTimes<double>(observations,
                              &pImpl->mObservations,
                              &pImpl->mReductionTime);
    }
    else
    {
        pImpl->mObservations = observations;
        pImpl->mReductionTime = 0;
    }
    pImpl->mWeights = weights;
}

bool OriginTime::haveArrivalTimes() const noexcept
{
    return !pImpl->mObservations.empty();
}

/// Travel times
void OriginTime::setTravelTimes(const std::vector<double> &travelTimes)
{
    if (travelTimes.empty())
    {
        throw std::invalid_argument("At least 1 travel time is required");
    }
    pImpl->mHaveOriginTime = false;
    pImpl->mTravelTimes = travelTimes;
}

bool OriginTime::haveTravelTimes() const noexcept
{
     return !pImpl->mTravelTimes.empty();
}

/// Origin time
double OriginTime::getTime() const
{
    if (!haveTime()){throw std::runtime_error("Origin time not yet computed");}
    return pImpl->mOriginTime;
}

bool OriginTime::haveTime() const noexcept
{
    return pImpl->mHaveOriginTime;
}

/// Reduce times?
void OriginTime::enableTimeReduction() noexcept
{
    pImpl->mReduceTimes = true;
}

void OriginTime::disableTimeReduction() noexcept
{
    pImpl->mReduceTimes = false;
}

bool OriginTime::reduceTimes() const noexcept
{
    return pImpl->mReduceTimes;
}

/// Optimize
void OriginTime::optimize()
{
    pImpl->mHaveOriginTime = false;
    if (!haveTravelTimes()){throw std::runtime_error("Travel times not set");}
    if (!haveArrivalTimes()){throw std::runtime_error("Arrival times not set");}
    if (pImpl->mObservations.size() != pImpl->mTravelTimes.size())
    {
        throw std::runtime_error("Travel times size != observations size");
    }
    auto norm = getNorm();
    if (norm == IOptimizer::Norm::LeastSquares)
    {
        pImpl->mOriginTime
            = ::optimizeOriginTimeLeastSquares(pImpl->mObservations,
                                               pImpl->mTravelTimes,
                                               pImpl->mWeights)
            + pImpl->mReductionTime;
    }
    else if (getNorm() == IOptimizer::Norm::L1)
    {
        pImpl->mOriginTime = ::optimizeOriginTimeL1(pImpl->mObservations,
                                                    pImpl->mTravelTimes,
                                                    pImpl->mWeights)
                           + pImpl->mReductionTime;
    }
    else if (getNorm() == IOptimizer::Norm::Lp)
    {
        pImpl->mOriginTime = ::optimizeOriginTimeLp(pImpl->mObservations,
                                                    pImpl->mTravelTimes,
                                                    pImpl->mWeights,
                                                    pImpl->mP,
                                                    getTimeWindow(),
                                                    getNumberOfIterations(),
                                                    getTolerance())
                            + pImpl->mReductionTime;
    }
    else
    {
        throw std::runtime_error("Unhandled norm");
    }
    pImpl->mHaveOriginTime = true;
}

