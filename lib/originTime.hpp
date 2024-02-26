#ifndef ORIGIN_TIME_HPP
#define ORIGIN_TIME_HPP
#include <iostream>
#include <iomanip>
#include <functional>
#include <limits>
#include <cmath>
#include <vector>
#include <algorithm>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "weightedMean.hpp"
#include "weightedMedian.hpp"
#include "optimizers/objectiveFunctions.hpp"
namespace
{

template<typename T>
[[nodiscard]] [[maybe_unused]]
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
[[nodiscard]] [[maybe_unused]]
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

/// @brief Optimizes for the origin time in an L1 optimization which is the
///        weighted median.
/// @param[in] observations   The times, UTC, of the arrivals in seconds since
///                           the epoch.
/// @param[in] estimates      The times, UTC, of the estimates in seconds since
///                           the epoch.
/// @param[in] weights        The corresponding weights.
/// @result The origin time, UTC, in seconds since the epoch to add to the 
///         estimates.
[[nodiscard]] [[maybe_unused]]
double optimizeOriginTimeL1(const std::vector<double> &observations,
                            const std::vector<double> &estimates,
                            const std::vector<double> &weights)
{
    if (observations.empty())
    {
        throw std::invalid_argument("At least 1 observation required");
    }
    if (observations.size() != estimates.size())
    {
        throw std::invalid_argument("observations.size() != estimates.size()");
    }
    if (observations.size() != weights.size())
    {
        throw std::invalid_argument("observations.size() != weights.size()");
    }
    auto residuals = ::computeResiduals(observations, estimates);
    std::vector<std::pair<double, int>> workSpace;
    return ::weightedMedian(residuals, weights, workSpace); 
}

[[nodiscard]] [[maybe_unused]]
double optimizeOriginTimeL1(const ULocator::Origin &origin,
                            const std::vector<double> &estimates)
{
    auto [observations, weights] = ::originToObservationsAndWeights(origin);
    std::vector<double> reducedObservations;
    double reductionTime;
    ::reduceTimes<double>(observations, &reducedObservations, &reductionTime);
    return ::optimizeOriginTimeL1(reducedObservations, estimates, weights)
          + reductionTime;
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
    if (observations.empty())
    {   
        throw std::invalid_argument("At least 1 observation required");
    }
    if (observations.size() != estimates.size())
    {
        throw std::invalid_argument("observations.size() != estimates.size()");
    }
    if (observations.size() != weights.size())
    {
        throw std::invalid_argument("observations.size() != weights.size()");
    }
    constexpr double zero{0}; 
    double sumOfWeights = std::accumulate(weights.begin(), weights.end(), zero);
#ifndef NDEBUG
    assert(sumOfWeights > 0);
#endif
    auto residuals = ::computeResiduals(observations, estimates);
    return ::weightedMean(residuals, weights);
}

[[nodiscard]] [[maybe_unused]]
double optimizeOriginTimeLeastSquares(const ULocator::Origin &origin,
                                      const std::vector<double> &estimates)
{
    auto [observations, weights] = ::originToObservationsAndWeights(origin);
    std::vector<double> reducedObservations;
    double reductionTime;
    ::reduceTimes<double>(observations, &reducedObservations, &reductionTime);
    return ::optimizeOriginTimeLeastSquares(reducedObservations,
                                            estimates, 
                                            weights) + reductionTime;
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
    if (std::abs(p - 2) < 1.e-14)
    {
        return ::optimizeOriginTimeLeastSquares(observations, estimates, 
                                                weights);
    }
    else if (std::abs(p - 1) < 1.e-14)
    {
        return ::optimizeOriginTimeL1(observations, estimates, weights);
    }
    else if (p < 1)
    {
        throw std::invalid_argument("p must be greater than 1");
    }
    if (observations.empty())
    {
        throw std::invalid_argument("At least 1 observation required");
    }
    if (observations.size() != estimates.size())
    {
        throw std::invalid_argument("observations.size() != estimates.size()");
    }
    if (observations.size() != weights.size())
    {
        throw std::invalid_argument("observations.size() != weights.size()");
    }
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
                                maxIterations);
    return 0.5*(bracket.first + bracket.second);
}

[[nodiscard]] [[maybe_unused]]
double optimizeOriginTimeLp(const ULocator::Origin &origin,
                            const std::vector<double> &estimates,
                            const double p = 1.5,
                            const double timeWindow = 200,
                            const int maxIterations = 1500,
                            const double tolerance = 1.e-5)
{
    auto [observations, weights] = ::originToObservationsAndWeights(origin);
    std::vector<double> reducedObservations;
    double reductionTime;
    ::reduceTimes<double>(observations, &reducedObservations, &reductionTime);
    return ::optimizeOriginTimeLp(reducedObservations,
                                  estimates,
                                  weights,
                                  p,
                                  timeWindow,
                                  maxIterations,
                                  tolerance) + reductionTime;
}

}
#endif
