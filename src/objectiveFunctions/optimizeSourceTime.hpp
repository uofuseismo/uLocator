#ifndef OPTIMIZE_SOURCE_TIME_HPP
#define OPTIMIZE_SOURCE_TIME_HPP
#include <vector>
#include <algorithm>
#include <numeric>
#ifndef NDEBUG
#include <cassert>
#endif
#include "../weightedMedian.hpp"
namespace
{
// Least-squares optimization -> best source time is mean
[[maybe_unused]]
[[nodiscard]]
double optimizeSourceTimeLeastSquares(const std::vector<double> &residuals,
                                      const std::vector<double> &weights,
                                      const double sumOfWeights)
{
#ifndef NDEBUG
    assert(residuals.size() == weights.size());
    assert(!weights.empty());
    assert(sumOfWeights > 0);
#endif
    double sumOfWeightedResiduals
        = std::inner_product(weights.begin(), weights.end(),
                             residuals.begin(), 0.0);
    return sumOfWeightedResiduals/sumOfWeights;
}
// L1 optimization -> best source time is weighted median
double optimizeSourceTimeL1(const std::vector<double> &residuals,
                            const std::vector<double> &weights,
                            std::vector<std::pair<double, int>> &workSpace)
{
#ifndef NDEBUG
    assert(residuals.size() == weights.size());
    assert(!weights.empty());
#endif
    return ::weightedMedian(residuals, weights, workSpace);
}
[[maybe_unused]]
[[nodiscard]]
double optimizeSourceTimeL1(const std::vector<double> &residuals,
                            const std::vector<double> &weights)
{
    std::vector<std::pair<double, int>> workSpace(weights.size());
    return optimizeSourceTimeL1(residuals, weights, workSpace);
}
}
#endif
