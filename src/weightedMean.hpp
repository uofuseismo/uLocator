#ifndef ULOCATOR_PRIVATE_WEIGHTED_MEAN
#define ULOCATOR_PRIVATE_WEIGHTED_MEAN
#include <algorithm>
#include <limits>
#include <numeric>
#ifndef NDEBUG
#include <cassert>
#endif
namespace
{
/// @param[in] residuals  The residuals in seconds.
/// @param[in] weights    The weights (in 1/seconds).
/// @result The weighted mean of a residual vector.
double weightedMean(const std::vector<double> &residuals,
                    const std::vector<double> &weights)
{
#ifndef NDEBUG
    assert(!residuals.empty());
    assert(!weights.empty());
#endif
    constexpr double zero{0};
    auto numerator = std::inner_product(residuals.begin(), residuals.end(),
                                        weights.begin(), zero);
    auto denominator = std::accumulate(weights.begin(), weights.end(), zero);
    denominator = std::max(denominator, std::numeric_limits<double>::epsilon());
    return numerator/denominator; 
}

}
#endif
