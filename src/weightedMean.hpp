#ifndef WEIGHTED_MEAN
#define WEIGHTED_MEAN
#include <algorithm>
#include <limits>
#include <numeric>

namespace
{

double weightedMean(const std::vector<double> &residuals,
                    const std::vector<double> &weights)
{
    constexpr double zero{0};
    auto numerator = std::inner_product(residuals.begin(), residuals.end(),
                                        weights.begin(), zero);
    auto denominator = std::accumulate(weights.begin(), weights.end(), zero);
    denominator = std::max(denominator, std::numeric_limits<double>::epsilon());
    return numerator/denominator; 
}

}
#endif
