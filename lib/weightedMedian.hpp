#ifndef ULOCATOR_WEIGHTED_MEDIAN_HPP
#define ULOCATOR_WEIGHTED_MEDIAN_HPP
#include <vector>
#include <numeric>
#ifndef NDEBUG
#include <cassert>
#endif
#include <algorithm>
namespace
{

/// @param[in] x  The vector residuals of which to find the weighted.
///               median.  This has units of seconds.
/// @param[in] weights  The corresponding weights.  This has units of
///                     1/s.
/// @param[in,out] workSpace  A workspace array.
/// @result The weighted median of the residual vector with units of
///         seconds.
[[nodiscard]]
double weightedMedian(const std::vector<double> &x,
                      const std::vector<double> &weights,
                      std::vector<std::pair<double, int>> &workSpace)
{
    constexpr double tol = std::numeric_limits<double>::epsilon()*100;
#ifndef NDEBUG
    assert(x.size() == weights.size());
#endif
    if (x.empty()){return 0;}
    if (x.size() == 1){return x[0];}
    if (workSpace.size() != x.size())
    {
        workSpace.resize(x.size());
    }
    auto nResiduals = static_cast<int> (x.size());
    // Copy x/weights into workspace
    bool weightsEqual{true};
    double weight0 = weights[0];
    double sumOfWeights = 0;
    for (int i = 0; i < nResiduals; ++i)
    {
#ifndef NDEBUG
        assert(weights[i] > 0);
#endif
        workSpace[i] = std::pair{x[i], i};
        sumOfWeights = sumOfWeights + weights[i];
        if (std::abs(weight0 - weights[i]) > tol)
        {
            weightsEqual = false;
        }
    }
    double sumOfWeightsInverse = 1./sumOfWeights;
    // (Arg)sort
    std::sort(workSpace.begin(), workSpace.end(),
             [](const std::pair<double, int> &lhs,
                const std::pair<double, int> &rhs)
             {
                 return lhs.first < rhs.first;
             });
    // Weights are equal -> this is just a median
    if (weightsEqual)
    {
        auto iHalf = static_cast<int> ((nResiduals + 1)/2) - 1;
        auto jHalf = workSpace[iHalf].second;
        if (nResiduals%2 == 1)
        {
            return x[jHalf];
        }
        else
        {
            auto j0 = jHalf;
            auto j1 = workSpace[iHalf + 1].second;
            return 0.5*(x[j0] + x[j1]);
        }
    }
    // Locate index where cumulative weight passes 1/2
    double wi2 = 0;
    int iHalf = nResiduals - 1; // Default is last
    for (int i = 0; i < nResiduals; ++i)
    {
        auto j = workSpace[i].second;
        auto wi = weights[j]*sumOfWeightsInverse;
        wi2 = wi2 + wi;
        workSpace[i] = std::pair {wi2 - 0.5*wi, j};
        if (workSpace[i].first >= 0.5)
        {
            iHalf = i;
            break;
        }
    }
    // Two special cases are the ends
    if (iHalf == 0 || iHalf == nResiduals - 1)
    {   
        auto jHalf = workSpace[iHalf].second;
        return x[jHalf];
    }   
    // Linearly interpolate the median
    auto jHalfUpper = workSpace.at(iHalf).second;
    auto jHalfLower = workSpace.at(iHalf - 1).second;
    auto x1 = x.at(jHalfUpper);
    auto x0 = x.at(jHalfLower);
    auto y1 = workSpace.at(iHalf).first;
    auto y0 = workSpace.at(iHalf - 1).first;
    return x[jHalfUpper] + ((0.5 - y1)/(y1 - y0))*(x1 - x0);
}

}
#endif
