#ifndef ULOCATOR_PAGMO_FITNESS_FUNCTIONS_HPP
#define ULOCATOR_PAGMO_FITNESS_FUNCTIONS_HPP
#include <string>
#include <initializer_list>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include <pagmo/problem.hpp>
#include <pagmo/types.hpp>
#include "../objectiveFunctions.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/station.hpp"

namespace
{

template<typename T>
[[nodiscard]] 
T leastSquares(const std::vector<T> &weights,
               const std::vector<T> &observations,
               const std::vector<T> &estimates)
{
     double sumSquaredOfResiduals{0};
     for (size_t i = 0; i < observations.size(); ++i)
     {
         auto weightedResidual
            = weights[i]*(observations[i] - estimates[i]);
         sumSquaredOfResiduals = sumSquaredOfResiduals
                               + weightedResidual*weightedResidual;
     }
     return sumSquaredOfResiduals;
}

template<typename T>
[[nodiscard]]
T l1(const std::vector<T> &weights,
     const std::vector<T> &observations,
     const std::vector<T> &estimates)
{
     double sumAbsoluteResiduals{0};
     for (size_t i = 0; i < observations.size(); ++i)
     {
         auto weightedResidual
            = weights[i]*std::abs(observations[i] - estimates[i]);
         sumAbsoluteResiduals = sumAbsoluteResiduals + weightedResidual;
     }
     return sumAbsoluteResiduals;
}

//----------------------------------------------------------------------------//

struct PagmoProblem3DAndTime : public ::Problem3DAndTime
{
    /// @result The objective function and the inequality constraint that 
    ///         requires the source be below the topography. 
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {
#ifndef NDEBUG
        assert(dv.size() == 4);
#endif
        // Fitness functions are to be minimized hence we make it negative
        auto fitness =-::Problem3DAndTime::logLikelihood(dv.size(), dv.data());
        return pagmo::vector_double {fitness};
    }
    /// @result The number of equality constraints
    pagmo::vector_double::size_type get_eic() const
    {
        return 0;
    }
    /// @result The number of inequality constraints.
    pagmo::vector_double::size_type get_nic() const
    {
        return 0;
    }
    /// @result The gradient {dT/dt0, dT/dx, dT/dy, dT/dz} of the objective
    ///         function.
    pagmo::vector_double gradient(const pagmo::vector_double &dv) const
    {
#ifndef NDEBUG
        assert(dv.size() == 4);
#endif
        pagmo::vector_double gradient(4);
        auto fitness =-::Problem3DAndTime::logLikelihoodAndGradient(
            dv.size(), dv.data(), gradient.data());
        std::transform(gradient.begin(), gradient.end(), gradient.begin(),
                       [](const double gi)
                       {
                           return -gi;
                       });
        return gradient;
    }
    /// @result The hard bounds on the search region.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
};

//----------------------------------------------------------------------------//

struct PagmoProblem2DAndTimeAndDepthAtFreeSurface : public
    Problem2DAndTimeAndDepthAtFreeSurface
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {   
#ifndef NDEBUG
        assert(dv.size() == 3); 
#endif
        auto fitness
            =-::Problem2DAndTimeAndDepthAtFreeSurface::logLikelihood(dv.size(),
                                                                     dv.data());
        return pagmo::vector_double {fitness};
    }   
    /// @result The number of equality constraints
    pagmo::vector_double::size_type get_eic() const
    {
        return 0;
    }
    /// @result The number of inequality constraints.
    pagmo::vector_double::size_type get_nic() const
    {
        return 0;
    }
    /// @result The gradient {dT/dt0, dT/dx, dT/dy} of the objective
    ///         function.
    pagmo::vector_double gradient(const pagmo::vector_double &dv) const
    {
#ifndef NDEBUG
        assert(dv.size() == 3); 
#endif
        pagmo::vector_double gradient(3);
        auto fitness
            =-::Problem2DAndTimeAndDepthAtFreeSurface::logLikelihoodAndGradient(
                 dv.size(), dv.data(), gradient.data());
        std::transform(gradient.begin(), gradient.end(), gradient.begin(),
                       [](const double gi) 
                       {
                           return -gi;
                       });
        return gradient;
    }
    /// @result The hard bounds on the search region.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
};

//----------------------------------------------------------------------------//

struct PagmoProblem2DAndTimeAndFixedDepth :
    public Problem2DAndTimeAndFixedDepth
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {   
#ifndef NDEBUG
        assert(dv.size() == 3); 
#endif
        auto fitness
            =-::Problem2DAndTimeAndFixedDepth::logLikelihood(dv.size(),
                                                             dv.data());
        return pagmo::vector_double {fitness};
    }   
    /// @result The number of equality constraints
    pagmo::vector_double::size_type get_eic() const
    {   
        return 0;
    }   
    /// @result The number of inequality constraints.
    pagmo::vector_double::size_type get_nic() const
    {   
        return 0;
    }
    /// @result The gradient {dT/dt0, dT/dx, dT/dy} of the objective
    ///         function.
    pagmo::vector_double gradient(const pagmo::vector_double &dv) const
    {
#ifndef NDEBUG
        assert(dv.size() == 3);
#endif
        pagmo::vector_double gradient(3);
        auto fitness
            =-::Problem2DAndTimeAndFixedDepth::logLikelihoodAndGradient(
                 dv.size(), dv.data(), gradient.data());
        std::transform(gradient.begin(), gradient.end(), gradient.begin(),
                       [](const double gi)
                       {
                           return -gi;
                       });
        return gradient;
    }
    /// @result The hard bounds on the search region.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    } 
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
};

}
#endif
