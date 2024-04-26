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
#include "optimizers/objectiveFunctions.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/station.hpp"

namespace
{

//----------------------------------------------------------------------------//

struct PagmoProblem3DAndTime : public ::Problem3DAndTime
{
    /// @result The objective function and the inequality constraint that 
    ///         requires the source be below the topography. 
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {
        double fitness = std::numeric_limits<double>::max();
#ifndef NDEBUG
        assert(static_cast<int> (dv.size()) == mParameters);
#endif
        try
        {
            // Fitness functions are to be minimized hence we make it negative
            fitness =-::Problem3DAndTime::logLikelihood(dv.size(), dv.data());
        }
        catch (const std::exception &e)
        {
            auto errorMessage = "Problem for 3D source with (t,x,y,z) = ("
                              + std::to_string(dv.at(0)) + ","
                              + std::to_string(dv.at(1)) + ","
                              + std::to_string(dv.at(2)) + ","
                              + std::to_string(dv.at(3))
                              + ").  Failed with: "  + std::string {e.what()};
            std::cerr << errorMessage << std::endl;
        }
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
        assert(dv.size() == mParameters);
#endif
        pagmo::vector_double gradient(mParameters);
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
    /// @brief Sets the search boundaries.
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 4");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 4");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }
    [[nodiscard]] pagmo::vector_double createInitialGuess(
        const ULocator::Origin &initialGuess,
        const ULocator::Position::IGeographicRegion &region)
    {
        if (mLowerBounds.empty())
        {
            throw std::runtime_error("Search boundaries not set");
        }
        // Initial guess
        pagmo::vector_double xStartingLocation(mParameters);
        for (int i = 0; i < mParameters; ++i)
        {
            xStartingLocation.at(i) = 0.5*(mLowerBounds.at(i)
                                         + mUpperBounds.at(i));
        }
        if (initialGuess.haveTime())
        {
            auto time = initialGuess.getTime() - mReductionTime;
            if (time >= mLowerBounds.at(0) &&
                time <= mUpperBounds.at(0))
            {
                xStartingLocation.at(0) = time;
            }
        }
        if (initialGuess.haveEpicenter())
        {
            auto epicenter = initialGuess.getEpicenter();
            auto [xGuess, yGuess] = region.geographicToLocalCoordinates(
                epicenter.getLatitude(), epicenter.getLongitude());
            if (xGuess >= mLowerBounds.at(1) &&
                xGuess <= mUpperBounds.at(1))
            {
                xStartingLocation.at(1) = xGuess;
            }
            if (yGuess >= mLowerBounds.at(2) &&
                yGuess <= mUpperBounds.at(2))
            {
                xStartingLocation.at(2) = yGuess;
            }
        }
        if (initialGuess.haveDepth())
        {
            auto depth = initialGuess.getDepth();
            if (depth >= mLowerBounds.at(3) &&
                depth <= mUpperBounds.at(3))
            {
                xStartingLocation.at(3) = depth;
            }
        }
        return xStartingLocation;
    }
    /// @result The location corresponding to a Pagmo result
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const pagmo::vector_double &xLocation,
                     const ULocator::Position::IGeographicRegion &region) const
    {
        ULocator::Origin origin;
        // Extract the origin information
        origin.setTime(mReductionTime + xLocation.at(0));
        double xSource{xLocation.at(1)};
        double ySource{xLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region.localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(xLocation.at(3));
        return origin;
    }
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
    int mParameters{4};
};

//----------------------------------------------------------------------------//

struct PagmoProblem2DAndTimeAndDepthAtFreeSurface : public
    Problem2DAndTimeAndDepthAtFreeSurface
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {
        double fitness = std::numeric_limits<double>::max();
#ifndef NDEBUG
        assert(dv.size() == mParameters); 
#endif
        try
        {
            fitness
                =-::Problem2DAndTimeAndDepthAtFreeSurface::logLikelihood(
                      dv.size(), dv.data());
        }
        catch (const std::exception &e)
        {
            auto errorMessage = "Problem for source at f.s. (t,x,y) = ("
                              + std::to_string(dv.at(0)) + ","
                              + std::to_string(dv.at(1)) + ","
                              + std::to_string(dv.at(2)) 
                              + ").  Failed with: "  + std::string {e.what()};
            std::cerr << errorMessage << std::endl;
        }
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
        assert(dv.size() == mParameters); 
#endif
        pagmo::vector_double gradient(mParameters);
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
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 3");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 3");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const pagmo::vector_double &xLocation,
                     const ULocator::Position::IGeographicRegion &region) const
    {   
        ULocator::Origin origin; 
        // Extract the origin information
        origin.setTime(mReductionTime + xLocation.at(0));
        double xSource{xLocation.at(1)};
        double ySource{xLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region.localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        double sourceDepth =-mTopography->evaluate(xSource, ySource);
        origin.setDepth(sourceDepth);
        return origin;
    }   
    /// @result The hard bounds on the search region.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    }
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
    int mParameters{3};
};

//----------------------------------------------------------------------------//

struct PagmoProblem2DAndTimeAndFixedDepth :
    public Problem2DAndTimeAndFixedDepth
{
    pagmo::vector_double fitness(const pagmo::vector_double &dv) const
    {   
        double fitness = std::numeric_limits<double>::max();
#ifndef NDEBUG
        assert(static_cast<int> (dv.size()) == mParameters); 
#endif
        try
        {
            fitness
               =-::Problem2DAndTimeAndFixedDepth::logLikelihood(dv.size(),
                                                                dv.data());
        }
        catch (const std::exception &e)
        {
            auto errorMessage = "Problem for source with (t,x,y,z) = ("
                              + std::to_string(dv.at(0)) + ","
                              + std::to_string(dv.at(1)) + ","
                              + std::to_string(dv.at(2)) + ","
                              + std::to_string(mDepth)
                              + ").  Failed with: "  + std::string {e.what()};
            std::cerr << errorMessage << std::endl;
        }
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
        assert(static_cast<int> (dv.size()) == mParameters);
#endif
        pagmo::vector_double gradient(mParameters);
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
    [[nodiscard]] ULocator::Origin
    locationToOrigin(const pagmo::vector_double &xLocation,
                     const ULocator::Position::IGeographicRegion &region) const
    {
        ULocator::Origin origin;
        // Extract the origin information
        origin.setTime(mReductionTime + xLocation.at(0));
        double xSource{xLocation.at(1)};
        double ySource{xLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region.localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(mDepth);
        return origin;
    }
    /// @result The hard bounds on the search region.
    std::pair<pagmo::vector_double, pagmo::vector_double> get_bounds() const
    {
        return std::pair {mLowerBounds, mUpperBounds};
    } 
    void setSearchBoundaries(
        const std::vector<double> &lowerBoundaries,
        const std::vector<double> &upperBoundaries)
    {
        if (static_cast<int> (lowerBoundaries.size()) != mParameters)
        {
            throw std::invalid_argument("lowerBoundaries size != 3");
        }
        if (lowerBoundaries.size() != upperBoundaries.size())
        {
            throw std::invalid_argument("upperBoundaries size != 3");
        }
        for (int i = 0; i < static_cast<int> (lowerBoundaries.size()); ++i)
        {
            if (lowerBoundaries[i] >= upperBoundaries[i])
            {
                throw std::invalid_argument("l[i] >= u[i]");
            }
        }
        mLowerBounds = lowerBoundaries;
        mUpperBounds = upperBoundaries;
    }
    pagmo::vector_double mLowerBounds;
    pagmo::vector_double mUpperBounds;
    int mParameters{3};
};

}
#endif
