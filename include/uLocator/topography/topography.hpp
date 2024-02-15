#ifndef ULOCATOR_TOPOGRAPHY_TOPOGRAPHY_HPP
#define ULOCATOR_TOPOGRAPHY_TOPOGRAPHY_HPP
namespace ULocator::Topography
{
/// @class ITopography "topography.hpp" "uLocator/topography/topography.hpp"
/// @brief This is the abstract base class that defines a topography function.
///        This evaluates the elevation, measured positive up from sea-level
///        in meters, at a point x, y - i.e., this evaluates E(x,y).
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class ITopography
{
public:
    /// @brief Destructor.
    virtual ~ITopography() = default;
    /// @result True indicates the topography was set.
    [[nodiscard]] virtual bool haveTopography() const noexcept = 0;
    /// @result The topographic elevation at a point.
    [[nodiscard]] virtual double evaluate(double, double) const = 0;
    /// @result The elevation at a point along with itsderivatives.
    [[nodiscard]] virtual double evaluate(double, double, double *, double *) const = 0;
    /// @result The minimum and maximum elevation with respect to sea-level
    ///         in meters.
    [[nodiscard]] virtual std::pair<double, double> getMinimumAndMaximumElevation() const = 0;
    /// @result The topographic elevation at a point.
    /// @param[in] x   The x ordinate at which to evaluate E(x,y).
    /// @param[in] y   The y ordinate at which to evaluate E(x,y).
    [[nodiscard]] virtual double operator()(double x, double y) const
    {
        return evaluate(x, y);
    }
    /// @result The topographic elevation at a point along with the derivatives.
    /// @param[in] x   The x ordinate at which to evaluate E(x,y).
    /// @param[in] y   The y ordinate at which to evaluate E(x,y).
    /// @param[out] dElevationDx  If this is not null then then is dE/dx 
    ///                           evaluated at (x, y).  This is unitless.
    /// @param[out] dElevationDy  If this is not null then then is dE/dy 
    ///                           evaluated at (x, y).  This is unitless.
    [[nodiscard]] virtual double operator()(double x, double y,
                                            double *dElevationDx, double *dElevationDy) const
    {
        return evaluate(x, y, dElevationDx, dElevationDy);
    }
};
}
#endif
