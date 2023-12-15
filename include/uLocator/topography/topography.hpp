#ifndef ULOCATOR_TOPOGRAPHY_TOPOGRAPHY_HPP
#define ULOCATOR_TOPOGRAPHY_TOPOGRAPHY_HPP
namespace ULocator::Topography
{
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
    [[nodiscard]] virtual double operator()(double x, double y) const
    {
        return evaluate(x, y);
    }
    [[nodiscard]] virtual double operator()(double x, double y,
                                            double *dElevationDx, double *dElevationDy) const
    {
        return evaluate(x, y, dElevationDx, dElevationDy);
    }
};
}
#endif
