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
    /// @result The topographic elevation at a point.
    [[nodiscard]] virtual double operator()(double latitude, double longitude) const
    {
        return evaluate(latitude, longitude);
    }
};
}
#endif
