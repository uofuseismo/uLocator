#ifndef ULOCATOR_TOPOGRAPHY_CONSTANT_HPP
#define ULOCATOR_TOPOGRAPHY_CONSTANT_HPP
#include <memory>
#include "uLocator/topography/topography.hpp"
namespace ULocator::Topography
{
/// @class Constant "constant.hpp" "uLocator/topography/constant.hpp"
/// @brief This is a constant topography function which is generally not useful
///        but fairly easy to use.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Constant : public ITopography
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Constant();
    /// @brief Copy constructor.
    /// @param[in] topography   The class from which to initialize this class.
    Constant(const Constant &topography);
    /// @brief Move constructor.
    /// @brief
    Constant(Constant &&topography) noexcept;
    /// @}

    /// @brief Sets the constant topography.
    /// @param[in] elevation  The constant topographic elevation in meters
    ///                       with respect to sea-level.  Here, elevation
    ///                       increases up into the sky.
    void set(double elevation) noexcept;
    /// @result True indicates the topography is set.
    [[nodiscard]] bool haveTopography() const noexcept final;
    /// @result The elevation with respect to sea-level in meters.  Note,
    /// @throws std::runtime_error if \c haveTopography() is false.
    [[nodiscard]] double evaluate(double, double) const final;
    /// @result The elevation with respect to sea-level in meters.  Note,
    /// @throws std::runtime_error if \c haveTopography() is false.
    [[nodiscard]] double evaluate(double, double,
                                  double *dElevationDx, double *dElevationDy) const final;
    /// @result The minimum and maximum elevation with respect to sea-level
    ///         in meters.
    [[nodiscard]] std::pair<double, double> getMinimumAndMaximumElevation() const final;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] topography  The topography to copy to this.
    /// @result A deep copy of the topography.
    Constant& operator=(const Constant &topography);
    /// @brief Move assignment.
    /// @param[in,out] topography  The topography whose memory will be moved 
    ///                            to this.  On exit, topography's behavior
    ///                            is undefined.
    /// @result The memory from topography moved to this.
    Constant& operator=(Constant &&topography) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    virtual ~Constant();
    /// @}
private:
    class ConstantImpl; 
    std::unique_ptr<ConstantImpl> pImpl;
};
}
#endif
