#ifndef ULOCATOR_OPTIMIZERS_NLOPT_DIRECT_HPP
#define ULOCATOR_OPTIMIZERS_NLOPT_DIRECT_HPP
#include <memory>
#include <uLocator/optimizers/optimizer.hpp>
#include <umps/logging/log.hpp>
namespace ULocator::Optimizers::NLOpt
{
/// @class DividedRectangles "dividedRectangles.hpp" "uLocator/optimizers/nlopt/dividedRectangles.hpp"
/// @brief The DIvided RECTangles (DIRECT) is an excellent global-optimization
///        scheme that iteratively refines a search space with smaller 
///        and smaller (hyper)-rectangles.  It's very similar to the Oct-tree 
///        search in NonLinLoc.  I typically use DIRECT to identify an initial
///        solution for a local-optimization algorithm.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class DividedRectangles : public ULocator::Optimizers::IOptimizer
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    DividedRectangles();
    /// @brief Constructor with a logger.
    explicit DividedRectangles(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @}

    /// @name Optimization Options
    /// @{

    /// @brief Sets the maximum number of objective function evaluations.
    /// @param[in] nEvaluations  The maximum number
    void setMaximumNumberOfObjectiveFunctionEvaluations(int nEvaluations);
    /// @result The model terminates if this number of objective function
    ///         evaluations is exceeded.
    /// @note If the model pert
    [[nodiscard]] int getMaximumNumberOfObjectiveFunctionEvaluations() const noexcept;

    /// @brief The algorithm will terminate if the event's hypocentral positions
    ///        are moving less than this tolerance. 
    /// @param[in] tolerance  The tolerance in meters.
    /// @bug I don't think NLOpt uses this.
    void setLocationTolerance(double tolerance);
    /// @result The location tolerance in meters.
    /// @bug I don't think NLOpt uses this.
    [[nodiscard]] double getLocationTolerance() const noexcept;

    /// @brief The algorithm will terminate if the event's origin times
    ///        are moving less than this tolerance.
    /// @param[in] tolerance  The tolerance in seconds.
    /// @bug I don't think NLOpt uses this.
    void setOriginTimeTolerance(double tolerance);
    /// @result The origin time tolerance in seconds.
    /// @bug I don't think NLOpt uses this.
    [[nodiscard]] double getOriginTimeTolerance() const noexcept;

    /// @brief Enables normalization of the simple box-boundaries.
    void enableNormalization() noexcept;
    /// @brief Disables normalization of the simple box-boundaries.
    void disableNormalization() noexcept;
    /// @result True indicates the simple box-boundaries will be normalized
    ///         prior to searching.
    /// @note By default normalization is set.
    [[nodiscard]] bool normalize() const noexcept;

    /// @brief Enables biasing towards local search.  This implements the
    ///        "A locally-biased form of the DIRECT algorithm" by 
    ///        Gablonsky and Kelley and is more efficient when the 
    ///        objective function does not have too many local minima.
    void enableLocallyBias() noexcept;
    /// @brief Disables biasing towards local search.  This implements
    ///        the classic "Lipshcitzian optimization without Lipschitz
    ///        constant" method of Jones et al.
    void disableLocallyBias() noexcept;
    /// @result Indicates whether or not to use the locally biased DIRECT
    ///         algorithm.
    [[nodiscard]] bool locallyBias() const noexcept;

    /// @brief Sets the origin time search window duration.
    /// @param[in] duration   The duration in seconds.  
    /// @throw std::invalid_argument if duration is not positive.
    void setOriginTimeSearchWindowDuration(double duration); 
    /// @result The origin time search window duration in seconds.
    [[nodiscard]] double getOriginTimeSearchWindowDuration() const noexcept;
    /// @}


    /// @brief Locates using the DIRECT method.
    /// @param[in] initialGuess     An initial guess for the origin.
    /// @param[in] locationProblem  The location problem.  
    /// @param[in] norm             The norm (misfit) in which to optimize.
    /// @note If an initial geuss is provided then only the depth will
    ///       be used.
    void locate(const ULocator::Origin &initialGuess,
                ULocator::Optimizers::IOptimizer::LocationProblem locationProblem,
                ULocator::Optimizers::IOptimizer::Norm norm = ULocator::Optimizers::IOptimizer::Norm::LeastSquares) final;
    /// @result True indicates the origin is available.
    [[nodiscard]] bool haveOrigin() const noexcept final;
    /// @result The number of objective function evaluations during
    ///         optimization.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] int getNumberOfObjectiveFunctionEvaluations() const;
    /// @result The number of gradient evaluations during optimization.
    /// @note This will be 0 because DIRECT Is derivative-free.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] int getNumberOfGradientEvaluations() const;
    /// @result The objective function at the optimum location.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] double getOptimalObjectiveFunction() const;

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~DividedRectangles() override;
    /// @}

    DividedRectangles& operator=(const DividedRectangles &) = delete;
    DividedRectangles(const DividedRectangles &) = delete;
private:
    class DividedRectanglesImpl;
    std::unique_ptr<DividedRectanglesImpl> pImpl;
};
}
#endif
