#ifndef ULOCATOR_OPTIMIZERS_NLOPT_STOGO_HPP
#define ULOCATOR_OPTIMIZERS_NLOPT_STOGO_HPP
#include <memory>
#include <uLocator/optimizers/optimizer.hpp>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator::Optimizers::NLOpt
{
/// @class StochasticGradientOptimization "dividedRectangles.hpp" "uLocator/optimizers/nlopt/dividedRectangles.hpp"
/// @brief The Stochastic Gradient Optimization (StoGo) is another excellent
///        global-optimization scheme that iteratively refines a search space
///        with smaller and smaller (hyper)-rectangles.  However, unlike DIRECT
///        it uses a local BFGS solver to search those hyper-rectangles.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class StochasticGradientOptimization : public ULocator::Optimizers::IOptimizer
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    StochasticGradientOptimization();
    /// @brief Constructor with a logger.
    explicit StochasticGradientOptimization(std::shared_ptr<UMPS::Logging::ILog> &logger);
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
    void setLocationTolerance(double tolerance);
    /// @result The location tolerance in meters.
    [[nodiscard]] double getLocationTolerance() const noexcept;

    /// @brief The algorithm will terminate if the event's origin times
    ///        are moving less than this tolerance.
    /// @param[in] tolerance  The tolerance in seconds.
    void setOriginTimeTolerance(double tolerance);
    /// @result The origin time tolerance in seconds.
    [[nodiscard]] double getOriginTimeTolerance() const noexcept;

    /// @brief Enables randomization in the optimization.  This seems to help.
    void enableRandomization() noexcept;
    /// @brief Disables randomization in the optimization.
    void disableRandomization() noexcept;
    /// @result True indicates there will be some stochasticity in the search.
    /// @note By default this is true.
    [[nodiscard]] bool randomize() const noexcept;

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
    ~StochasticGradientOptimization() override;
    /// @}

    StochasticGradientOptimization& operator=(const StochasticGradientOptimization &) = delete;
    StochasticGradientOptimization(const StochasticGradientOptimization &) = delete;
private:
    class StochasticGradientOptimizationImpl;
    std::unique_ptr<StochasticGradientOptimizationImpl> pImpl;
};
}
#endif
