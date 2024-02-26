#ifndef ULOCATOR_OPTIMIZERS_NLOPT_MMA_HPP
#define ULOCATOR_OPTIMIZERS_NLOPT_MMA_HPP
#include <memory>
#include <uLocator/optimizers/optimizer.hpp>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator::Optimizers::NLOpt
{
/// @class MethodOfMovingAsymptotes "methodOfMovingAsymtotes.hpp" "uLocator/optimizers/nlopt/methodOfMovingAsymtotes.hpp"
/// @brief The Method of Moving Asymtotes (MMA) is a local, gradient-based
///        optimizer.  If it is provided a feasible starting condition then it
///        is supposed to converge to a corresponding local optimum.
///        I typically use MMA to fine-tune earthquake locations generated
///        by a globally search method like the DIRECT algorithm. 
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class MethodOfMovingAsymptotes : public ULocator::Optimizers::IOptimizer
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    MethodOfMovingAsymptotes();
    /// @brief Constructor with a logger.
    explicit MethodOfMovingAsymptotes(std::shared_ptr<UMPS::Logging::ILog> &logger);
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

    /// @brief Sets the origin time search window duration.
    /// @param[in] duration   The duration in seconds.  
    /// @throw std::invalid_argument if duration is not positive.
    void setOriginTimeSearchWindowDuration(double duration); 
    /// @result The origin time search window duration in seconds.
    [[nodiscard]] double getOriginTimeSearchWindowDuration() const noexcept;
    /// @}

    /// @name Refined Search Region
    /// @{

    /// @brief Sets the search extent in X.
    /// @param[in] extentInX   The lower and upper extent to search in x in
    ///                        meters.  This should be a subset of the region's
    ///                        extent.
    void setExtentInX(const std::pair<double, double> &extentInX);
    /// @result The extent to search in x.  By default this is the region's
    ///         extent.
    [[nodiscard]] std::pair<double, double> getExtentInX() const;
    /// @brief Sets the search extent in Y.
    /// @param[in] extentInY   The lower and upper extent to search in y in
    ///                        meters.  This should be a subset of the region's
    ///                        extent.
    void setExtentInY(const std::pair<double, double> &extentInY);
    /// @result The extent to search in y.  By default this is the region's
    ///         extent.
    [[nodiscard]] std::pair<double, double> getExtentInY() const;
    /// @}


    /// @brief Locates using the MMA method.
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
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] int getNumberOfGradientEvaluations() const;
    /// @result The objective function at the optimum location.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] double getOptimalObjectiveFunction() const;

    /// @name Destructors
    /// @{

    /// @brief Destructor.
    ~MethodOfMovingAsymptotes() override;
    /// @}

    MethodOfMovingAsymptotes& operator=(const MethodOfMovingAsymptotes &) = delete;
    MethodOfMovingAsymptotes(const MethodOfMovingAsymptotes &) = delete;
private:
    class MethodOfMovingAsymptotesImpl;
    std::unique_ptr<MethodOfMovingAsymptotesImpl> pImpl;
};
}
#endif
