#ifndef ULOCATOR_OPTIMIZERS_PAGMO_PARTICLE_SWARM_HPP
#define ULOCATOR_OPTIMIZERS_PAGMO_PARTICLE_SWARM_HPP
#include <memory>
#include <vector>
#include <uLocator/optimizers/optimizer.hpp>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator
{
 class Origin;
}
namespace ULocator::Optimizers::Pagmo
{
class ParticleSwarm : public ULocator::Optimizers::IOptimizer
{
public:
    /// @brief Constructor.
    ParticleSwarm();
    /// @brief Constructor with a logger.
    explicit ParticleSwarm(std::shared_ptr<UMPS::Logging::ILog> &logger);
 
    /// @brief Destructor.
    virtual ~ParticleSwarm();

    /// @brief Sets the number of particles.
    /// @throws std::invalid_argument if nParticles is not positive.
    void setNumberOfParticles(int particles);
    /// @result The number of particles.
    [[nodiscard]] int getNumberOfParticles() const noexcept;
    /// @brief Sets the number of generations to evolve the particle
    ///        population.
    /// @throws std::invalid_argument if nGenerations is not positive.
    void setNumberOfGenerations(int nGenerations);
    /// @result The number of generations.
    [[nodiscard]] int getNumberOfGenerations() const noexcept;

    /// @brief Sets the origin time search window duration.
    /// @param[in] duration   The duration in seconds.  
    /// @throw std::invalid_argument if duration is not positive.
    void setOriginTimeSearchWindowDuration(double duration); 
    /// @result The origin time search window duration in seconds.
    [[nodiscard]] double getOriginTimeSearchWindowDuration() const noexcept;
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
    /// @brief Sets the depth search extent in meters.
    /// @param[in] extentInZ  The lower and upper extent to search in z in
    ///                       meters.  This is with respect to the free surface.
    void setExtentInZ(const std::pair<double, double> &extentInZ);
    /// @result The extent to search in z.  By default this is -2000, 65000.
    [[nodiscard]] std::pair<double, double> getExtentInZ() const;

    /// @brief Locates using the particle swarm method.
    /// @param[in] initialGuess     An initial guess for the origin.
    /// @param[in] locationProblem  The location problem.  
    /// @param[in] norm             The norm (misfit) in which to optimize.
    /// @note If an initial geuss is provided then only the depth will
    ///       be used.
    void locate(const ULocator::Origin &initialGuess,
                ULocator::Optimizers::IOptimizer::LocationProblem locationProblem,
                ULocator::Optimizers::IOptimizer::Norm norm = ULocator::Optimizers::IOptimizer::Norm::LeastSquares) final;
    /// @result The objective function corresponding to the best solution.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] double getOptimalObjectiveFunction() const;
    /// @result The number of objective function evaluations.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] int getNumberOfObjectiveFunctionEvaluations() const;
    /// @result True indicates the location phase was successful and there
    ///         is an origin. 
    [[nodiscard]] bool haveOrigin() const noexcept final;

    ParticleSwarm(const ParticleSwarm &) = delete;
    ParticleSwarm(ParticleSwarm &&) noexcept = delete;
private:
    class ParticleSwarmImpl;
    std::unique_ptr<ParticleSwarmImpl> pImpl;
};
}
#endif
