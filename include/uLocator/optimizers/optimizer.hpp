#ifndef ULOCATOR_OPTIMIZERS_OPTIMIZER_HPP
#define ULOCATOR_OPTIMIZERS_OPTIMIZER_HPP
#include <memory>
#include <vector>
#include <umps/logging/log.hpp>
namespace ULocator
{
 class Arrival;
 class Origin;
 class TravelTimeCalculatorMap;
 namespace Position 
 {
  class IGeographicRegion;
 }
 namespace Topography
 {
  class ITopography;
 }
}
namespace ULocator::Optimizers
{
/// @brief This some common/core functionality for an optimizer.
class IOptimizer
{
public:
    /// @brief Defines the location problem to solve.
    enum class LocationProblem
    {
        ThreeDimensionsAndTime,  /*!< This is the most general problem which
                                      involves solving for the origin time
                                      and earthquake hypocenter. */
        FixedToFreeSurfaceAndTime, /*!< This solves for the origin time
                                        and epicenter while restricting the
                                        event to the free surface.  This is
                                        useful for quarry blasts. */
        FixedDepthAndTime /*!< This solves for the origin time and epicenter
                               while fixing the event to a given depth.  This
                               can be useful during an initial, global search,
                               location phase, then refining with
                               a local search using ThreeDimensionsAndtime. */
    };
    /// @brief Defines
    enum class Norm
    {
        LeastSquares, /*!< Standard least-squares */
        L1, /*!< The L1 norm.  */
        Lp  /*!< An L_p norm.  Note, you should use the specialized 
                 LeastSquares for Lp with p = 2 and L1 with p = 1. */ 
    };
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    IOptimizer();
    /// @brief Constructor with a logger.
    explicit IOptimizer(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @}

    /// @brief Locates the event.
    /// @param[in] locationProblem  Defines the location heuristic.
    /// @param[in] initialGuess     This provides an initial depth and/or
    ///                             epicentral estimate.
    /// @param[in] norm             Defines the norm to optimize.
    virtual void locate(const ULocator::Origin &initialGuess,
                        LocationProblem locationProblem,
                        Norm norm = Norm::LeastSquares) = 0;
    /// @brief Locates the event but assumes nothing about the event location.
    /// @param[in] locationProblem  Defines the location heuristic.
    /// @param[in] norm             Defines the norm to optimize.
    virtual void locate(LocationProblem locationProblem,
                        Norm norm = Norm::LeastSquares);
    /// @brief Locates an event afixed to the free surface.
    virtual void locateEventAtFreeSurface(Norm norm = Norm::LeastSquares);
    /// @brief Locates an event at a fixed depth.
    virtual void locateEventWithFixedDepth(double depth,
                                           Norm norm = Norm::LeastSquares);
    /// @brief Evaluates the objective function at a given location.
    virtual double evaluateObjectiveFunction(const ULocator::Origin &origin,
                                             Norm norm = Norm::LeastSquares) const;
    /// @name Region
    /// @{
 
    /// @param[in] region  This defines the search region as a simple rectangle 
    ///                    and the mapping from this simple rectangle to and
    ///                    from geographic coordinates.
    virtual void setGeographicRegion(const Position::IGeographicRegion &region);
    /// @result The geographic region.
    /// @throws std::runtime_error \c haveGeographicRegion() is false.
    [[nodiscard]] virtual std::unique_ptr<Position::IGeographicRegion> getGeographicRegion() const;
    /// @result True indicates the region was set.
    [[nodiscard]] virtual bool haveGeographicRegion() const noexcept;
    /// @}

    /// @name Travel Time Calculator
    /// @{

    virtual void setTravelTimeCalculatorMap(std::unique_ptr<ULocator::TravelTimeCalculatorMap> &&calculatorMap);
    virtual bool haveTravelTimeCalculatorMap() const noexcept; 
    [[nodiscard]] virtual const ULocator::TravelTimeCalculatorMap *getTravelTimeCalculatorMap() const;
    [[nodiscard]] std::unique_ptr<ULocator::TravelTimeCalculatorMap> releaseTravelTimeCalculatorMap();
    /// @}

    /// @name Topography
    /// @{

    /// @brief Sets the topography calculator.
    /// @param[in,out] topography  The topography calculator to move to this.
    virtual void setTopography(std::unique_ptr<ULocator::Topography::ITopography> &&topography);
    /// @result The topography calculator.
    virtual const ULocator::Topography::ITopography *getTopography() const;
    /// @result The topography calculator.
    [[nodiscard]] virtual std::unique_ptr<ULocator::Topography::ITopography> releaseTopography();
    /// @result True indicates the topography was set.
    [[nodiscard]] virtual bool haveTopography() const noexcept;
    /// @}

    /// @name Arrivals 
    /// @{

    /// @param[in] arrivals  The seismic phase arrivals for this event.
    /// @throws std::runtime_error if \c haveTravelTimeCalculatorMap() is false.
    /// @throws std::invalid_argument if there are no arrivals.
    virtual void setArrivals(const std::vector<ULocator::Arrival> &arrivals);
    /// @result A constant reference to the arrivals.
    [[nodiscard]] virtual const std::vector<ULocator::Arrival> &getArrivalsReference() const noexcept;
    /// @}

    /// @name Origin
    /// @{

    /// @result The optimized-for origin.
    /// @throws std::runtime_error if \c haveOrigin() is false.
    [[nodiscard]] virtual ULocator::Origin getOrigin() const;
    /// @result True indicates the origin was set.
    [[nodiscard]] virtual bool haveOrigin() const noexcept;
    /// @}
 
    /// @brief Destructor.
    virtual ~IOptimizer();

    IOptimizer(const IOptimizer &) = delete;
    IOptimizer& operator=(const IOptimizer &) = delete; 
protected:
    /// @brief This utility sets the origin information.
    /// @param[in] origin   The event origin.
    virtual void setOrigin(const ULocator::Origin &origin);
    /// @brief This utility sets the origin information.
    /// @param[in] origin   The event origin.
    virtual void setOrigin(ULocator::Origin &&origin);
private:
    class IOptimizerImpl;
    std::unique_ptr<IOptimizerImpl> pImpl;
};
}
#endif
