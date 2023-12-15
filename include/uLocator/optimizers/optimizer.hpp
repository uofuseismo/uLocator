#ifndef ULOCATOR_OPTIMIZERS_OPTIMIZER_HPP
#define ULOCATOR_OPTIMIZERS_OPTIMIZER_HPP
#include <memory>
#include <umps/logging/log.hpp>
namespace ULocator
{
 class Arrival;
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
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    IOptimizer();
    /// @brief Constructor with a logger.
    explicit IOptimizer(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @}

    /// @brief Locates the event.
    virtual void locate() = 0;

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
 
    /// @brief Destructor.
    virtual ~IOptimizer();

    IOptimizer(const IOptimizer &) = delete;
    IOptimizer& operator=(const IOptimizer &) = delete; 
private:
    class IOptimizerImpl;
    std::unique_ptr<IOptimizerImpl> pImpl;
};
}
#endif
