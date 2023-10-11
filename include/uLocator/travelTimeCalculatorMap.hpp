#ifndef ULOCATOR_TRAVEL_TIME_CALCULATOR_MAP_HPP
#define ULOCATOR_TRAVEL_TIME_CALCULATOR_MAP_HPP
#include <vector>
#include <string>
#include <memory>
namespace ULocator
{
class Station;
class ITravelTimeCalculator;
namespace Position
{
 class WGS84;
}
}
namespace ULocator
{
class TravelTimeCalculatorMap
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    TravelTimeCalculatorMap();
    /// @brief Move constructor.
    TravelTimeCalculatorMap(TravelTimeCalculatorMap &&calculators) noexcept;
    /// @}

    /// @brief Adds a travel time calculator.
    void insert(const Station &station,
                const std::string &phase,
                std::unique_ptr<const ITravelTimeCalculator> &&calculator);
    /// @param[in] stationPhase  The station/phase name pair.
    /// @result True indicates the map contains the given station, phase pair.
    [[nodiscard]] bool contains(const Station &station, const std::string &phase) const;

    /// @result The number of travel time calculators.
    [[nodiscard]] int size() const noexcept; 

    /// @name Operators
    /// @{

    /// @brief Move assignment operator.
    TravelTimeCalculatorMap& operator=(TravelTimeCalculatorMap &&calculators) noexcept;
    /// @result A pointer to the travel time calculator for this
    ///         station/phase pair. 
    const ULocator::ITravelTimeCalculator *at(const Station &station, const std::string &stationPhase) const;
    /// @}

    [[nodiscard]] double evaluate(const Station &station, const std::string &phase,
                                  const Position::WGS84 &epicenter, double depth,
                                  bool applyCorrection = true) const;
    [[nodiscard]] double evaluate(const Station &station, const std::string &phase,
                                  const Position::WGS84 &epicenter, double depth,
                                  double *dtdx, double *dtdy, double *dtdz,
                                  bool applyCorrection = true) const;

    void evaluate(const std::vector<std::pair<Station, std::string>> &stationPhase,
                  const Position::WGS84 &epicenter, double depth,
                  std::vector<double> *travelTimes,
                  bool applyCorrection = true) const;
    void evaluate(const std::vector<std::pair<Station, std::string>> &stationPhase,
                  const Position::WGS84 &epicenter, double depth,
                  std::vector<double> *travelTimes,
                  std::vector<double> *dtdx, std::vector<double> *dtdy, std::vector<double> *dtdz,
                  bool applyCorrection = true) const;

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~TravelTimeCalculatorMap();
    /// @}

    
    TravelTimeCalculatorMap(const TravelTimeCalculatorMap &) = delete;
    TravelTimeCalculatorMap& operator=(const TravelTimeCalculatorMap &) = delete;
private:
    class TravelTimeCalculatorMapImpl;
    std::unique_ptr<TravelTimeCalculatorMapImpl> pImpl;
};
}
#endif
