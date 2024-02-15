#ifndef ULOCATOR_TRAVEL_TIME_CALCULATOR_MAP_HPP
#define ULOCATOR_TRAVEL_TIME_CALCULATOR_MAP_HPP
#include <vector>
#include <string>
#include <memory>
namespace ULocator
{
class Station;
class ITravelTimeCalculator;
}
namespace ULocator
{
/// @class TravelTimeCalculatorMap "travelTimeCalculatorMap.hpp" "uLocator/travelTimeCalculatorMap.hpp"
/// @brief A travel time calculator map is a collection of travel time
///        calculators that takes a (station, phase) pair and  evaluates
///        the corresponding travel time for a source at a given position.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class TravelTimeCalculatorMap
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    TravelTimeCalculatorMap();
    /// @brief Move constructor.
    /// @param[in,out] calculators  The calculator map from which to initialize
    ///                             this class.  On exit, the behavior is of
    ///                             calculators is undefined.
    TravelTimeCalculatorMap(TravelTimeCalculatorMap &&calculators) noexcept;
    /// @}

    /// @brief Adds a travel time calculator.
    /// @param[in] station     The station information (primarily the network
    ///                        code and station name).
    /// @param[in] phase       The phase - e.g., P or S.
    /// @param[in] calculator  The travel time calculator for this
    ///                        station/phase pair.  Note, calculator's 
    ///                        behavior will be undefined after this class
    ///                        takes ownership.
    void insert(const Station &station,
                const std::string &phase,
                std::unique_ptr<const ITravelTimeCalculator> &&calculator);
    /// @param[in] station   The station information (network code and name).
    /// @param[in] phase     The phase.
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
                                  double originTime, double xSource, double ySource, double zSource,
                                  bool applyCorrection = true) const;
    [[nodiscard]] double evaluate(const Station &station, const std::string &phase,
                                  double originTime, double xSource, double ySource, double zSource,
                                  double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
                                  bool applyCorrection = true) const;
    void evaluate(const std::vector<std::pair<Station, std::string>> &stationPhase,
                  double originTime, double xSource, double ySource, double zSource,
                  std::vector<double> *travelTimes,
                  bool applyCorrection = true) const;
    void evaluate(const std::vector<std::pair<Station, std::string>> &stationPhase,
                  double originTime, double xSource, double ySource, double zSource,
                  std::vector<double> *travelTimes,
                  std::vector<double> *dtdt0,
                  std::vector<double> *dtdx,
                  std::vector<double> *dtdy,
                  std::vector<double> *dtdz,
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
