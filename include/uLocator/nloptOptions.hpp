#ifndef ULOCATOR_NLOPT_OPTIONS_HPP
#define ULOCATOR_NLOPT_OPTIONS_HPP
#include <memory>
#include <vector>
#include <array>
namespace ULocator
{
class Topography;
class TravelTimeCalculatorMap;
class Quarry;
}
namespace ULocator
{
class NLOptOptions
{
public:
    enum class Optimizer
    {
        DIRECT,    /*!< DI(vidided)RECT angles.  This is a global optimization
                        utility that is useful for initially guessing an initial
                        solution.  This is most analogous to a */
        DIRECT_L,  /*!< DIRECT L.  This adds a slight bias to the DIRECT
                        refinement algorithm that seems to make it work better. */
        STOGO,     /*!< Stochastic Gradient Optimization.  This performs a 
                        global search like DIRECT but intelligently investigates
                        a cell with gradient information. */
        COBYLA,    /*!< Local search method.  This can handle non-linear
                        topographic constraints.  For objective functions with 
                        continuous second derivatives try BOBYQA. */
        BOBYQA     /*!< Local search method.  Unlike COBYLA this creates
                        a local quadratic approximation of the objective
                        function. */
    };
    enum class Region
    {
        Unknown = 0,
        Utah,
        Yellowstone 
    };
    enum class ObjectiveFunction
    {
        LeastSquares,
        L1,
        LP,
        DoubleDifferenceL1,
        DoubleDifferenceLeastSquares
    };
public:
    NLOptOptions();
    explicit NLOptOptions(Region region);
    NLOptOptions(const NLOptOptions &options);
    NLOptOptions(NLOptOptions &&options) noexcept;

    void setObjectiveFunction(ObjectiveFunction objectiveFunction) noexcept;
    [[nodiscard]] ObjectiveFunction getObjectiveFunction() const noexcept;

    void setUTMZone(int zone, bool northernHemisphere = true);
    [[nodiscard]] int getUTMZone() const noexcept;
    [[nodiscard]] bool isNorthernHemisphere() const noexcept;

    void setDefaultElevation(double elevation);
    [[nodiscard]] double getDefaultElevation() const noexcept;

    void setLatitudeBoundaries(const std::array<double, 2> &boundaries);
    [[nodiscard]] std::array<double, 2> getLatitudeBoundaries() const;
    [[nodiscard]] bool haveLatitudeBoundaries() const noexcept;

    void setLongitudeBoundaries(const std::array<double, 2> &boundaries);
    [[nodiscard]] std::array<double, 2> getLongitudeBoundaries() const;
    [[nodiscard]] bool haveLongitudeBoundaries() const noexcept;

    void setDepthBoundaries(const std::array<double, 2> &boundaries);
    [[nodiscard]] std::array<double, 2> getDepthBoundaries() const;
    [[nodiscard]] bool haveDepthBoundaries() const noexcept;
    

    /// @brief When the change in the model parameters is less than this amount
    ///        then the optimization will terminate.
    void setAbsoluteModelTolerance(double tolerance);
    /// @result The model stopping criteria.  By default this is 1.e-4 which
    ///         translates to about a few meters in latitude, longitude, and depth.
    [[nodiscard]] double getAbsoluteModelTolerance() const noexcept;
    /// @brief When the change in the model parameters is less than this amount
    ///        then the initial optimization will terminate.
    /// @throws std::invalid_argument if this is not positive.
    void setInitialAbsoluteModelTolerance(double tolerance);
    /// @result The initial model stopping criteria.  By default this is 1.e-4
    ///         which translates to about a few meters in latitude, longitude,
    ///         and depth.
    [[nodiscard]] double getInitialAbsoluteModelTolerance() const noexcept;

    /// @brief The initial search depth for earthquakes relative to sea-level in meters.
    void setInitialEarthquakeSearchDepth(double setDepth);
    /// @result The initial search depth in meters.
    [[nodiscard]] double getInitialEarthquakeSearchDepth() const noexcept;
    /// @brief The initial search depth relative to sea-level in meters.
    void setInitialQuarryBlastSearchDepth(double setDepth);
    /// @result The initial search depth in meters.  By default this is -1500.
    [[nodiscard]] double getInitialQuarryBlastSearchDepth() const noexcept;

    /// @result The refined latitude search will be within +/- this distance
    ///         specified in meters.
    [[nodiscard]] double getLatitudeRefinement() const noexcept;
    /// @result The refined longitude search will be within +/- this distance
    ///         specified in meters.
    [[nodiscard]] double getLongitudeRefinement() const noexcept;

    [[nodiscard]] std::vector<double> getRefinementSearchDepths() const noexcept;

    /// @brief The global search methods will stop after this many function
    ///        evaluations.
    void setMaximumNumberOfFunctionEvaluations(int nEvaluations);
    /// @result The maximum number of global search function evaluations.
    ///         By default this is 10000.
    [[nodiscard]] int getMaximumNumberOfFunctionEvaluations() const noexcept;
    /// @brief The initial global search methods will stop after this many
    ///        function evaluations.
    void setInitialMaximumNumberOfFunctionEvaluations(int nEvaluations);
    /// @result The maximum number of global search function evaluations.
    ///         By default this is 8000.
    [[nodiscard]] int getInitialMaximumNumberOfFunctionEvaluations() const noexcept;
    
    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~NLOptOptions();

    [[nodiscard]] std::vector<Quarry> getQuarries() const;

    NLOptOptions& operator=(const NLOptOptions &options);
    NLOptOptions& operator=(NLOptOptions &&options) noexcept;
private:
    class NLOptOptionsImpl;
    std::unique_ptr<NLOptOptionsImpl> pImpl;
};
}
#endif
