#ifndef ULOCATOR_NLOPT_HPP
#define ULOCATOR_NLOPT_HPP
#include <umps/logging/log.hpp>
#include <vector>
#include <memory>
namespace ULocator
{
 class Arrival;
 class Origin;
 class NLOptOptions;
 class Topography;
 class TravelTimeCalculatorMap;
}
namespace ULocator
{
class NLOpt
{
public:
    enum class SourceDepthConstraint
    {
        Free,
        FixedToFreeSurface,
        Fixed
    };
    enum class InitialGuessStrategy
    {
        InternalHeuristic,
        SpecifyHypocenter 
    };
public:
    NLOpt();
    explicit NLOpt(std::shared_ptr<UMPS::Logging::ILog> &logger);

    void setOptions(const NLOptOptions &options);

    void setTravelTimeCalculatorMap(std::unique_ptr<const TravelTimeCalculatorMap> &&calculators);
    [[nodiscard]] std::unique_ptr<const TravelTimeCalculatorMap> releaseTravelTimeCalculatorMap();

    void setTopography(std::unique_ptr<const Topography> &&topgraphy);
    [[nodiscard]] std::unique_ptr<const Topography> releaseTopography();


//    void initialize(const NLOptOptions &options);
    [[nodiscard]] bool isInitialized() const noexcept;



    void setArrivals(const std::vector<Arrival> &arrivals);
    void clearArrivals() noexcept;
    [[nodiscard]] bool haveArrivals() const noexcept;

    void locateQuarryBlast();
    void locateEarthquake();
    void locateQuarryBlast(const Origin &initialGuess);
    void locateEarthquake(const Origin &initialGuess);
    [[nodiscard]] double evaluateLoss(const ULocator::Origin &origin) const;
    [[nodiscard]] Origin predict(const ULocator::Origin &origin, const bool applyCorrection = true);
    //void locate(const SourceDepthConstraint sourceConstraint = SourceDepthConstraint::Free,
    //            double sourceDepth = 6400);
    [[nodiscard]] Origin getOrigin() const;
    [[nodiscard]] bool haveOrigin() const noexcept;

    [[nodiscard]] int getNumberOfObjectiveFunctionEvaluations() const noexcept;
    [[nodiscard]] int getNumberOfGradientEvaluations() const noexcept;

    void clear() noexcept;
    ~NLOpt();

    
    NLOpt(const NLOpt &) = delete;
    NLOpt& operator=(const NLOpt &) = delete;
private:
    class NLOptImpl;
    std::unique_ptr<NLOptImpl> pImpl;
};
}
#endif
