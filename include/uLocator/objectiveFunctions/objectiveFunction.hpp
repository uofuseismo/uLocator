#ifndef ULOCATOR_OBJECTIVE_FUNCTIONS_OBJECTIVE_FUNCTION_HPP
#define ULOCATOR_OBJECTIVE_FUNCTIONS_OBJECTIVE_FUNCTION_HPP
#include <memory>
#include <map>
#include <string>
#include <vector>
#include <functional>
#include <Eigen/Dense>
namespace UMPS::Logging
{
 class ILog;
}
namespace ULocator
{
 class Arrival;
 class Origin;
 class TravelTimeCalculatorMap;
 class Topography;
}
namespace ULocator
{
enum class InversionStrategy 
{
    LatitudeLongitudeDepth,
    LatitudeLongitudeFixedDepth,
    LatitudeLongitudeDepthAtFreeSurface
};
}
namespace ULocator::ObjectiveFunctions
{
class IObjectiveFunction
{
public:
    /// @brief Constructor.
    IObjectiveFunction(); 
    explicit IObjectiveFunction(std::shared_ptr<UMPS::Logging::ILog> &logger);
    /// @brief Destructor.
    virtual ~IObjectiveFunction();
    /// @brief Sets the inversion strategy.
    virtual void setInversionStrategy(const InversionStrategy strategy);
    [[nodiscard]] virtual bool haveInversionStrategy() const noexcept;
    virtual void setUTMZone(int utmZone, bool north = true);
    [[nodiscard]] virtual int getUTMZone() const noexcept;
    virtual void toggleApplyCorrection(bool applyCorrection) noexcept;
    [[nodiscard]] bool applyCorrection() const noexcept;

    void setTravelTimeCalculatorMap(std::unique_ptr<const ULocator::TravelTimeCalculatorMap> &&travelTimeCalculatorMap);
    [[nodiscard]] bool haveTravelTimeCalculatorMap() const noexcept;
    [[nodiscard]] std::unique_ptr<const ULocator::TravelTimeCalculatorMap> releaseTravelTimeCalculatorMap();

    void setFixedDepthCallbackFunction(const std::function<double ()> &callback);
    [[nodiscard]] bool haveFixedDepthCallbackFunction() const noexcept;

    void setTopographyCallbackFunction(
         const std::function<double(const double latitude, const double longitude)> &callback);
    [[nodiscard]] bool haveTopographyCallbackFunction() const noexcept;

    void setLikelihoodCallbackFunction(
         const std::function<double(const Eigen::VectorXd &estimates)> &callback);
    [[nodiscard]] bool haveLikelihoodCallbackFunction() const noexcept;

    void setComputeSourceTimeCallbackFunction(
        const std::function<double(const Eigen::VectorXd &estimates)> &callback);
    [[nodiscard]] bool haveComputeSourceTimeCallbackFunction() const noexcept;

    void setGradientCallbackFunction(
        const std::function<void (const Eigen::MatrixXd &jacobian,
                                  const Eigen::VectorXd &estimates,
                                  Eigen::VectorXd &gradient)> &callback);
    [[nodiscard]] bool haveGradientCallbackFunction() const noexcept;

    void setJacobianCallbackFunction(
        const std::function<void (const std::map<std::string, std::vector<double>> &sensitivities,
                                  Eigen::MatrixXd &jacobian)> &callback);
    [[nodiscard]] bool haveJacobianCallbackFunction() const noexcept;
//    [[nodiscard]] std::unique_ptr<ULocator::Topo

    void toggleComputeSourceTimeInFAndG(bool computeSourceTime) noexcept;
    void toggleUseFiniteDifference(
         const bool useFiniteDifference, const double perturbation = 0.5);

  
    virtual void setArrivals(const std::vector<Arrival> &arrivals);
    [[nodiscard]] virtual const std::vector<Arrival>& getArrivalsReference() const final;
    [[nodiscard]] virtual std::vector<Arrival> getArrivals() const;
    [[nodiscard]] virtual int getNumberOfArrivals() const; 
    virtual void clearArrivals() noexcept;

    [[nodiscard]] virtual double evaluateLoss(const Origin &origin) const;
    [[nodiscard]] virtual ULocator::Origin predict(const Origin &origin) const;

    void resetCounters() noexcept;
    [[nodiscard]] virtual int getNumberOfObjectiveFunctionEvaluations() const noexcept;
    [[nodiscard]] virtual int getNumberOfGradientEvaluations() const noexcept;

    //virtual void setObservations(const Eigen::VectorXd &observations);

    /// @brief Sets the observations.  Here the weights will be unity.
    //virtual void setObservations(const std::vector<double> &observations);
    /// @brief Sets the observerations and weights.
    //virtual void setObservations(const std::vector<double> &observations,
    //                             const std::vector<double> &weights);
    /// @brief Sets the current model.
    //void setModel(const std::vector<double> &x);
    //virtual void setModel(const int n, const double x[]) = 0;
    /// @result The sum of the weights.
    //[[nodiscard]] virtual double getSumOfWeights() const;
    /// @result True indicates the observations were set.
    //[[nodiscard]] virtual bool haveObservations() const noexcept;
    /// @result The objective function evaluated at x. 
    [[nodiscard]] virtual double f(const int n, const double x[]) const;// = 0;
    /// @result The objective function and the gradient evaluated at x.
    [[nodiscard]] virtual double g(const int n, const double x[], double g[]) const;// = 0;
    /// @result The number of model parameters in this optimization.
    [[nodiscard]] virtual int getNumberOfModelParameters() const;
    /// @param[in] n          The number of model parameters.  This should match
    ///                       \c getNumberOfModelParameters().
    /// @param[in] x          The model - e.g., (utmx, utmy, depth).
    /// @param[out] gradient  If not nullptr then this is the gradient.  In this
    ///                       case this is an array of length [n].
    /// @result The objective function evaluated at the given x. 
    virtual double operator()(unsigned int n, const double *x, double *gradient) = 0;
    //[[nodiscard]] virtual double getLastComputedOriginTime() const noexcept;
private:
    class IObjectiveFunctionImpl;
    std::unique_ptr<IObjectiveFunctionImpl> pImpl;
};
}
#endif
