#ifndef ULOCATOR_OBJECTIVE_FUNCTIONS_LEAST_SQUARES_HPP
#define ULOCATOR_OBJECTIVE_FUNCTIONS_LEAST_SQUARES_HPP
#include <memory>
#include "uLocator/objectiveFunctions/objectiveFunction.hpp"
namespace ULocator::ObjectiveFunctions
{
class LeastSquares : public IObjectiveFunction
{
public:
    LeastSquares();
    explicit LeastSquares(std::shared_ptr<UMPS::Logging::ILog> &logger);
    //void initialize(); 
    void setInversionStrategy(InversionStrategy strategy) noexcept override;
    void setArrivals(const std::vector<Arrival> &arrivals) override;
    double operator()(unsigned int, const double *, double *) override;
    //double f(int n, const double x[]) const override;
    //double g(int n, const double x[], double g[]) const override;
    /// Destructor
    ~LeastSquares() override;
private:
    class LeastSquaresImpl;
    std::unique_ptr<LeastSquaresImpl> pImpl;
};
}
#endif
