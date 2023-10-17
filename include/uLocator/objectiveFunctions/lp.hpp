#ifndef ULOCATOR_OBJECTIVE_FUNCTIONS_LP_HPP
#define ULOCATOR_OBJECTIVE_FUNCTIONS_LP_HPP
#include <memory>
#include "uLocator/objectiveFunctions/objectiveFunction.hpp"
namespace ULocator::ObjectiveFunctions
{
class LP : public IObjectiveFunction
{
public:
    LP(double p = 1.5);
    explicit LP(std::shared_ptr<UMPS::Logging::ILog> &logger, double p = 1.5);
    //void initialize(); 
    void setInversionStrategy(InversionStrategy strategy) noexcept override;
    void setArrivals(const std::vector<Arrival> &arrivals) override;
    double operator()(unsigned int, const double *, double *) override;
    ~LP() override;
private:
    class LPImpl;
    std::unique_ptr<LPImpl> pImpl;
};
}
#endif
