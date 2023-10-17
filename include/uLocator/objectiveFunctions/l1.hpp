#ifndef ULOCATOR_OBJECTIVE_FUNCTIONS_L1_HPP
#define ULOCATOR_OBJECTIVE_FUNCTIONS_L1_HPP
#include <memory>
#include "uLocator/objectiveFunctions/objectiveFunction.hpp"
namespace ULocator::ObjectiveFunctions
{
class L1 : public IObjectiveFunction
{
public:
    L1();
    explicit L1(std::shared_ptr<UMPS::Logging::ILog> &logger);
    //void initialize(); 
    void setInversionStrategy(InversionStrategy strategy) noexcept override;
    void setArrivals(const std::vector<Arrival> &arrivals) override;
    double operator()(unsigned int, const double *, double *) override;
    ~L1() override;
private:
    class L1Impl;
    std::unique_ptr<L1Impl> pImpl;
};
}
#endif
