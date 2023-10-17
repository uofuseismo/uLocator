#ifndef ULOCATOR_OBJECTIVE_FUNCTIONS_DOUBLE_DIFFERENCE_L1_HPP
#define ULOCATOR_OBJECTIVE_FUNCTIONS_DOUBLE_DIFFERENCE_L1_HPP
#include <memory>
#include "uLocator/objectiveFunctions/objectiveFunction.hpp"
namespace ULocator::ObjectiveFunctions
{
class DoubleDifferenceL1 : public IObjectiveFunction
{
public:
    DoubleDifferenceL1();
    explicit DoubleDifferenceL1(std::shared_ptr<UMPS::Logging::ILog> &logger);
    //void initialize(); 
    void setInversionStrategy(InversionStrategy strategy) override;
    void setArrivals(const std::vector<Arrival> &arrivals) override;
    double operator()(unsigned int, const double *, double *) override;
    ~DoubleDifferenceL1() override;
private:
    class DoubleDifferenceL1Impl;
    std::unique_ptr<DoubleDifferenceL1Impl> pImpl;
};
}
#endif
