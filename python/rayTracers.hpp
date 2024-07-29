#ifndef ULOCATOR_RAY_TRACERS_HPP
#define ULOCATOR_RAY_TRACERS_HPP
#include <string>
#include <memory>
#include <uLocator/travelTimeCalculator.hpp>
#include <uLocator/uussRayTracer.hpp>
namespace ULocator::Python
{
 class Station;
}

namespace ULocator::Python
{
class ITravelTimeCalculator
{
public:
/*
    [[nodiscard]] virtual double evaluate(double t0, double x, double y, double z,
                                          double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
                                          bool applyCorrection) const = 0;
    [[nodiscard]] virtual double evaluate(double t0, double x, double y, double z,
                                          bool applyCorrection) const
    {
        return evaluate(t0, x, y, z,
                        nullptr, nullptr, nullptr, nullptr,
                        applyCorrection); 
    }
*/
    virtual ~ITravelTimeCalculator() = default;
    //virtual std::unique_ptr<ULocator::ITravelTimeCalculator> getBaseClass() const = 0;
};

namespace RayTracers
{

class Utah : public ULocator::Python::ITravelTimeCalculator
{
public:
    Utah(const ULocator::Python::Station &station,
         const ULocator::UUSSRayTracer::Phase phase);
    ~Utah() override;
    [[nodiscard]] double
       computeArrivalTime(double t0, double x, double y, double z,
                         bool applyCorrection = true) const;
    [[nodiscard]] std::tuple<double, double, double, double, double>
       computeArrivalTimeAndDerivatives(double t0, double x, double y, double z,
                                        bool applyCorrection = true) const;
    //[[nodiscard]] std::unique_ptr<ULocator::ITravelTimeCalculator> getBaseClass() const;
/*
    [[nodiscard]] virtual double operator()(double x, double y,
                                            double *dElevationDx, double *dElevationDy) const
    {
        return evaluate(x, y, dElevationDx, dElevationDy);
    }
*/
private:
    std::unique_ptr<ULocator::UUSSRayTracer> pImpl{nullptr};
};

class YNP : public ULocator::Python::ITravelTimeCalculator
{
public:
    YNP(const ULocator::Python::Station &station,
         const ULocator::UUSSRayTracer::Phase phase);
    ~YNP() override;
    [[nodiscard]] double
       computeArrivalTime(double t0, double x, double y, double z,
                          bool applyCorrection = true) const;
    [[nodiscard]] std::tuple<double, double, double, double, double>
       computeArrivalTimeAndDerivatives(double t0, double x, double y, double z,
                                        bool applyCorrection = true) const;
    //[[nodiscard]] std::unique_ptr<ULocator::ITravelTimeCalculator> getBaseClass() const;
/*
    [[nodiscard]] virtual double operator()(double x, double y,
                                            double *dElevationDx, double *dElevationDy) const
    {
        return evaluate(x, y, dElevationDx, dElevationDy);
    }
*/
private:
    std::unique_ptr<ULocator::UUSSRayTracer> pImpl{nullptr};
};


void initialize(pybind11::module &m);
}

}
#endif

