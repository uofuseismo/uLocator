#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <uLocator/uussRayTracer.hpp>
#include <uLocator/station.hpp>
#include "dataStructures.hpp"
#include "position.hpp"
#include "rayTracers.hpp"

using namespace ULocator::Python::RayTracers;

Utah::Utah(const ULocator::Python::Station &station,
           const ULocator::UUSSRayTracer::Phase phase) :
    pImpl(std::make_unique<ULocator::UUSSRayTracer> (station.getBaseClassReference(),
                                                     phase,
                                                     ULocator::UUSSRayTracer::Region::Utah,
                                                     nullptr))
{
}

Utah::~Utah() = default;

double Utah::computeArrivalTime(const double t0,
                                const double x, const double y, const double z,
                                bool applyCorrection) const
{
    if (pImpl)
    {
        return pImpl->evaluate(t0, x, y, z,
                               nullptr, nullptr, nullptr, nullptr,
                               applyCorrection);
    }
    else
    {
        throw std::runtime_error("Travel time calculator not initialized");
    }
}

std::tuple<double, double, double, double, double> 
Utah::computeArrivalTimeAndDerivatives(
    const double t0,
    const double x, const double y, const double z,
    const bool applyCorrection) const
{
    if (pImpl)
    { 
        double dtdt0, dtdx, dtdy, dtdz;
        auto travelTime = pImpl->evaluate(t0, x, y, z,
                                          &dtdt0, &dtdx, &dtdy, &dtdz,
                                          applyCorrection);
        return std::tuple {travelTime, dtdt0, dtdx, dtdy, dtdz};
    }
    else
    {   
        throw std::runtime_error("Travel time calculator not initialized");
    }
}

//----------------------------------------------------------------------------//

YNP::YNP(const ULocator::Python::Station &station,
         const ULocator::UUSSRayTracer::Phase phase) :
    pImpl(std::make_unique<ULocator::UUSSRayTracer> (station.getBaseClassReference(),
                                                     phase,
                                                     ULocator::UUSSRayTracer::Region::YNP,
                                                     nullptr))
{
}

YNP::~YNP() = default;

double YNP::computeArrivalTime(const double t0, 
                               const double x, const double y, const double z,
                               bool applyCorrection) const
{
    if (pImpl)
    {   
        return pImpl->evaluate(t0, x, y, z,
                               nullptr, nullptr, nullptr, nullptr,
                               applyCorrection);
    }   
    else
    {   
        throw std::runtime_error("Travel time calculator not initialized");
    }   
}

std::tuple<double, double, double, double, double>
YNP::computeArrivalTimeAndDerivatives(
    const double t0,
    const double x, const double y, const double z,
    const bool applyCorrection) const
{
    if (pImpl)
    {
        double dtdt0, dtdx, dtdy, dtdz;
        auto travelTime = pImpl->evaluate(t0, x, y, z,
                                          &dtdt0, &dtdx, &dtdy, &dtdz,
                                          applyCorrection);
        return std::tuple {travelTime, dtdt0, dtdx, dtdy, dtdz};
    }
    else
    {
        throw std::runtime_error("Travel time calculator not initialized");
    }
}


//----------------------------------------------------------------------------//

void ULocator::Python::RayTracers::initialize(pybind11::module &m) 
{
    pybind11::module rayTracerModule = m.def_submodule("RayTracers");
    rayTracerModule.attr("__doc__") = "Defines layer tracers for computing travel times";

    pybind11::class_<ULocator::Python::ITravelTimeCalculator> (rayTracerModule, "TravelTimeCalculator");


    // Enums
    pybind11::enum_<ULocator::UUSSRayTracer::Region> (rayTracerModule, "Region")
        .value("Utah", ULocator::UUSSRayTracer::Region::Utah,
               "Utah region.")
        .value("YNP",  ULocator::UUSSRayTracer::Region::YNP,
               "Yellowstone national park.");
    pybind11::enum_<ULocator::UUSSRayTracer::Phase> (rayTracerModule, "Phase")
        .value("P", ULocator::UUSSRayTracer::Phase::P,
               "Predict a first-arriving P phase.")
        .value("S", ULocator::UUSSRayTracer::Phase::S,
               "Predict a first-arriving S phase.");

    ///----------------------------------Utah--------------------------------///
    pybind11::class_<ULocator::Python::RayTracers::Utah,
                     ULocator::Python::ITravelTimeCalculator> utah(rayTracerModule, "Utah");
    utah.def(pybind11::init<const ULocator::Python::Station,
                            const ULocator::UUSSRayTracer::Phase> ());
    utah.def("compute_arrival_time",
             &Utah::computeArrivalTime,
             "Given the origin time, local x and y epicenter, and depth in meters this computes the arrival time at the station in seconds.");
    utah.def("compute_arrival_time_and_derivatives",
             &Utah::computeArrivalTimeAndDerivatives,
             "Given the origin time, local x and y epicenter, and depth in meters this computes the (1) arrival time at the station in second (2) derivative of the travel-time, (3) the derivative of the x position with respect to time, (4) the derivative of the y position with respect to time, (5) and the derivative of the depth respect to time.");
 
    ///------------------------------------YNP-------------------------------///
    pybind11::class_<ULocator::Python::RayTracers::YNP,
                     ULocator::Python::ITravelTimeCalculator> ynp(rayTracerModule, "YNP");
    ynp.def(pybind11::init<const ULocator::Python::Station,
                           const ULocator::UUSSRayTracer::Phase> ());
    ynp.def("compute_arrival_time",
            &YNP::computeArrivalTime,
            "Given the origin time, local x and y epicenter, and depth in meters this computes the arrival time at the station in seconds.");
    ynp.def("compute_arrival_time_and_derivatives",
            &YNP::computeArrivalTimeAndDerivatives,
            "Given the origin time, local x and y epicenter, and depth in meters this computes the (1) arrival time at the station in second (2) derivative of the travel-time, (3) the derivative of the x position with respect to time, (4) the derivative of the y position with respect to time, (5) and the derivative of the depth respect to time.");



}

