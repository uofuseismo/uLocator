#include "uLocator/travelTimeCalculator.hpp"

using namespace ULocator;

double ITravelTimeCalculator::evaluate(
    const double t0, const double x, const double y, const double z,
    const bool applyCorrection) const
{
    return evaluate(t0, x, y, z,
                    nullptr, nullptr, nullptr, nullptr,
                    applyCorrection);
}

double ITravelTimeCalculator::evaluate(
    const double x, const double y, const double z,
    const bool applyCorrection) const
{
    constexpr double zero{0};
    return evaluate(zero, x, y, z,
                    nullptr, nullptr, nullptr, nullptr,
                    applyCorrection);
}

double ITravelTimeCalculator::evaluate(
    const double x, const double y, const double z,
    double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    constexpr double zero{0};
    return evaluate(zero, x, y, z,
                    nullptr, dtdx, dtdy, dtdz,
                    applyCorrection);
}

double ITravelTimeCalculator::operator()(
    const double x, const double y, const double z,
    const bool applyCorrection) const
{
    return evaluate(x, y, z, applyCorrection);
}

double ITravelTimeCalculator::operator()(
    const double t0, const double x, const double y, const double z,
    const bool applyCorrection) const
{
    return evaluate(t0, x, y, z, applyCorrection);
}
