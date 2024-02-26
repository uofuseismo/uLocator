#ifndef PRIVATE_RAY_LENGTH_HPP
#define PRIVATE_RAY_LENGTH_HPP
#include <cmath>
#include "point2d.hpp"
namespace
{
[[nodiscard]] [[maybe_unused]]
double computeLength(const EikonalXX::Ray::Point2D &a,
                     const EikonalXX::Ray::Point2D &b) 
{
    auto dx = b.getPositionInX() - a.getPositionInX();
    auto dz = b.getPositionInZ() - a.getPositionInZ();
    return std::hypot(dx, dz);
}
}
#endif
