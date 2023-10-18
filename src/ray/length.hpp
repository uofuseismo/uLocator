#ifndef PRIVATE_RAY_LENGTH_HPP
#define PRIVATE_RAY_LENGTH_HPP
#include <cmath>
#include "eikonalxx/ray/point2d.hpp"
#include "eikonalxx/ray/point3d.hpp"
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
[[nodiscard]] [[maybe_unused]]
double computeLength(const EikonalXX::Ray::Point3D &a,
                     const EikonalXX::Ray::Point3D &b) 
{
    auto dx = b.getPositionInX() - a.getPositionInX();
    auto dy = b.getPositionInY() - a.getPositionInY();
    auto dz = b.getPositionInZ() - a.getPositionInZ();
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}
}
#endif
