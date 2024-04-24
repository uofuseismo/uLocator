#ifndef ULOCATOR_EIKONALXX_RAY_POINT_2D_HPP
#define ULOCATOR_EIKONALXX_RAY_POINT_2D_HPP
#include <ostream>
#include <memory>
namespace EikonalXX::Ray
{
/// @class Point2D "point2d.hpp" "eikonalxx/ray/point2d.hpp"
/// @brief Defines a point along a ray path.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Point2D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Point2D();
    /// @brief Copy constructor.
    /// @param[in] point  The point from which to initialize this class.
    Point2D(const Point2D &point);
    /// @brief Move constructor.
    /// @param[in,out] point  The point from which to initialize this class. 
    ///                       On exit, point's behavior is undefined.
    Point2D(Point2D &&point) noexcept;
    /// @brief Constructs a point from an (x, z) position.
    /// @param[in] x  The position in the model in x.  This has units of meters.
    /// @param[in] z  The position in the model in z.  This has units of meters.
    Point2D(double x, double z);
    /// @}

    /// @name X Position
    /// @{

    /// @param[in] x  The position in the model in x.  This has units of meters.
    void setPositionInX(double x) noexcept;
    /// @result The position in the model in x (meters).
    [[nodiscard]] double getPositionInX() const;
    /// @result True indicates the x position was set.
    [[nodiscard]] bool havePositionInX() const noexcept;
    /// @}

    /// @name Z Position
    /// @{

    /// @param[in] z  The position in the model in z.  This has units of meters.
    void setPositionInZ(double z) noexcept;
    /// @result The position in the model in z (meters).
    [[nodiscard]] double getPositionInZ() const;
    /// @result True indicates the z position was set.
    [[nodiscard]] bool havePositionInZ() const noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] point  The point to copy to this.
    /// @result A deep copy of the input point.
    Point2D& operator=(const Point2D &point);
    /// @brief Move assignment operator.
    /// @param[in,out] point  The point whose memory will be moved to this.
    ///                       On exit, point's behavior is undefined.
    /// @result The memory from point moved to this.
    Point2D& operator=(Point2D &&point) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Point2D();
    /// @}
private:
    friend void swap(Point2D &lhs, Point2D &rhs) noexcept;
private:
    class Point2DImpl;
    std::unique_ptr<Point2DImpl> pImpl;
};
void swap(Point2D &lhs, Point2D &rhs) noexcept;
std::ostream& operator<<(std::ostream &os, const Point2D &point);
}
#endif 
