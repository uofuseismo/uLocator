#ifndef ULOCATOR_EIKONALXX_RAY_SEGMENT_2D_HPP
#define ULOCATOR_EIKONALXX_RAY_SEGMENT_2D_HPP
#include <ostream>
#include <memory>
namespace EikonalXX::Ray
{
class Point2D;
/// @class Segment2D "segment2d.hpp" "eikonalxx/ray/segment2d.hpp"
/// @brief Defines a ray segment in a 2D model.  A segment is a line between
///        a start and end point.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Segment2D
{
public:
    /// @name Constructors
    /// @{

    /// @brief Default constructor.
    Segment2D();
    /// @brief Copy constructor.
    /// @param[in] segment  The segment class from which to initialize this
    ///                     class.
    Segment2D(const Segment2D &segment);
    /// @brief Move constructor.
    /// @param[in,out] segment  The segment class from which to initialize this
    ///                         class. On exit, segment's behavior is undefined.
    Segment2D(Segment2D &&segment) noexcept;
    /// @} 

    /// @name Geometry
    /// @{

    /// @brief Defines the start and end point of the ray segment.
    /// @param[in] startAndEndPoint  startAndEndPoint.first is the ray segment's
    ///                              start point and startAndEndPoint.second is
    ///                              the ray segment's end point.
    /// @throws std::invalid_argument if start or end point do not have 
    ///         x and z positions.
    void setStartAndEndPoint(const std::pair<Point2D, Point2D> &startAndEndPoint);
    /// @result The start point of the ray. 
    /// @throws std::runtime_error if \c haveStartAndEndPoint() is false.
    [[nodiscard]] Point2D getStartPoint() const;
    /// @result The end point of the ray.
    /// @throws std::runtime_error if \c haveStartAndEndPoint() is false.
    [[nodiscard]] Point2D getEndPoint() const;
    /// @result True indicates the ray start and end point were set.
    [[nodiscard]] bool haveStartAndEndPoint() const noexcept;
    /// @brief Interchanges the start and end point.
    /// @throws std::runtime_error if \c haveStartAndEndPoint() is false.
    void reverse();
    /// @result The ray segment's length in meters.
    /// @throws std::runtime_error if \c haveStartAndEndPoint() is false.
    [[nodiscard]] double getLength() const;
    /// @}

    /// @name Velocity
    /// @{

    /// @brief Sets the velocity along the ray segment.
    /// @param[in] velocity  The velocity in m/s.
    /// @throws std::invalid_argument if this is not positive.
    void setVelocity(double velocity); 
    /// @brief Sets the slowness along the ray segment.
    /// @param[in] slowness  The slowness in s/m.
    /// @throws std::invalid_argument if this is not positive.
    void setSlowness(double slowness);
    /// @result The velocity (m/s) along the ray.
    /// @throws std::runtime_error if \c haveVelocity() is false.
    [[nodiscard]] double getVelocity() const;
    /// @result The slowness (s/m) along the ray.
    /// @throws std::runtime_error if \c haveVelocity() is false.
    [[nodiscard]] double getSlowness() const;
    /// @result True indicates the velocity was set.
    [[nodiscard]] bool haveVelocity() const noexcept;
    /// @}

    /// @name Travel Time
    /// @{

    /// @result The travel time along the ray segment in seconds.
    /// @throws std::runtime_error if \c haveVelocity() or
    ///         \c haveStartAndEndPoint() is false.
    [[nodiscard]] double getTravelTime() const;
    /// @}

    /// @name Cell Index
    /// @{

    /// @brief In inversion problems this is the cell through which the ray
    ///        passes.  This is useful when building that large, sparse
    ///        modeling matrix.
    void setVelocityModelCellIndex(int index);
    /// @result The velocity model cell index.
    [[nodiscard]] int getVelocityModelCellIndex() const;
    /// @result True indicates the velocity model cell index was set.
    [[nodiscard]] bool haveVelocityModelCellIndex() const noexcept;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] segment  The segment to copy to this.
    /// @result A deep copy of the input segment.
    Segment2D& operator=(const Segment2D &segment);
    /// @brief Move assignment operator.
    /// @param[in,out] segment  The segment whose memory will be moved to this.
    ///                         On exit, segment's behavior is undefined.
    /// @result The memory from segment moved to this.
    Segment2D& operator=(Segment2D &&segment) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Rests the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~Segment2D();
    /// @}
private:
    friend void swap(Segment2D &lhs, Segment2D &rhs) noexcept;
private:
    class Segment2DImpl;
    std::unique_ptr<Segment2DImpl> pImpl;
};
std::ostream& operator<<(std::ostream &os, const Segment2D &segment);
}
#endif
