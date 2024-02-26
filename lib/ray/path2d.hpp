#ifndef ULOCATOR_EIKONALXX_RAY_PATH_2D_HPP
#define ULOCATOR_EIKONALXX_RAY_PATH_2D_HPP
#include <memory>
#include <vector>
namespace EikonalXX::Ray
{
 class Segment2D;
}
namespace EikonalXX::Ray
{
/// @class Path2D "path2d.hpp" "eikonalxx/ray/path2d.hpp"
/// @brief Defines a ray path in two dimensions.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Path2D
{
private:
    using PathType = std::vector<Segment2D>;
public:
    using iterator = typename PathType::iterator;
    using const_iterator = typename PathType::const_iterator;
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Path2D();
    /// @brief Copy constructor.
    /// @param[in] path  The path from which to initialize this class.
    Path2D(const Path2D &path);
    /// @brief Move constructor.
    /// @param[in,out] path  The path from which to initialize this class.
    ///                      On exit, path's behavior is undefined.
    Path2D(Path2D &&path) noexcept;
    /// @}

    /// @name Option 1 For Path Creation
    /// @{

    /// @brief Opens the path construction process.
    void open();
    /// @result True indicates the path is open for construction.
    [[nodiscard]] bool isOpen() const noexcept;
    /// @brief Appends to the path.
    /// @param[in] segment  The segment to append to the path.
    /// @throws std::invalid_argument if the segment does not have a start/end
    ///         point and a velocity.  Additionally, if this is the not the
    ///         source segment then the segment's start point must coincide
    ///         with the previous segment's end point.
    /// @throws std::runtime_error if \c isOpen() is false.
    void append(const Segment2D &segment);
    /// @brief Appends to the path.
    /// @param[in,out] segment  The segment to append to the path.
    ///                         On exit, segment's behavior is undefined.
    void append(Segment2D &&segment);
    /// @brief Closes the path construction process.
    /// @throws std::runtime_error if \c isOpen() is false.
    void close();
    /// @}

    /// @name Option 2 For Path Creation
    /// @{

    /// @brief Sets the ray path.
    /// @param[in] segments  The segments to set.
    /// @throws std::invalid_argument if the i'th segment's end point does
    ///         not equal the i+1'th segment's start point, or the velocity
    ///         is not set on a segment.
    void set(const std::vector<Segment2D> &segments);
    /// @brief Sets the ray path.
    /// @param[in,out] segments  The segments to set.  On exit, segments
    ///                          behavior is undefined.
    /// @throws std::invalid_argument if the i'th segment's end point does
    ///         not equal the i+1'th segment's start point, or the velocity
    ///         is not set on a segment.
    void set(std::vector<Segment2D> &&segments);
    /// @brief Reverses the ray path.
    /// @throws std::runtime_error if the ray path construction is open.
    void reverse();
    /// @}

    /// @result True indicates that there are no segments.
    [[nodiscard]] bool empty() const noexcept; 
    /// @result The number of segments. 
    [[nodiscard]] size_t size() const noexcept;
    /// @result The total travel time in seconds along the ray.
    [[nodiscard]] double getTravelTime() const;
    /// @result The total length of the ray in meters.
    [[nodiscard]] double getLength() const;
    /// @result The takeoff angle in degrees measured positive up from nadir. 
    /// @throws std::runtime_error if \c size() is 0.
    [[nodiscard]] double getTakeOffAngle() const;

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Path2D();
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] path   The path to copy to this.
    /// @result A deep copy of the path. 
    Path2D &operator=(const Path2D &path);
    /// @brief Move assignment operator.
    /// @result The memory from path moved to this.
    Path2D &operator=(Path2D &&path) noexcept;
    /// @result A reference to the first ray segment in the path.
    iterator begin();
    /// @result A constant reference to the first ray segment in the path.
    const_iterator begin() const noexcept;
    /// @result A constant reference to the first ray segment in the path.
    const_iterator cbegin() const noexcept;

    /// @result A reference to the last ray segment in the path. 
    iterator end(); 
    /// @result A reference to the last ray segment in the path.
    const_iterator end() const noexcept;
    /// @result A reference to the last ray segment in the path.
    const_iterator cend() const noexcept;

    /// @param[in] index  The index of the desired segment.
    /// @result A reference to the segment at the given position.
    /// @throws std::invalid_argument if this is out of bounds.
    Segment2D& at(size_t index);
    /// @param[in] index  The index of the desired segment.
    /// @result A reference to the segment at the given position.
    /// @throws std::invalid_argument if this is out of bounds.
    const Segment2D& at(size_t index) const;
    /// @}
private:
    class Path2DImpl;
    std::unique_ptr<Path2DImpl> pImpl;
};
}
#endif
