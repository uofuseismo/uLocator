#ifndef ULOCATOR_POSITION_YNP_HPP
#define ULOCATOR_POSITION_YNP_HPP
#include <memory>
#include <uLocator/position/geographicRegion.hpp>
namespace ULocator::Position
{
/// @brief Defines the YNP geographic region for optimization.
///        This region is approximately -114.75 to -108.25 degrees
///        longitude and approximately 36.25 to 43 degrees latitude.
/// @copyright Ben Baker (University of YNP) distributed under the MIT license.
class YNP : public IGeographicRegion
{
public:
    /// @name Constructors.
    /// @{

    /// @brief Constructor.
    YNP();
    /// @brief Copy constructor.
    /// @param[in] ynp  The YNP region from which to initialize this class.
    YNP(const YNP &ynp);
    /// @brief Move constructor.
    /// @param[in,out] ynp  The YNP region from which to initialize this class.
    ///                     On exit, YNP's behavior is undefined.
    YNP(YNP &&ynp) noexcept;
    /// @}

    /// @name Destructors.
    /// @{

    /// @brief Destructor.
    ~YNP();
    /// @}

    /// @result Converts the given latitude, longitude pair to a local x, y coordinate in meters.
    [[nodiscard]] std::pair<double, double> geodeticToLocalCoordinates(double latitude, double longitude) const final;
    /// @result Converts the given local (x, y) coordinates in metres to latitude and longitude in degrees.
    [[nodiscard]] std::pair<double, double> localToGeodeticCoordinates(double x, double y) const final;

    /// @result The minimum and maximum x extent of the YNP search region in 
    ///         meters in the local coordinate system.
    [[nodiscard]] std::pair<double, double> getExtentInX() const noexcept final;
    /// @result The minimum and maximum y extent of the YNP search region in
    ///         meters in the local coordinate system. 
    [[nodiscard]] std::pair<double, double> getExtentInY() const noexcept final;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] ynp  The YNP region to copy to this.
    /// @result A deep copy of the YNP boundaries.
    YNP& operator=(const YNP &ynp);
    /// @brief Move assignment.
    /// @param[in,out] ynp  The YNP region whose memory will be moved to this.
    ///                     On exit, YNP's behavior is undefined.
    /// @result The memory from the YNP boundaries moved to this.
    YNP& operator=(YNP &&ynp) noexcept;
    /// @}
private:
    class YNPImpl;
    std::unique_ptr<YNPImpl> pImpl;
};
}
#endif