#ifndef ULOCATOR_TOPOGRAPHY_GRIDDED_HPP
#define ULOCATOR_TOPOGRAPHY_GRIDDED_HPP
#include <memory>
#include "uLocator/topography/topography.hpp"
namespace ULocator::Position
{
 class IGeographicRegion;
}
namespace ULocator::Topography
{
/// @class Gridded "gridded.hpp" "uLocator/topography/gridded.hpp"
/// @brief This is a topography function interpolator that operates on
///        a regular grid.  Effectively, it takes a latitude, longitude pair
///        and returns the linearly interpolated elevation e.g.,
///        elevation = topograhy(latitude, longitude)
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Gridded : public ITopography
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Gridded();
    /// @brief Copy constructor.
    /// @param[in] topography   The class from which to initialize this class.
    Gridded(const Gridded &topography);
    /// @brief Move constructor.
    /// @brief
    Gridded(Gridded &&topography) noexcept;
    /// @}

    /// @brief Loads the topography model from an HDF5 file.  Since the
    ///        topography comes from Etopo we get (lat, lon, elevation) 
    ///        tuples so this routine will convert the topographic model
    ///        to (x, y, elevation).
    /// @param[in] hdf5FileName  The name of the HDF5 file with the topogrpahy.
    /// @param[in] region        The geographic region to transform 
    ///                          from (lat, lon) to local (x, y) coordinates.
    /// @throws std::invalid_argument if the file does not exist or is not
    ///         properly defined.
    /// @throws std::runtime_error the code was not compiled with HDF5. 
    void load(const std::string &hdf5FileName, const ULocator::Position::IGeographicRegion &region);
    /// @brief Sets the topography.
    /// @param[in] nx          The number of x grid points.
    /// @param[in] x           The positions of the x grid points in meters.
    /// @param[in] ny          The number of y grid points.
    /// @param[in] y           The positions of the y grid points in meters.
    /// @param[in] nGrid       The number of grid points.  This must equal nx*ny.
    /// @param[in] elevation   The elevation in meters with respect to
    ///                        sea-level at the specified (x, y).
    ///                        This is an array whose dimension is
    ///                        [nx x ny] with leading nx.
    template<typename U>
    void set(int nx, const U x[],
             int ny, const U y[],
             int nGrid, const U elevation[]);
    /// @result True indicates the topography is set.
    [[nodiscard]] bool haveTopography() const noexcept final;
    /// @result The elevation with respect to sea-level in meters.  Note,
    /// @param[in] x    The x position in meters.
    /// @param[in] y    The y position in meters.
    /// @throws std::runtime_error if \c haveTopography() is false.
    [[nodiscard]] double evaluate(double x, double y) const final;
    /// @result The elevation with respect to sea-level in meters.  Note,
    /// @param[in] x    The x position in meters.
    /// @param[in] y    The y position in meters.
    /// @param[out] dElevationDx  The change in elevation with respect to x.
    /// @param[out] dElevationDy  The change in elevation with respect to y. 
    /// @throws std::runtime_error if \c haveTopography() is false.
    [[nodiscard]] double evaluate(double x, double y, double *dElevationDx, double *dElevationDy) const final;
    /// @result The minimum and maximum elevation with respect to sea-level
    ///         in meters.
    [[nodiscard]] std::pair<double, double> getMinimumAndMaximumElevation() const final;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] topography  The topography to copy to this.
    /// @result A deep copy of the topography.
    Gridded& operator=(const Gridded &topography);
    /// @brief Move assignment.
    /// @param[in,out] topography  The topography whose memory will be moved 
    ///                            to this.  On exit, topography's behavior
    ///                            is undefined.
    /// @result The memory from topography moved to this.
    Gridded& operator=(Gridded &&topography) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    virtual ~Gridded();
    /// @}
private:
    class GriddedImpl; 
    std::unique_ptr<GriddedImpl> pImpl;
};
}
#endif
