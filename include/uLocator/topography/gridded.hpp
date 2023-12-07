#ifndef ULOCATOR_TOPOGRAPHY_GRIDDED_HPP
#define ULOCATOR_TOPOGRAPHY_GRIDDED_HPP
#include <memory>
#include "uLocator/topography/topography.hpp"
namespace ULocator::Topography
{
/// @class Gridded "topgraphy.hpp" "uLocator/topography.hpp"
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

    /// @brief Loads the topography model from an HDF5 file.
    /// @param[in] hdf5FileName  The name of the HDF5 file with the topogrpahy.
    /// @throws std::invalid_argument if the file does not exist or is not
    ///         properly defined.
    /// @throws std::runtime_error the code was not compiled with HDF5. 
    void load(const std::string &hdf5FileName);
    /// @brief Sets the topography.
    /// @param[in] nLatitudes   The number of latitudes in the grid.
    /// @param[in] latitudes    The latitudes in degrees in the grid.  This
    ///                         is an array whose dimension is [nLatitudes].
    /// @param[in] nLongitudes  The number of longitudes in teh grid.  This
    ///                         is an array whose dimension is [nLongitudes].
    /// @param[in] elevation    The elevation in meters with respect to
    ///                         sea-level at the specified (latitude,longitude).
    ///                         This is an array whose dimension is
    ///                         [nLatitudes x nLongitudes] with leading
    ///                         dimension nLongitudes.
    template<typename U>
    void set(int nLatitudes,  const U latitudes[],
             int nLongitudes, const U longitudes[],
             int nGrid, const U elevation[]);
    /// @result True indicates the topography is set.
    [[nodiscard]] bool haveTopography() const noexcept final;
    /// @result The elevation with respect to sea-level in meters.  Note,
    /// @param[in] latitude   The latitude in degrees.
    /// @param[in] longitude  The longitude in degrees.
    /// @throws std::invalid_argument if the latitude is not in the
    ///         range [-90,90].
    /// @throws std::runtime_error if \c haveGridded() is false.
    [[nodiscard]] double evaluate(double latitude, double longitude) const final;

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
