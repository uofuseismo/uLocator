#ifndef ULOCATOR_TOPOGRAPHY_HPP
#define ULOCATOR_TOPOGRAPHY_HPP
#include <memory>
namespace ULocator
{
class Topography
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    Topography();
    /// @brief Copy constructor.
    /// @param[in] topography   The class from which to initialize this class.
    Topography(const Topography &topography);
    /// @brief Move constructor.
    /// @brief
    Topography(Topography &&topography) noexcept;
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
    [[nodiscard]] bool haveTopography() const noexcept;
    /// @result The elevation
    [[nodiscard]] double evaluate(double latitude, double longitude) const;

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] topography  The topography to copy to this.
    /// @result A deep copy of the topography.
    Topography& operator=(const Topography &topography);
    /// @brief Move assignment.
    /// @param[in,out] topography  The topography whose memory will be moved 
    ///                            to this.  On exit, topography's behavior
    ///                            is undefined.
    /// @result The memory from topography moved to this.
    Topography& operator=(Topography &&topography) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Topography();
    /// @}
private:
    class TopographyImpl; 
    std::unique_ptr<TopographyImpl> pImpl;
};
}
#endif
