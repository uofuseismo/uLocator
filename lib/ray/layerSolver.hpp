#ifndef ULOCATOR_EIKONALXX_RAY_LAYER_SOLVER_HPP
#define ULOCATOR_EIKONALXX_RAY_LAYER_SOLVER_HPP
#include <memory>
namespace EikonalXX::Ray
{
class Path2D;
}
namespace EikonalXX::Ray
{
/// @brief A simple utility class for computing rays for a first arrival
///        in 1D, isotropic, layered media.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class LayerSolver
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    LayerSolver();
    /// @brief Move constructor.
    LayerSolver(LayerSolver &&solver) noexcept;
    /// @}
 
    /// @brief Sets a layer cake model of the form:
    ///        Interface 0 ---------------------------- (Free Surface)
    ///                    |        Velocity 0        |
    ///        Interface 1 ----------------------------
    ///                    |        Velocity 1        |
    ///                    ----------------------------
    ///                                 .
    ///                                 .
    ///                                 .
    ///        Interface n --------------------------- (Underlying half space)
    ///                    |        Velocity n       |
    ///                    ---------------------------
    ///
    /// @param[in] interfaces     The interfaces for each layer.
    /// @param[in] velocityModel  The velocities in each layer.
    /// @note Velocity inversions are not yet considered.
    void setVelocityModel(const std::vector<double> &interfaces,
                          const std::vector<double> &velocityModel);
    /// @result True indicates the velocity model was set.
    [[nodiscard]] bool haveVelocityModel() const noexcept;

    /// @brief Sets the source depth.
    /// @param[in] depth  The source depth in meters.
    /// @throws std::invalid_argument if the source is not in the layer cake
    ///         model.
    /// @throws std::runtime_error if \c haveVelocityModel() is false.
    void setSourceDepth(const double depth);
    /// @result The source depth in meters.
    [[nodiscard]] double getSourceDepth() const;
    /// @result True indicates the source depth was set.
    [[nodiscard]] bool haveSourceDepth() const noexcept;

    /// @brief Sets the station offset, in meters, with the station at the
    ///        free surface.
    void setStationOffset(double offset);
    /// @brief Sets the station offset, in meters, and the station depth
    ///        in meters.
    void setStationOffsetAndDepth(double offset, double depth);
    [[nodiscard]] bool haveStationOffsetAndDepth() const noexcept;

    /// @brief Shoots rays through the model at the given take-off angle.
    /// @param[in] takeOffAngle  The take-off angle at which to shoot in
    ///                          degrees measured positive up from nadir.
    /// @param[in] keepOnlyHits  If true, then only ray paths that `hit'
    ///                          i.e., finish within some offset of the
    ///                          station will be retained.
    /// @result The ray paths shot at this angle.
    [[nodiscard]] std::vector<Path2D> shoot(double takeOffAngle, bool keepOnlyHits = false);

    /// @brief Traces rays from the source to the receiver through the 1D
    ///        medium.
    /// @param[in] doReflections   If true then the layer reflections will also
    ///                            be tabulated.  However, these are never the
    ///                            first arriving energy so this may be of
    ///                            limited interest.
    void solve(bool doReflections = false);
 
    /// @result The ray paths from the source-to-receiver sorted in ascending
    ///         order.
    /// @throws std::runtime_error if \c haveRayPaths() is false.
    [[nodiscard]] std::vector<Path2D> getRayPaths() const;
    [[nodiscard]] const std::vector<Path2D> &getRayPathsReference() const;
    /// @result True indicates the ray path is available.
    [[nodiscard]] bool haveRayPaths() const noexcept;

    /// @name Destructors
    /// @{

    /// @brief Reset the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~LayerSolver();
    /// @}

    LayerSolver& operator=(LayerSolver &&solver) noexcept;

    LayerSolver(const LayerSolver &) = delete;
    LayerSolver& operator=(const LayerSolver &) = delete;
private:
    class LayerSolverImpl;
    std::unique_ptr<LayerSolverImpl> pImpl;
};
}
#endif
