#ifndef ULOCATOR_CORRECTIONS_SOURCE_SPECIFIC_HPP
#define ULOCATOR_CORRECTIONS_SOURCE_SPECIFIC_HPP
#include <memory>
namespace ULocator::Corrections
{
/// @class SourceSpecific "sourceSpecific.hpp" "uLocator/corrections/sourceSpecific.hpp"
/// @brief For a given station and phase these corrections memorize the
///        residuals for given sources and use these residuals to estimate
///        a correction using a K-Nearest Neighbor technique.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class SourceSpecific
{
public:
    /// @brief Optimization algorithm.
    enum class OptimizationAlgorithm
    {
         KDTree = 0, /*!< Use a KD-tree.  This is way faster than brute froce. */
         BruteForce  /*!< Compute nearest neighbors with brute force. */
    };
    /// @brief This defines the interpolation scheme used during model
    ///        evaluation.
    enum class EvaluationMethod : int
    {
         WeightedAverage = 0,  /*!< Takes the weighted average of the 
                                    k-nearest neighbors where the weights
                                    are the pick qualities. */
         MaximumDistanceWeighted = 1, /*!< Takes the weighted average of the
                                           k-nearest neighbors, however, if
                                           a neighbor exceeds a certain 
                                           distance then its contribution is
                                           zero. */
         InverseDistanceWeighted = 2  /*!< Uses inverse distance weighted of 
                                           the k-nearest neighbors and weights
                                           based on pick qualities.  
                                           Additionally, if a neighbor exceeds
                                           a certain distance then the 
                                           contribution is zero. */
    };
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    SourceSpecific();
    /// @brief Copy constructor.
    /// @param[in] correction  The correction class from which to initialize
    ///                        this class.
    SourceSpecific(const SourceSpecific &correction);
    /// @brief Move constructor.
    /// @param[in,out] correction  The correction class from which to initialize
    ///                            this class.  On exit, correction's behavior
    ///                            is undefined.
    SourceSpecific(SourceSpecific &&correction) noexcept;
    /// @}

    /// @name Station and Phase Identification
    /// @{

    /// @brief Defines the station name and phase to which this correction
    ///        was computed.
    /// @param[in] network  The network code - e.g., UU.
    /// @param[in] station  The station name - e.g., FORK.
    /// @param[in] phase    The seismic phase - e.g., P or S.
    void setStationNameAndPhase(const std::string &network,
                                const std::string &station,
                                const std::string &phase);
    /// @result True indicates the station and phase were set.
    [[nodiscard]] bool haveStationNameAndPhase() const noexcept;
    /// @result The network code.
    /// @throws std::runtime_error if \c haveStationNameAndPhase() is false.
    [[nodiscard]] std::string getNetwork() const;
    /// @result The station name.
    /// @throws std::runtime_error if \c haveStationNameAndPhase() is false.
    [[nodiscard]] std::string getStation() const;
    /// @result The seismic phase.
    /// @throws std::runtime_error if \c haveStationNameAndPhase() is false.
    [[nodiscard]] std::string getPhase() const;
    /// @}

    /// @name Hyper-Parameters
    /// @{

    /// @brief Sets the number of neighbors in the KNN interpolation scheme.
    /// @param[in] numberOfNeighbors  The number of neighbors in the KNN
    ///                               training scheme.
    /// @throws std::invalid_argument if this is not positive.
    void setNumberOfNeighbors(int numberOfNeighbors);
    /// @result The number of neighbors in the KNN scheme.
    /// @note If the model is trained, saved, then loaded then the number of
    ///      neighbors will also be set.
    [[nodiscard]] int getNumberOfNeighbors() const noexcept; 

    /// @brief Sets the maximum distance beyond which a source's residual
    ///        will not contribute to the interpolation.  For example,
    ///        if there are five nearest neighbors, but the distance of two
    ///        neighbors exceed this distance, then the residuals of those
    ///        two neighbors will not be used for interpolation.
    /// @param[in] maximumDistance  The maximum distance in meters.
    /// @note This is only relevant if the evaluation method uses 
    ///       EvaluationMethod::MaximumDistanceWeighted or 
    ///       EvaluationMethod::InverseDistanceWeighted weighting.
    /// @throws std::invalid_argument if this is not positive.
    void setMaximumDistance(double maximumDistance); 
    /// @result The maximum cutoff distance in meters when the interpolation
    ///         method is EvaluationMethod::MaximumDistanceWeighted or
    ///         EvaluationMethod::InverseDistanceWeighted.
    /// @note If the model is trained using the MaximumDistanceWeighted or
    ///       InverseDistanceWeighted evaluation method, saved, then loaded
    ///       then this will be set.
    /// @throws std::runtime_error if \c haveMaximumDistance() is false.
    [[nodiscard]] double getMaximumDistance() const;
    /// @result True indicates the maximum distance was set.
    [[nodiscard]] bool haveMaximumDistance() const noexcept;

    /// @brief Sets the evaluation method.
    /// @param[in] evalulationMethod  The evaluation method.
    void setEvaluationMethod(EvaluationMethod evaluationMethod) noexcept;
    /// @result The evaluation method.
    /// @note If the model is trained, saved, then loaded then the evaluation
    ///       method will be set corresponding to that training procedure.
    [[nodiscard]] EvaluationMethod getEvaluationMethod() const noexcept;

    /// @brief Optimization algorithm
    void setOptimizationAlgorithm(OptimizationAlgorithm algorithm) noexcept;
    /// @result The optimization algorithm.
    [[nodiscard]] OptimizationAlgorithm getOptimizationAlgorithm() const noexcept;
    /// @}
     
    /// @name Training
    /// @{

    /// @brief Fits the KNN model.
    /// @param[in] xSourcePositions    The x position of the sources in meters.
    /// @param[in] ySourcePositions    The y position of the sources in meters.
    /// @param[in] zSourcePositions    The z position of the sources in meters.
    /// @param[in] observedArrivalTimes   The observed arrival times in seconds 
    ///                                   since the epoch.
    /// @param[in] predictedArrivalTimes  The predicted arrival tiems in seconds
    ///                                   since the epoch.
    /// @param[in] weights  The weights on the observations.  Typically 1 is
    ///                     maximal and zero 0 is minimal.  
    /// @throws std::invalid_argument if the input vectors are of unequal size,
    ///         are empty, or all the weights are zero.
    void train(const std::vector<double> &xSourcePositions,
               const std::vector<double> &ySourcePositions,
               const std::vector<double> &zSourcePositions,
               const std::vector<double> &observedArrivalTimes,
               const std::vector<double> &predictedArrivalTimes,
               const std::vector<double> &weights);
    /// @result True indicates the model is trained.
    [[nodiscard]] bool haveModel() const noexcept;
    /// @}

    /// @name Saving/Loading Model
    /// @{

    /// @brief Saves the model to an HDF5 file.
    /// @param[in] fileName  The name of the HDF5 file.
    /// @throws std::runtime_error if \c haveModel() is false
    ///         \c haveStationNameAndPhase() is false, or \c compiledWithHDF5()
    ///         is false.
    void save(const std::string &fileName) const;
    /// @result True indicates the library was compiled with HDF5 support. 
    [[nodiscard]] static bool compiledWithHDF5() noexcept; 

    /// @brief Loads the KNN model from an HDF5 archive.
    /// @param[in] fileName  The name of the HDF5 archive.
    /// @throws std::runtime_error if \c haveStationNameAndPhase() is false or
    ///         \c compiledWithHDF5() is false.
    /// @throws std::invalid_argument if the HDF5 file doesn't exist or 
    ///         contain the correction.
    void load(const std::string &fileName);
    /// @}

    /// @name Model Evaluation
    /// @{

    /// @param[in] xSource  The x position of the source in meters.
    /// @param[in] ySource  The y position of the source in meters.
    /// @param[in] zSource  The z position of the source in meters.
    /// @param[in] predictedTime  The predicted travel time in seconds.
    /// @result The `corrected' predicted travel time in seconds - i.e.,
    ///         correctedTime = predictedTime + KNN(xSource, ySource, zSource).
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double evaluate(double xSource,
                                  double ySource,
                                  double zSource,
                                  double predictedTime = 0) const;
    /// @param[in] xSource  The x position of the source in meters.
    /// @param[in] ySource  The y position of the source in meters.
    /// @param[in] zSource  The z position of the source in meters.
    /// @param[out] dtdt0   The derivative of the correction with respect to
    ///                     the origin time.
    /// @param[out] dtdx    The derivative of the correction with respect to
    ///                     the source x position.
    /// @param[out] dtdy    The derivative of the correction with respect to
    ///                     the source y position.
    /// @param[out] dtdz    The derivative of the correction with respect to
    ///                     the source z position.
    /// @param[in] predictedTime  The predicted travel time in seconds.
    /// @result The `corrected' predicted travel time in seconds - i.e.,
    ///         correctedTime = predictedTime + getCorrection().
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double evaluate(double xSource, 
                                  double ySource,
                                  double zSource,
                                  double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
                                  double predictedTime = 0) const;
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] correction  The correction to copy to this.
    /// @result A deep copy of the correction. 
    SourceSpecific& operator=(const SourceSpecific &correction);
    /// @brief Move assignment.
    /// @param[in,out] correction  The correction whose memory will be moved to
    ///                            this.  On exit, correction's behavior is
    ///                            undefined.
    /// @result The memory from the correction moved to this.
    SourceSpecific& operator=(SourceSpecific &&correction) noexcept;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Default destructor.
    ~SourceSpecific();
    /// @}
private:
    class SourceSpecificImpl;
    std::unique_ptr<SourceSpecificImpl> pImpl;
};
}
#endif
