#ifndef ULOCATOR_SOURCE_SPECIFIC_STATION_CORRECTION_HPP
#define ULOCATOR_SOURCE_SPECIFIC_STATION_CORRECTION_HPP
#include <memory>
namespace ULocator
{
class SourceSpecificStationCorrection
{
public:
    enum class EvaluationMethod
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
    SourceSpecificStationCorrection();
    SourceSpecificStationCorrection(SourceSpecificStationCorrection &&correction) noexcept;
    SourceSpecificStationCorrection& operator=(SourceSpecificStationCorrection &&correction) noexcept;

    void initialize(const std::string &network,
                    const std::string &station,
                    const std::string &phase);
    [[nodiscard]] bool isInitialized() const noexcept;
    /// @brief Fits the KNN model.
    /// @param[in] observedArrivalTimes   The observed arrival times in seconds 
    ///                                   since the epoch.
    /// @param[in] predictedArrivalTiems  The predicted arrival tiems in seconds
    ///                                   since the epoch.
    /// @param[in] weights  The weights on the observations.  Typically 1 is
    ///                     maximal and zero 0 is minimal.  
    /// @param[in] maximumDistance  The maximum distance in meters.  If this
    ///                             is negative then the default distance will
    ///                             be used.
    /// @throws std::invalid_argument if the input vectors are of unequal size,
    ///         are empty, or all the weights are zero.
    void train(const std::vector<double> &latitudes,
               const std::vector<double> &longitudes,
               const std::vector<double> &depths,
               const std::vector<double> &observedArrivalTimes,
               const std::vector<double> &predictedArrivalTimes,
               const std::vector<double> &weights,
               int nNeighbors = 5,
               double maximumDistance =-1);
    [[nodiscard]] int getNumberOfNeighbors() const;
    [[nodiscard]] double getMaximumDistance() const;
    /// @result True indicates the model is trained.
    [[nodiscard]] bool haveModel() const noexcept;
    /// @brief Saves the model to an HDF5 file.
    /// @param[in] fileName  The name of the HDF5 file.
    /// @throws std::runtime_error if \c isInitialized() is false or
    ///         \c haveModel() is false. 
    void save(const std::string &fileName) const;
    /// @brief Loads the KNN model from an HDF5 archive.
    /// @param[in] fileName  The name of the HDF5 archive.
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @throws std::invalid_argument if the HDF5 file doesn't exist or 
    ///         contain the training data.
    void load(const std::string &fileName);

    [[nodiscard]] double evaluate(double latitude,
                                  double longitude,
                                  double depth,
                                  EvaluationMethod method = EvaluationMethod::InverseDistanceWeighted) const;
    [[nodiscard]] std::vector<double> evaluate(const std::vector<double> &latitudes,
                                               const std::vector<double> &longitudes,
                                               const std::vector<double> &depths,
                                               EvaluationMethod method = EvaluationMethod::InverseDistanceWeighted) const;

    ~SourceSpecificStationCorrection();

    SourceSpecificStationCorrection& operator=(const SourceSpecificStationCorrection &) = delete;
    SourceSpecificStationCorrection(const SourceSpecificStationCorrection &) = delete;
private:
    class SourceSpecificStationCorrectionImpl;
    std::unique_ptr<SourceSpecificStationCorrectionImpl> pImpl;
};
}
#endif
