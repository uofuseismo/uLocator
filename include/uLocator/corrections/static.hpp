#ifndef ULOCATOR_CORRECTIONS_STATIC_HPP
#define ULOCATOR_CORRECTIONS_STATIC_HPP
#include <memory>
namespace ULocator::Corrections
{
/// @class Static "static.hpp" "uLocator/corrections/static.hpp"
/// @brief Defines a static correction for a given NETWORK.STATION and
///        phase.  This is a constant value that is always added to the
///        the predicted travel time.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class Static
{
public:
    /// @brief Defines the static correction computation strategy.
    enum class Method
    {
        Mean,   /*!< Weighted mean correction.  This is appropriate for 
                     least-squares. */
        Median  /*!< Weighted median correction.  This is apprporiate for
                     L1-norm inversions. */
    }; 
public:
    /// @name Constructors
    /// @{
 
    /// @brief Constructor.
    Static();
    /// @brief Copy constructor.
    /// @param[in] correction  The correction from which to initialize this
    ///                        class.
    Static(const Static &correction);
    /// @brief Move constructor.
    /// @param[in] correction  The correction from which to initialize this
    ///                        class.  On exit, correction's behavior is
    ///                        undefined.
    Static(Static &&correction) noexcept;
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

    /// @name Defining the Correction
    /// @{

    /// @brief Sets the static correction.
    /// @param[in] correction  The correction in seconds.  This is the time
    ///                        to add to the estimated travel time.
    void setCorrection(double correction) noexcept;
    /// @result The static correction in seconds to add to the estimated time.
    /// @throw std::runtime_error if \c haveCorrection() is false.
    [[nodiscard]] double getCorrection() const;
    /// @result True indicates the static correction was set.
    [[nodiscard]] bool haveCorrection() const noexcept;

    /// @brief Computes the static correction assuming all arrivals are evenly
    ///        weighted.
    /// @param[in] observedArrivalTimes   The observed arrival times in seconds 
    ///                                   since the epoch.
    /// @param[in] predictedArrivalTimes  The predicted arrival times in seconds
    ///                                   since the epoch.
    /// @param[in] method   Defines the data fitting method.
    /// @throws std::invalid_argument if the input vectors are of unequal size
    ///         or are empty.
    void train(const std::vector<double> &observedArrivalTimes,
               const std::vector<double> &predictedArrivalTimes,
               Method method);
    /// @brief Computes the static correction for unevenly weighted arrivals.
    /// @param[in] observedArrivalTimes   The observed arrival times in seconds 
    ///                                   since the epoch.
    /// @param[in] predictedArrivalTimes  The predicted arrival times in seconds
    ///                                   since the epoch.
    /// @param[in] weights  The weights on the observations.  Typically 1 is
    ///                     maximal and zero 0 is minimal.  
    /// @param[in] method   Defines the data fitting method.
    /// @throws std::invalid_argument if the input vectors are of unequal size,
    ///         are empty, or all the weights are zero.
    void train(const std::vector<double> &observedArrivalTimes,
               const std::vector<double> &predictedArrivalTimes,
               const std::vector<double> &weights,
               Method method); 
    /// @}

    /// @name Saving/Loading Model
    /// @{

    /// @brief Saves the static correction to an HDF5 archive.
    /// @param[in] fileName  The name of the HDF5 archive.
    /// @throws std::runtime_error if \c haveCorrection() is false
    ///         \c haveStationName() is false, or \c compiledWithHDF5() is
    ///         false.
    void save(const std::string &fileName) const;
    /// @brief Loads the static correction from an HDF5 archive.
    /// @param[in] fileName  The name of the HDF5 archive.
    /// @throws std::runtime_error if \c haveStationNameAndPhase() is false or
    ///         \c compiledWithHDF5() is false.
    /// @throws std::invalid_argument if the HDF5 file doesn't exist or 
    ///         contain the correction.
    void load(const std::string &fileName);
    /// @result True indicates the library was compiled with HDF5 support. 
    [[nodiscard]] static bool compiledWithHDF5() noexcept; 
    /// @}

    /// @name Model Evaluation
    /// @{

    /// @param[in] predictedTime  The predicted travel time in seconds.
    /// @result The `corrected' predicted travel time in seconds - i.e.,
    ///         correctedTime = predictedTime + getCorrection().
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double evaluate(double predictedTime) const;
    /// @param[in] predictedTime  The predicted travel time in seconds.
    /// @param[out] dtdt0   The derivative of the correction with respect to
    ///                     the origin time.
    /// @param[out] dtdx    The derivative of the correction with respect to
    ///                     the source x position.
    /// @param[out] dtdy    The derivative of the correction with respect to
    ///                     the source y position.
    /// @param[out] dtdy    The derivative of the correction with respect to
    ///                     the source z position.
    /// @result The `corrected' predicted travel time in seconds - i.e.,
    ///         correctedTime = predictedTime + getCorrection().
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double evaluate(double predictedTime, double *dtdt0, double *dtdx, double *dtdy, double *dtdz) const;
    /// @}

    /// @name Operators
    /// @{
 
    /// @result A deep copy of the static correction.
    Static& operator=(const Static &correction);
    /// @result The memory from correction moved to this.
    Static& operator=(Static &&correction) noexcept;
    /// @param[in] predictedTime  The predicted travel time in seconds.
    /// @result The `corrected' predicted travel time in seconds - i.e.,
    ///         correctedTime = predictedTime + getCorrection().
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double operator()(double predictedTime) const;
    /// @param[in] predictedTime  The predicted travel time in seconds.
    /// @param[out] dtdt0   The derivative of the correction with respect to
    ///                     the origin time.
    /// @param[out] dtdx    The derivative of the correction with respect to
    ///                     the source x position.
    /// @param[out] dtdy    The derivative of the correction with respect to
    ///                     the source y position.
    /// @param[out] dtdy    The derivative of the correction with respect to
    ///                     the source z position.
    /// @result The `corrected' predicted travel time in seconds - i.e.,
    ///         correctedTime = predictedTime + getCorrection().
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double operator()(double predictedTime, double *dtdt0, double *dtdx, double *dtdy, double *dtdz) const;
    /// @}

    /// @name Destructors
    /// @{

    /// @brief Resets the class and releases memory.
    void clear() noexcept;
    /// @brief Destructor.
    ~Static();
    /// @}
private:
    class StaticImpl;
    std::unique_ptr<StaticImpl> pImpl;
};
}
#endif
