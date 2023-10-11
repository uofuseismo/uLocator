#ifndef ULOCATOR_STATIC_CORRECTION_HPP
#define ULOCATOR_STATIC_CORRECTION_HPP
#include <memory>
namespace ULocator
{
class StaticCorrection
{
public:
    enum class Method
    {
        Mean,   /*!< Weighted mean correction.  This is appropriate for 
                     least-squares. */
        Median  /*!< Weighted median correction.  This is apprporiate for
                     L1-norm inversions. */
    }; 
public:
    StaticCorrection();
    StaticCorrection(const StaticCorrection &correction);
    StaticCorrection(StaticCorrection &&correction) noexcept;

    void initialize(const std::string &network,
                    const std::string &station,
                    const std::string &phase);
    void initialize(const std::string &network,
                    const std::string &station,
                    const std::string &phase,
                    double correction);
    /// @result True indicates the class is initialized. 
    [[nodiscard]] bool isInitialized() const noexcept;

    /// @brief Sets the static correction in seconds.
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
    /// @param[in] predictedArrivalTiems  The predicted arrival tiems in seconds
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
    /// @param[in] predictedArrivalTiems  The predicted arrival tiems in seconds
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
    /// @brief Saves the static correction to an HDF5 archive.
    /// @param[in] fileName  The name of the HDF5 archive.
    /// @throws std::runtime_error if \c haveCorrection() is \c isInitialized()
    ///         is false.  
    void save(const std::string &fileName) const;
    /// @brief Loads the static correction from an HDF5 archive.
    /// @param[in] fileName  The name of the HDF5 archive.
    /// @throws std::runtime_error if \c isInitialized() is false.
    /// @throws std::invalid_argument if the HDF5 file doesn't exist or 
    ///         contain the correction.
    void load(const std::string &fileName);

    /// @result Evaluates the model - i.e., returns the `corrected' prediction
    ///         time.
    [[nodiscard]] double evaluate(double predictedTime) const;

    /// @result A deep copy of the static correction.
    StaticCorrection& operator=(const StaticCorrection &correction);
    /// @result The memory from correction moved to this.
    StaticCorrection& operator=(StaticCorrection &&correction) noexcept;
    /// @brief Evalutes the static correction (in seconds) to add to the
    ///        the estimates.
    /// @throws std::runtime_error if \c haveCorrection() is false. 
    [[nodiscard]] double operator()(double predictedTime) const;

    /// @brief Destructor.
    ~StaticCorrection();

private:
    class StaticCorrectionImpl;
    std::unique_ptr<StaticCorrectionImpl> pImpl;
};
}
#endif
