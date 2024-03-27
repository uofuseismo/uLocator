#ifndef ULOCATOR_OPTIMIZERS_ORIGIN_TIME_HPP
#define ULOCATOR_OPTIMIZERS_ORIGIN_TIME_HPP
#include <memory>
#include <uLocator/optimizers/optimizer.hpp>
namespace ULocator
{
class Origin;
}
namespace ULocator::Optimizers
{
/// @class OriginTime "originTime.hpp" "uLocator/optimizers/originTime.hpp"
/// @brief Optimizes origin times.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class OriginTime
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    OriginTime();
    /// @brief Copy constructor.
    OriginTime(const OriginTime &originTime);
    /// @brief Move constructor.
    OriginTime(OriginTime &&originTime) noexcept;
    /// @}
 
    /// @brief Sets the norm.   
    /// @param[in] norm  The norm in which to optimize.
    /// @param[in] p     For an Lp norm optimization this is the value of p.
    ///                  This is ignored if norm is not Lp.
    void setNorm(IOptimizer::Norm norm, double p = 1.5);
    /// @result The norm in which to optimize the origin time.
    [[nodiscard]] IOptimizer::Norm getNorm() const noexcept;

    /// @brief For a p-norm optimization - if the bracketing window's width is
    ///        less than this value, in seconds, then we will terminate the
    ///        the search.
    /// @param[in] tolerance  The window tolerance in the golden section search.
    ///                       Note, if this is negative then the tolerance
    ///                       check is disabled and the goldern search will
    ///                       continue to its maximum number of iterations.
    void setTolerance(double tolerance);
    /// @result The tolerance for the window.
    [[nodiscard]] double getTolerance() const noexcept;
    /// @brief For a p-norm optimization - the golden section search will
    ///        terminate after this many iterations.
    /// @param[in] nIterations  The maximum number of golden section search
    ///                         iterations to perform.
    /// @throws std::invalid_argument if this is not positive.
    void setNumberOfIterations(int nIterations);
    /// @result The number of iterations.
    [[nodiscard]] int getNumberOfIterations() const noexcept;

    /// @brief When using the golden section search, we create an initial,
    ///        window from this many seconds to before the first arrival
    ///        up until the initial arrival time.  This is a tuning
    ///        variable that can be estimated by looking at the maximum
    ///        first arrival - origin time in your catalog then adding 
    ///        in a safety tolerance.
    /// @param[in] window  The number of seconds before the first arrival
    ///                    to begin the golden section bracket.
    /// @throws std::invalid_argument if this is not positive. 
    void setTimeWindow(double window);
    /// @result The initial window, in seconds before the first arrival, to
    ///         search with the golden section. 
    [[nodiscard]] double getTimeWindow() const noexcept;

    /// @brief This is a general strategy when dealing with epochal times
    ///        to make the smallest arrival time the `0'.  This can help
    ///        with adding meaning to the tolerance when dealing with an
    ///        iterative method.
    void enableTimeReduction() noexcept;
    /// @brief Calculations will be performed with reducing arrival times. 
    void disableTimeReduction() noexcept;
    /// @result True indicates calculations will be performed with reduced
    ///         arrival times.
    /// @note The default is true.
    [[nodiscard]] bool reduceTimes() const noexcept;

    /// @brief Sets the arrival times.  The weights will be lifted from the
    ///        arrival class.
    /// @param[in] arrivals  The arrival times and weights to set.
    void setArrivalTimes(const std::vector<ULocator::Arrival> &arrivals);
    /// @brief Sets the arrival times and weights.
    /// @param[in] arrivalTimes  The arrival times, UTC seconds since the epoch,
    ///                          to set.
    /// @param[in] weights       The corresponding weights.
    /// @throws std::invalid_argument if arrivalTimes.size() != weights.sizse()
    ///         or any weight is not positive or there are no arrival times.
    void setArrivalTimes(const std::vector<double> &arrivalTimes,
                         const std::vector<double> &weights);
    /// @result True indicates the arrival times have been set.
    [[nodiscard]] bool haveArrivalTimes() const noexcept;

    /// @param[in] travelTimes  The travel time, in seconds, corresponding to
    ///                         each arrival.
    void setTravelTimes(const std::vector<double> &travelTimes);
    /// @result True indicates the travel times have been set.
    [[nodiscard]] bool haveTravelTimes() const noexcept;

    /// @brief Optimizes the origin time for the given times.
    /// @param[in] arrivalTimes  The observed arrival times in UTC seconds since
    ///                          the epoch.
    /// @param[in] travelTimes   The source-to-receiver travel time in seconds.
    /// @param[in] weights       The weights associated with each observation.
    /// @throws std::invalid_argument if lengths of the arrays are inconsistent
    ///         or any of the weights are negative or there are an insufficient
    ///         number of observations. 
    void optimize();
   
    /// @result The origin time in seconds since the epoch.  This constant
    ///         should be added to the estimates to get the estimated
    ///         arrival times.
    [[nodiscard]] double getTime() const; 
    /// @result True indicates the origin time was computed.
    [[nodiscard]] bool haveTime() const noexcept;

    /// @name Destructors
    /// @{

    /// @brief Resets the class.
    void clear() noexcept;
    /// @brief Destructor.
    ~OriginTime();
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment.
    /// @param[in] originTime  The originTime class to copy to this.
    /// @result A deep copy of originTime.
    OriginTime& operator=(const OriginTime &originTime);
    /// @brief Move assignment.
    /// @param[in,out] originTime  The originTime class whose memory will be
    ///                            moved to this.  On exit, originTime's 
    ///                            behavior is undefined.
    /// @result The memory from originTime moved to this.
    OriginTime& operator=(OriginTime &&originTime) noexcept;
    /// @}
private:
    class OriginTimeImpl;
    std::unique_ptr<OriginTimeImpl> pImpl;
};
}
#endif
