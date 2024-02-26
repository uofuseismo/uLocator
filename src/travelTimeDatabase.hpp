#ifndef UTAH_TRAVEL_TIME_DATABASE_HPP
#define UTAH_TRAVEL_TIME_DATABASE_HPP
#include <vector>
#include <mutex>
#include <tuple>
#include <limits>
#include <string>
#include <map>
#include <filesystem>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#ifndef NDEBUG
#include <cassert>
#endif
#include "utahQuarries.hpp"
#include "originTime.hpp"
#include "optimizers/objectiveFunctions.hpp"
#include "uLocator/corrections/static.hpp"
#include "uLocator/corrections/sourceSpecific.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/uussRayTracer.hpp"
#include "uLocator/travelTimeCalculator.hpp"
#include "uLocator/station.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/optimizers/optimizer.hpp"

namespace
{

/// @brief A known location is simply an (x,y,z) location in the local coordinates.
class IKnownLocation
{
public:
    /// @result The local x coordinates in meters.
    [[nodiscard]] virtual double x() const = 0;
    /// @result The local y coordinates in meters.
    [[nodiscard]] virtual double y() const = 0;
    /// @result The local z coordinates in meters.  This increases + down.
    [[nodiscard]] virtual double z() const = 0;
    virtual ~IKnownLocation() = default;
};

/// @brief A known quarry is the (x,y,z) of a Utah quarry in the Utah region.
class KnownUtahQuarry : public ::IKnownLocation
{
public:
    explicit KnownUtahQuarry(const ULocator::Position::UtahQuarry &quarry) :
        mQuarry(quarry)
    {
    }
    [[nodiscard]] double x() const override
    {
        return mQuarry.getLocalCoordinates().first;
    }
    [[nodiscard]] double y() const override
    {
        return mQuarry.getLocalCoordinates().second;
    }
    [[nodiscard]] double z() const override
    {
        return -mQuarry.getElevation();
    }
    ~KnownUtahQuarry() override = default; 
    ULocator::Position::UtahQuarry mQuarry;
};

/// @class UtahPoint "utahPoint.hpp" "uLocator/position/utahPoint.hpp"
/// @brief This is a geographic point in the Utah region.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license. 
class UtahPoint : public ULocator::Position::IGeographicPoint
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    UtahPoint() :
        ULocator::Position::IGeographicPoint(ULocator::Position::UtahRegion {})
    {
    }
    /// @brief Copy constructor.
    /// @param[in] point   The point from which to initialize this point.
    UtahPoint(const UtahPoint &point) :
        ULocator::Position::IGeographicPoint(ULocator::Position::UtahRegion {})
    {
        *this = point;
    }
    /// @brief Move constructor.
    /// @param[in,out] point  The point from which to initialize this point.
    ///                       On exit, point's behavior is undefined.
    UtahPoint(UtahPoint &&point) noexcept :
        ULocator::Position::IGeographicPoint(ULocator::Position::UtahRegion {})
    {
        *this = std::move(point);
    }
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] point   The point to copy to this.
    /// @result A deep copy of the input point.
    UtahPoint& operator=(const UtahPoint &point)
    {
        if (&point == this){return *this;}
        if (point.haveGeographicCoordinates())
        {
            auto [latitude, longitude] = point.getGeographicCoordinates();
            setGeographicCoordinates(latitude, longitude);
        }
        return *this;
    } 
    /// @brief Move assignment operator.
    /// @param[in,out] point  The point whose memory will be moved to this.
    ///                       On exit, point's behavior is undefined.
    /// @result The memory from point moved to this.
    UtahPoint& operator=(UtahPoint &&point) noexcept
    {   
        if (&point == this){return *this;}
        if (point.haveGeographicCoordinates())
        {
            auto [latitude, longitude] = point.getGeographicCoordinates();
            setGeographicCoordinates(latitude, longitude);
        }
        return *this;
    }
    /// @}

    /// @result A copy of the point.
    [[nodiscard]] std::unique_ptr<ULocator::Position::IGeographicPoint> clone() const override
    {
        std::unique_ptr<ULocator::Position::IGeographicPoint> result
            = std::make_unique<::UtahPoint> (*this);
        return result;
    };
    /// @brief Destructor.
    ~UtahPoint() override = default;
};

/// @class YNPPoint "ynpPoint.hpp" "uLocator/position/ynpPoint.hpp"
/// @brief This is a geographic point in the Yellowstone region.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license. 
class YNPPoint : public ULocator::Position::IGeographicPoint
{
public:
    /// @name Constructors
    /// @{

    /// @brief Constructor.
    YNPPoint() :
        ULocator::Position::IGeographicPoint(ULocator::Position::YNPRegion {}) 
    {   
    }
    /// @brief Copy constructor.
    /// @param[in] point   The point from which to initialize this point.
    YNPPoint(const YNPPoint &point) :
        ULocator::Position::IGeographicPoint(ULocator::Position::YNPRegion {}) 
    {
        *this = point;
    }
    /// @brief Move constructor.
    /// @param[in,out] point  The point from which to initialize this point.
    ///                       On exit, point's behavior is undefined.
    YNPPoint(YNPPoint &&point) noexcept :
        ULocator::Position::IGeographicPoint(ULocator::Position::YNPRegion {}) 
    {
        *this = std::move(point);
    } 
    /// @}

    /// @name Operators
    /// @{

    /// @brief Copy assignment operator.
    /// @param[in] point   The point to copy to this.
    /// @result A deep copy of the input point.
    YNPPoint& operator=(const YNPPoint &point)
    {
        if (&point == this){return *this;}
        if (point.haveGeographicCoordinates())
        {
            auto [latitude, longitude] = point.getGeographicCoordinates();
            setGeographicCoordinates(latitude, longitude);
        }
        return *this;
    }   
    /// @brief Move assignment operator.
    /// @param[in,out] point  The point whose memory will be moved to this.
    ///                       On exit, point's behavior is undefined.
    /// @result The memory from point moved to this.
    YNPPoint& operator=(YNPPoint &&point) noexcept
    {
        if (&point == this){return *this;}
        if (point.haveGeographicCoordinates())
        {
            auto [latitude, longitude] = point.getGeographicCoordinates();
            setGeographicCoordinates(latitude, longitude);
        }
        return *this;
    }   
    /// @}

    /// @result A copy of the point.
    [[nodiscard]] std::unique_ptr<ULocator::Position::IGeographicPoint> clone() const override
    {
        std::unique_ptr<ULocator::Position::IGeographicPoint> result
            = std::make_unique<::YNPPoint> (*this);
        return result;
    };  
    /// @brief Destructor.
    ~YNPPoint() override = default;
};

/// @brief A known Utah event is the (x,y,z) of an event in the Utah region. 
class KnownUtahEvent : public ::IKnownLocation
{
public:
    KnownUtahEvent(const KnownUtahEvent &event)
    {
        *this = event;
    }
    KnownUtahEvent(KnownUtahEvent &&event) noexcept
    {
        *this = std::move(event);
    }
    KnownUtahEvent(const double latitude, const double longitude, const double depth)
    {
        mPoint->setGeographicCoordinates(latitude, longitude);
        mDepth = depth;
    }
    [[nodiscard]] double x() const override
    {
        return mPoint->getLocalCoordinates().first;
    }  
    [[nodiscard]] double y() const override
    {
        return mPoint->getLocalCoordinates().second;
    }
    [[nodiscard]] double z() const override
    {
        return mDepth;
    }
    KnownUtahEvent &operator=(const KnownUtahEvent &event)
    {
        if (&event == this){return *this;}
        mPoint = std::make_unique<::UtahPoint> (*event.mPoint);
        mDepth = event.mDepth;
        return *this;
    }
    KnownUtahEvent &operator=(KnownUtahEvent &&event) noexcept
    {
        if (&event == this){return *this;}
        mPoint = std::move(event.mPoint);
        mDepth = event.mDepth;
        return *this;
    }
    [[nodiscard]] std::pair<double, double> getGeographicCoordinates() const
    {
        return mPoint->getGeographicCoordinates();
    } 
    ~KnownUtahEvent() override = default;
private:
    std::unique_ptr<::UtahPoint> mPoint{std::make_unique<::UtahPoint> ()};
    double mDepth{0};
};

/// @brief A known YNP event is the (x,y,z) of an event in the Yellowstone region. 
class KnownYNPEvent : public ::IKnownLocation
{
public:
    KnownYNPEvent(const KnownYNPEvent &event)
    {
        *this = event;
    }
    KnownYNPEvent(KnownYNPEvent &&event) noexcept
    {
        *this = std::move(event);
    }
    KnownYNPEvent(const double latitude, const double longitude, const double depth)
    {
        mPoint->setGeographicCoordinates(latitude, longitude);
        mDepth = depth;
    }
    [[nodiscard]] double x() const override
    {
        return mPoint->getLocalCoordinates().first;
    }  
    [[nodiscard]] double y() const override
    {
        return mPoint->getLocalCoordinates().second;
    }
    [[nodiscard]] double z() const override
    {
        return mDepth;
    }
    [[nodiscard]] std::pair<double, double> getGeographicCoordinates() const
    {
        return mPoint->getGeographicCoordinates();
    }
    KnownYNPEvent &operator=(const KnownYNPEvent &event)
    {   
        if (&event == this){return *this;}
        mPoint = std::make_unique<::YNPPoint> (*event.mPoint);
        mDepth = event.mDepth;
        return *this;
    }   
    KnownYNPEvent &operator=(KnownYNPEvent &&event) noexcept
    {   
        if (&event == this){return *this;}
        mPoint = std::move(event.mPoint);
        mDepth = event.mDepth;
        return *this;
    }
    ~KnownYNPEvent() override = default;
private:
    std::unique_ptr<::YNPPoint> mPoint{std::make_unique<::YNPPoint> ()};
    double mDepth{0};
};

std::string stationPhaseToString(
    const ULocator::Station &station,
    const std::string &phase)
{
    return station.getNetwork() + "." + station.getName() + "." + phase;
}

std::string quarryStationPhaseToString(
    const ULocator::Position::UtahQuarry &quarry,
    const ULocator::Station &station,
    const std::string &phase)
{
    auto stationName = station.getNetwork() + "." + station.getName();
    return quarry.getName() + "." + ::stationPhaseToString(station, phase);
}

/// @class TravelTimeDatabase
/// @brief A travel time database stores travel times from a known location to
///        a station, phase pair.  This particularly instance is dynamic - i.e,.
///        station, phase pairs can be added as they are observed.
/// @copyright Ben Baker (University of Utah) distributed under the MIT license.
class TravelTimeDatabase
{
//public:
//    using VectorType = std::vector<std::vector<double>>;
public:
    /// @brief Defines a precomputed travel time.
    struct PrecomputedTime
    {
        double travelTimeUncorrected{0}; /*!< The uncorrected travel time in seconds. */
        double dtdt0Uncorrected{0};  /*!< The derivative of the uncorrected travel time with respect to the origin.  This is unitless. */
        double dtdxUncorrected{0};   /*!< The derivative of the uncorrected travel time with respect to x.  This has units of s/m. */
        double dtdyUncorrected{0};   /*!< The derivative of the uncorrected travel time with respect to y.  This has units of s/m. */
        double dtdzUncorrected{0};   /*!< The derivative of the uncorrected travel time with respect to z.  This has units of s/m. */
        double travelTime{0}; /*!< The corrected travel time in seconds. */
        double dtdt0{0};  /*!< The derivative of the corrected travel time with respect to the origin.  This is unitless. */
        double dtdx{0};   /*!< The derivative of the corrected travel time with respect to x.  This has units of s/m. */
        double dtdy{0};   /*!< The derivative of the corrected travel time with respect to y.  This has units of s/m. */
        double dtdz{0};   /*!< The derivative of the corrected travel time with respect to z.  This has units of s/m. */
    };
public:
    TravelTimeDatabase(ULocator::UUSSRayTracer::Region regionIdentifier,
                       std::unique_ptr<ULocator::Position::IGeographicRegion> &&region) :
        mRegion(std::move(region)),
        mRegionIdentifier(regionIdentifier),
        mLogger(std::make_shared<UMPS::Logging::StandardOut> ())
    {
    }
    /// @brief Sets the known locations.
    /// @param[in] knownLocations  The known locations on which to build 
    ///                            the database.
    void setKnownLocations(
        std::vector<std::unique_ptr<::IKnownLocation>> &&knownLocations)
    {
        for (const auto &knownLocation : knownLocations)
        {
            if (knownLocation == nullptr)
            {
                throw std::invalid_argument("Null known location");
            }
        }
        mKnownLocations = std::move(knownLocations);
    }
    /// @brief Adds a station/phase pair to the database.
    void addStationPhasePair(
        const std::pair<ULocator::Station, std::string> &stationPhase)
    {
        _addStationPhasePair(stationPhase);
    }
    /// @brief Evaluates the travel times from the known locations to the
    ///        arrivals attached to this origin.
    /// @param[in] origin  The origin with the arrivals (station, phase) pairs
    ///                    for which to compute travel times.
    /// @param[in] applyCorrection  True indicates the travel times corrections
    ///                             should be supplied.
    [[nodiscard]] std::vector<std::vector<double>>
        evaluate(const ULocator::Origin &origin,
                 const bool applyCorrection = true) const
    {
        const auto &arrivals = origin.getArrivalsReference(); 
        return evaluate(arrivals, applyCorrection);
    }
    /// @brief Evaluates the travel times from the known locations to the
    ///        arrivals.  If an arrival is not in the database then it will
    ///        be added.
    /// @param[in] arrivals         The arrivals (station, phase) for which to
    ///                             compute travel times.
    /// @param[in] applyCorrection  True indicates the travel times corrections
    ///                             should be supplied.
    [[nodiscard]] std::vector<std::vector<double>>
         evaluate(const std::vector<ULocator::Arrival> &arrivals,
                  const bool applyCorrection = true) const
    {
        // Convert the arrivals to station/phase pairs and see if we need
        // to do any updates for new stations
        std::vector<std::pair<ULocator::Station, std::string>> stationPhases;
        for (const auto &arrival : arrivals)
        {
            stationPhases.push_back( std::pair {arrival.getStation(),
                                                arrival.getPhase()} );
            auto stationPhase
                = ::stationPhaseToString(stationPhases.back().first,
                                         stationPhases.back().second);
            if (!mPhaseTravelTimeDatabase.contains(stationPhase))
            {
                std::lock_guard<std::mutex> lockGuard(mMutex);
                {
                _addStationPhasePair(stationPhases.back());
                }
            }
        }
        // Allocate space
        auto nLocations = getNumberOfLocations();
        std::vector<std::vector<double>> estimateTimes;
        estimateTimes.resize(nLocations);
        for (int iLocation = 0; iLocation < nLocations; ++iLocation)
        {
            estimateTimes[iLocation].resize(arrivals.size(), 0);
        }
        // Tabulate the arrival times
        for (int iArrival = 0;
             iArrival < static_cast<int> (stationPhases.size());
             ++iArrival)
        {
            auto key = ::stationPhaseToString(stationPhases[iArrival].first,
                                              stationPhases[iArrival].second);
            const auto &precomputedTravelTimes
                = mPhaseTravelTimeDatabase[key]; //.second;
#ifndef NDEBUG
            assert(precomputedTravelTimes.size() == mKnownLocations.size());
#endif
            if (applyCorrection)
            {
                for (int iLocation = 0; iLocation < nLocations; ++iLocation)
                {
                    estimateTimes[iLocation][iArrival]
                        = precomputedTravelTimes[iLocation].travelTime;
                }
            }
            else
            {
                for (int iLocation = 0; iLocation < nLocations; ++iLocation)
                {
                    estimateTimes[iLocation][iArrival]
                        = precomputedTravelTimes[iLocation]
                         .travelTimeUncorrected;
                }
            }
        }
        return estimateTimes;
    }
    /// @brief This finds the known location that best explains the arrivals
    ///        contains on this origin.  This works by computing an optimal
    ///        origin time for each location then scoring the corresponding
    ///        objective function.
    /// @param[in] origin   The origin with the arrivals.
    /// @param[in] norm     The norm in which to optimize the origin time. 
    /// @param[in] p        The p norm if using an Lp norm.
    /// @param[in] applyCorrection  True indicates the corrected travel
    ///                             times should be used.
    /// @result [best known location index,
    ///          corresponding origin,
    ///          corresponding objective function]
    ///         of the best-fitting location in the database.
    [[nodiscard]] std::tuple<int, ULocator::Origin, double>
        findBestLocation(const ULocator::Origin &origin,
                         const ULocator::Optimizers::IOptimizer::Norm norm,
                         const double p = 1.5,
                         const double timeWindow = 150,
                         const bool applyCorrection = true)
    {
        // Extract arrivals
        auto arrivals = origin.getArrivals();
        auto nArrivals = static_cast<int> (arrivals.size());
        if (nArrivals < 2)
        {
            throw std::runtime_error("At least two arrivals required");
        }
        // Compute travel times
        auto estimateTimes = evaluate(origin, applyCorrection); 
        // Find the quarry with the smallest loss
        auto nLocations = getNumberOfLocations();
        std::vector<std::pair<double, double>> objectiveFunctions(nLocations);
        for (int i = 0; i < nLocations; ++i)
        {
            auto travelTimes = estimateTimes[i];
            double originTime;
            if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
            {
                originTime
                    = ::optimizeOriginTimeLeastSquares(origin, travelTimes);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
            {
                originTime
                    = ::optimizeOriginTimeL1(origin, travelTimes);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
            {
                originTime = ::optimizeOriginTimeLp(origin, travelTimes,
                                                    p, timeWindow);
            } 
#ifndef NDEBUG
            else
            {
                assert(false);
            }
#endif
            std::transform(travelTimes.begin(), travelTimes.end(),
                           travelTimes.begin(), [&](const auto &value)
                           {
                              return value + originTime;
                           });
            if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
            {
                objectiveFunctions[i].first
                    = ::leastSquares(origin, travelTimes,
                                     ::Measurement::Standard);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
            {
                objectiveFunctions[i].first
                    = ::l1(origin, travelTimes, ::Measurement::Standard);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
            {
                objectiveFunctions[i].first
                    = ::lp(origin, travelTimes, p, ::Measurement::Standard);
            }
            objectiveFunctions[i].second = originTime;
        }
        // Find the smallest objective function
        auto it = std::min_element(objectiveFunctions.begin(),
                                   objectiveFunctions.end(),
                                   [](const auto &lhs, const auto &rhs)
                                   {
                                      return lhs.first < rhs.first;
                                   });
        int bestIndex
             = static_cast<int> (std::distance(objectiveFunctions.begin(), it));
        auto bestObjectiveFunction = objectiveFunctions.at(bestIndex).first;
        auto originTime = objectiveFunctions.at(bestIndex).second;
        // Build the origin
        ULocator::Origin bestOrigin{origin}; 
        bestOrigin.setTime(originTime);
        auto [latitude, longitude] 
           = mRegion->localToGeographicCoordinates(
                mKnownLocations.at(bestIndex)->x(),
                mKnownLocations.at(bestIndex)->y());
        bestOrigin.setEpicenter(
             ULocator::Position::WGS84 {latitude, longitude});
        bestOrigin.setDepth(mKnownLocations.at(bestIndex)->z());
        // Tabulate the residuals
        for (int i = 0; i < nArrivals; ++i)
        {
            auto observedTime = arrivals.at(i).getTime();
            auto estimateTime = estimateTimes[bestIndex].at(i) + originTime;
            arrivals[i].setResidual(observedTime - estimateTime);
        }
        bestOrigin.setArrivals(arrivals);
        return std::tuple {bestIndex, bestOrigin, bestObjectiveFunction};
    }
    /// @brief Sets the source-specific corrections HDF5 file
    void setSourceSpecificCorrectionsFile(const std::filesystem::path &file)
    {
        if (!std::filesystem::exists(file))
        {
            throw std::invalid_argument("SSSC file " + file.string()
                                      + " does not exist");
        }
        mSourceSpecificCorrectionsFile = file;
    }
    /// @brief Sets the static corrections HDF5 file
    void setStaticCorrectionsFile(const std::filesystem::path &file)
    {
        if (!std::filesystem::exists(file))
        {
            throw std::invalid_argument("Static corrections file "
                                      + file.string() + " does not exist");
        }
        mStaticCorrectionsFile = file;
    } 
    /// @result The number of known locations.
    int getNumberOfLocations() const noexcept
    {
        return static_cast<int> (mKnownLocations.size());
    }
    /// @result The number of station/phase pairs in the database.
    /// @note The memory usage is
    ///       ~ getNumberOfStationPhasePairs() x getNumberOfLocations() x 8.
    int getNumberOfStationPhasePairs() const noexcept
    {
        return static_cast<int> (mPhaseTravelTimeDatabase.size());
    }
/*
    /// @result The travel times from the i'th known location to all arrivals.
    const std::vector<double> &operator[](const size_t i) const
    {
        return mEstimateTimes[i];
    }
    /// @result The travel times from the i'th known location to all arrivals. 
    const std::vector<double> &operator()(const size_t i) const
    {
        return mEstimateTimes.at(i);
    } 
*/
/*
    VectorType::iterator begin() {return mEstimateTimes.begin();}
    VectorType::iterator end()   {return mEstimateTimes.end();}
    VectorType::const_iterator cbegin() const {return mEstimateTimes.cbegin();}
    VectorType::const_iterator cend() const   {return mEstimateTimes.cend();}
*/
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::unique_ptr<ULocator::Position::IGeographicRegion> mRegion{nullptr};
    mutable std::map<std::string, std::vector<PrecomputedTime>>
        mPhaseTravelTimeDatabase;
    mutable std::mutex mMutex;
    std::vector<std::unique_ptr<::IKnownLocation>> mKnownLocations;
    mutable std::vector<double> mBestResiduals;
    ULocator::UUSSRayTracer::Region mRegionIdentifier;
    std::filesystem::path mStaticCorrectionsFile;
    std::filesystem::path mSourceSpecificCorrectionsFile;
private:
    /// @brief Adds a station/phase pair to the database.
    void _addStationPhasePair(
        const std::pair<ULocator::Station, std::string> &stationPhase) const
    {
        constexpr double originTime{0};
        const auto &station = stationPhase.first;
        auto [xStation, yStation] = station.getLocalCoordinates();
        auto phase = stationPhase.second;
        auto key = ::stationPhaseToString(station, phase);
        if (phase == "P")
        {
            ULocator::Corrections::Static staticCorrection;
            staticCorrection.setStationNameAndPhase(station.getNetwork(),
                                                    station.getName(),
                                                    "P");
            ULocator::Corrections::SourceSpecific sourceSpecific;
            sourceSpecific.setStationNameAndPhase(station.getNetwork(),
                                                  station.getName(),
                                                  "P");
            if (std::filesystem::exists(mStaticCorrectionsFile))
            {
                try
                {
                    staticCorrection.load(mStaticCorrectionsFile);
                }
                catch (const std::exception &e)
                {
                    mLogger->debug(e.what());
                }
            }
            if (std::filesystem::exists(mSourceSpecificCorrectionsFile))
            {
                try
                {
                    sourceSpecific.load(mSourceSpecificCorrectionsFile);
                }
                catch (const std::exception &e)
                {
                    mLogger->debug(e.what());
                }
            }
            ULocator::UUSSRayTracer 
                pRayTracer(station,
                           ULocator::UUSSRayTracer::Phase::P,
                           mRegionIdentifier,
                           std::move(staticCorrection),
                           std::move(sourceSpecific),
                           nullptr);
            std::vector<PrecomputedTime> travelTimes;
            travelTimes.reserve(mKnownLocations.size());
            for (const auto &knownLocation : mKnownLocations)
            {
                auto xSource = knownLocation->x();
                auto ySource = knownLocation->y();
                auto zSource = knownLocation->z();
                PrecomputedTime precomputedTime;
                precomputedTime.travelTimeUncorrected
                     = pRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0Uncorrected,
                                           &precomputedTime.dtdxUncorrected,
                                           &precomputedTime.dtdyUncorrected,
                                           &precomputedTime.dtdzUncorrected,
                                           false);
                precomputedTime.travelTime
                     = pRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0,
                                           &precomputedTime.dtdx,
                                           &precomputedTime.dtdy,
                                           &precomputedTime.dtdz, 
                                           true);
                //std::cout << precomputedTime.travelTime << " " << precomputedTime.travelTimeUncorrected << std::endl;
                travelTimes.push_back(std::move(precomputedTime));
            }
            mPhaseTravelTimeDatabase.insert(
                std::pair {key, travelTimes} );
        }        
        else if (phase == "S")
        {
            ULocator::Corrections::Static staticCorrection;
            staticCorrection.setStationNameAndPhase(station.getNetwork(),
                                                    station.getName(),
                                                    "S");
            ULocator::Corrections::SourceSpecific sourceSpecific;
            sourceSpecific.setStationNameAndPhase(station.getNetwork(),
                                                  station.getName(),
                                                  "S");
            if (std::filesystem::exists(mStaticCorrectionsFile))
            {
                try
                {
                    staticCorrection.load(mStaticCorrectionsFile);
                }
                catch (const std::exception &e)
                {
                    mLogger->debug(e.what());
                }
            }
            if (std::filesystem::exists(mSourceSpecificCorrectionsFile))
            {
                try
                {
                    sourceSpecific.load(mSourceSpecificCorrectionsFile);
                }
                catch (const std::exception &e)
                {
                    mLogger->debug(e.what());
                }
            }
            ULocator::UUSSRayTracer 
                sRayTracer(station,
                           ULocator::UUSSRayTracer::Phase::S,
                           mRegionIdentifier,
                           std::move(staticCorrection),
                           std::move(sourceSpecific),
                           nullptr);
            std::vector<PrecomputedTime> travelTimes;
            travelTimes.reserve(mKnownLocations.size());
            for (const auto &knownLocation : mKnownLocations)
            {
                auto xSource = knownLocation->x();
                auto ySource = knownLocation->y();
                auto zSource = knownLocation->z();
                PrecomputedTime precomputedTime;
                precomputedTime.travelTimeUncorrected
                     = sRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0Uncorrected,
                                           &precomputedTime.dtdxUncorrected,
                                           &precomputedTime.dtdyUncorrected,
                                           &precomputedTime.dtdzUncorrected,
                                           false);
                precomputedTime.travelTime
                     = sRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0,
                                           &precomputedTime.dtdx,
                                           &precomputedTime.dtdy,
                                           &precomputedTime.dtdz,
                                           true);
                travelTimes.push_back(std::move(precomputedTime));
            }
            mPhaseTravelTimeDatabase.insert(
                std::pair {key, travelTimes} );
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#else
        else
        {
            throw std::invalid_argument("Unhandled phase: " + phase);
        }
#endif
    }
};

class UtahQuarryBlastTravelTimeDatabase : public TravelTimeDatabase
{
public:
    UtahQuarryBlastTravelTimeDatabase() :
        TravelTimeDatabase(ULocator::UUSSRayTracer::Region::Utah,
                           ULocator::Position::UtahRegion {}.clone())
    {
        auto nQuarries = getNumberOfQuarries();
        std::vector<std::unique_ptr<::IKnownLocation>> knownLocations;
        for (int i = 0; i < nQuarries; ++i)
        {
            knownLocations.push_back(std::make_unique<::KnownUtahQuarry> (mQuarries[i]));
        }
        setKnownLocations(std::move(knownLocations));
    }
    int getNumberOfQuarries() const noexcept
    {
        return static_cast<int> (mQuarries.size());
    }
    [[nodiscard]] std::tuple<int, ULocator::Origin, double>
        findBestQuarry(const ULocator::Origin &origin,
                       const ULocator::Optimizers::IOptimizer::Norm norm,
                       const double p = 1.5,
                       const double timeWindow = 150,
                       const bool applyCorrection = true)
    {
        return findBestLocation(origin, norm, p, timeWindow, applyCorrection);
    }
    ~UtahQuarryBlastTravelTimeDatabase() = default;
    [[nodiscard]] ULocator::Position::UtahQuarry operator[](const size_t index) const
    {    
        return mQuarries[index];
    }    
    [[nodiscard]] ULocator::Position::UtahQuarry at(const size_t index) const
    {
        return mQuarries.at(index);
    }
private:
    std::vector<ULocator::Position::UtahQuarry> mQuarries{::getUtahQuarries()};
};

#include "knownEvents.hpp"
class UtahEventTravelTimeDatabase : public TravelTimeDatabase
{
public:
    UtahEventTravelTimeDatabase() :
        TravelTimeDatabase(ULocator::UUSSRayTracer::Region::Utah,
                           ULocator::Position::UtahRegion {}.clone())
    {
        auto nEvents = getNumberOfEvents();
        std::vector<std::unique_ptr<::IKnownLocation>> knownLocations;
        for (int i = 0; i < nEvents; ++i)
        {
            knownLocations.push_back(std::make_unique<::KnownUtahEvent> (mEvents[i]));
        }
        setKnownLocations(std::move(knownLocations));
    }   
    int getNumberOfEvents() const noexcept
    {
        return static_cast<int> (mEvents.size());
    }   
    [[nodiscard]] std::tuple<int, ULocator::Origin, double>
        findBestEvent(const ULocator::Origin &origin,
                      const ULocator::Optimizers::IOptimizer::Norm norm,
                      const double p = 1.5,
                      const double timeWindow = 150,
                      const bool applyCorrection = true)
    {
        return findBestLocation(origin, norm, p, timeWindow, applyCorrection);
    }
    ~UtahEventTravelTimeDatabase() = default;
    [[nodiscard]] ::KnownUtahEvent operator[](const size_t index) const
    {
        return mEvents[index];
    }
    [[nodiscard]] ::KnownUtahEvent operator()(const size_t index) const
    {
        return mEvents.at(index);
    }
private:
    std::vector<::KnownUtahEvent> mEvents{::utahKnownEvents()};
};

class YNPEventTravelTimeDatabase : public TravelTimeDatabase
{
public:
    YNPEventTravelTimeDatabase() :
        TravelTimeDatabase(ULocator::UUSSRayTracer::Region::YNP,
                           ULocator::Position::YNPRegion {}.clone())
    {   
        auto nEvents = getNumberOfEvents();
        std::vector<std::unique_ptr<::IKnownLocation>> knownLocations;
        for (int i = 0; i < nEvents; ++i)
        {
            knownLocations.push_back(std::make_unique<::KnownYNPEvent> (mEvents[i]));
        }
        setKnownLocations(std::move(knownLocations));
    }   
    int getNumberOfEvents() const noexcept
    {   
        return static_cast<int> (mEvents.size());
    }   
    [[nodiscard]] std::tuple<int, ULocator::Origin, double>
        findBestEvent(const ULocator::Origin &origin,
                      const ULocator::Optimizers::IOptimizer::Norm norm,
                      const double p = 1.5,
                      const double timeWindow = 150,
                      const bool applyCorrection = true)
    {   
        return findBestLocation(origin, norm, p, timeWindow, applyCorrection);
    }   
    ~YNPEventTravelTimeDatabase() = default;
    [[nodiscard]] ::KnownYNPEvent operator[](const size_t index) const
    {    
        return mEvents[index];
    }    
    [[nodiscard]] ::KnownYNPEvent operator()(const size_t index) const
    {
        return mEvents.at(index);
    }
private:
    std::vector<::KnownYNPEvent> mEvents{::ynpKnownEvents()};
};

/*
class UtahQuarryBlastTravelTimeDatabase
{
public:
    using VectorType = std::vector<std::vector<double>>;
public:
    struct PrecomputedTime
    {
        double travelTimeUncorrected{0}; 
        double dtdt0Uncorrected{0};
        double dtdxUncorrected{0};
        double dtdyUncorrected{0};
        double dtdzUncorrected{0};
        double travelTime{0};
        double dtdt0{0};
        double dtdx{0};
        double dtdy{0};
        double dtdz{0};
    };
public:
    UtahQuarryBlastTravelTimeDatabase()
    {
        mQuarries = ::getUtahQuarries(); 
    }
    // Compute stationPhase pair
    void addStationPhasePair(
        const std::pair<ULocator::Station, std::string> &stationPhase)
    {
        constexpr double originTime{0};
        const auto &station = stationPhase.first;
        auto [xStation, yStation] = station.getLocalCoordinates();
        auto phase = stationPhase.second;
        auto key = ::stationPhaseToString(station, phase);
        if (phase == "P")
        {
            ULocator::UUSSRayTracer 
                pRayTracer(station,
                           ULocator::UUSSRayTracer::Phase::P,
                           ULocator::UUSSRayTracer::Region::Utah,
                           nullptr);
            std::vector<PrecomputedTime> travelTimes;
            travelTimes.reserve(mQuarries.size());
            for (const auto &quarry : mQuarries)
            {
                auto [xSource, ySource] =  quarry.getLocalCoordinates();
                auto zSource =-quarry.getElevation();
                PrecomputedTime precomputedTime;
                precomputedTime.travelTimeUncorrected
                     = pRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0Uncorrected,
                                           &precomputedTime.dtdxUncorrected,
                                           &precomputedTime.dtdyUncorrected,
                                           &precomputedTime.dtdzUncorrected,
                                           false);
                precomputedTime.travelTime
                     = pRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0,
                                           &precomputedTime.dtdx,
                                           &precomputedTime.dtdy,
                                           &precomputedTime.dtdz, 
                                           true);
                travelTimes.push_back(std::move(precomputedTime));
            }
            mQuarryPhaseTravelTimeDatabase.insert(
                std::pair {key, travelTimes} );
        }        
        else if (phase == "S")
        {
            ULocator::UUSSRayTracer 
                sRayTracer(station,
                           ULocator::UUSSRayTracer::Phase::S,
                           ULocator::UUSSRayTracer::Region::Utah,
                           nullptr);
            std::vector<PrecomputedTime> travelTimes;
            travelTimes.reserve(mQuarries.size());
            for (const auto &quarry : mQuarries)
            {
                auto [xSource, ySource] =  quarry.getLocalCoordinates();
                auto zSource =-quarry.getElevation();
                PrecomputedTime precomputedTime;
                precomputedTime.travelTimeUncorrected
                     = sRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0Uncorrected,
                                           &precomputedTime.dtdxUncorrected,
                                           &precomputedTime.dtdyUncorrected,
                                           &precomputedTime.dtdzUncorrected,
                                           false);
                precomputedTime.travelTime
                     = sRayTracer.evaluate(originTime,
                                           xSource, ySource, zSource,
                                           &precomputedTime.dtdt0,
                                           &precomputedTime.dtdx,
                                           &precomputedTime.dtdy,
                                           &precomputedTime.dtdz,
                                           true);
                travelTimes.push_back(std::move(precomputedTime));
            }
            mQuarryPhaseTravelTimeDatabase.insert(
                std::pair {key, travelTimes} );
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#else
        else
        {
            throw std::invalid_argument("Unhandled phase: " + phase);
        }
#endif
    }
    void evaluate(const ULocator::Origin &origin,
                  const bool applyCorrection = true)
    {
        const auto &arrivals = origin.getArrivalsReference(); 
        evaluate(arrivals, applyCorrection);
    }
    /// @result The result is a [nQuarries x nArrivals] matrix of travel times.
    void evaluate(const std::vector<ULocator::Arrival> &arrivals,
                  const bool applyCorrection = true)
    {
        // Convert the arrivals to station/phase pairs and see if we need
        // to do any updates for new stations
        std::vector<std::pair<ULocator::Station, std::string>> stationPhases;
        for (const auto &arrival : arrivals)
        {
            stationPhases.push_back( std::pair {arrival.getStation(),
                                                arrival.getPhase()} );
            auto stationPhase
                = ::stationPhaseToString(stationPhases.back().first,
                                         stationPhases.back().second);
            if (!mQuarryPhaseTravelTimeDatabase.contains(stationPhase))
            {
                addStationPhasePair(stationPhases.back());
            }
        }
        // Allocate space
        auto nQuarries = static_cast<int> (mQuarries.size());
        mEstimateTimes.resize(nQuarries);
        for (int iQuarry = 0; iQuarry < nQuarries; ++iQuarry)
        {
            mEstimateTimes[iQuarry].resize(arrivals.size(), 0);
        }
        // Tabulate the arrival times
        for (int iArrival = 0;
             iArrival < static_cast<int> (stationPhases.size());
             ++iArrival)
        {
            auto key = ::stationPhaseToString(stationPhases[iArrival].first,
                                              stationPhases[iArrival].second);
            const auto &precomputedTravelTimes
                = mQuarryPhaseTravelTimeDatabase[key]; //.second;
#ifndef NDEBUG
            assert(precomputedTravelTimes.size() == mQuarries.size());
#endif
            if (applyCorrection)
            {
                for (int iQuarry = 0; iQuarry < nQuarries; ++iQuarry)
                {
                    mEstimateTimes[iQuarry][iArrival]
                        = precomputedTravelTimes[iQuarry].travelTime;
                }
            }
        }
    }
    [[nodiscard]] std::tuple<double, ULocator::Position::UtahQuarry, double>
        findBestQuarry(const ULocator::Origin &origin,
                       const ULocator::Optimizers::IOptimizer::Norm norm,
                       const double p = 1.5,
                       const bool applyCorrection = true)
    {
        // Extract arrivals
        const auto &arrivals = origin.getArrivalsReference();
        auto nArrivals = static_cast<int> (arrivals.size());
        if (nArrivals < 2)
        {
            throw std::runtime_error("At least two arrivals required");
        }
        std::vector<double> observations(nArrivals);
        std::vector<double> weights(nArrivals);
        for (int i = 0; i < nArrivals; ++i)
        {
            observations[i] = arrivals[i].getTime();
            weights[i] = 1./arrivals[i].getStandardError();
        }
        // Compute travel times
        evaluate(origin, applyCorrection); 
        // Find the quarry with the smallest loss
        double originTime{0};
        auto nQuarries = getNumberOfQuarries();
        std::vector<double>
            objectiveFunctions(nQuarries, std::numeric_limits<double>::max());
        for (int i = 0; i < nQuarries; ++i)
        {
            auto travelTimes = mEstimateTimes[i];
            if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
            {
                originTime
                    = ::optimizeOriginTimeLeastSquares(origin, travelTimes);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
            {
                originTime
                    = ::optimizeOriginTimeL1(origin, travelTimes);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
            {
                originTime
                    = ::optimizeOriginTimeLp(origin, travelTimes, p);
            } 
#ifndef NDEBUG
            else
            {
                assert(false);
            }
#endif
            std::transform(travelTimes.begin(), travelTimes.end(),
                           travelTimes.begin(), [&](const auto &value)
                           {
                              return value + originTime;
                           });
            if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
            {
                objectiveFunctions[i]
                    = ::leastSquares(origin, travelTimes,
                                     ::Measurement::Standard);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
            {
                objectiveFunctions[i]
                    = ::l1(origin, travelTimes, ::Measurement::Standard);
            }
            else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
            {
                objectiveFunctions[i]
                    = ::lp(origin, travelTimes, p, ::Measurement::Standard);
            }

        }
        // Find the smallest objective function
        auto it = std::min_element(objectiveFunctions.begin(),
                                   objectiveFunctions.end());
        size_t idx = std::distance(objectiveFunctions.begin(), it); 
        if (objectiveFunctions.at(idx) == std::numeric_limits<double>::max())
        {
            throw std::runtime_error("Could not find best quarry");
        }
        return std::tuple {objectiveFunctions.at(idx), mQuarries.at(idx), originTime};
    }
    int getNumberOfQuarries() const noexcept
    {
        return static_cast<int> (mQuarries.size());
    }
    int getNumberOfStationPhasePairs() const noexcept
    {
        return static_cast<int> (mQuarryPhaseTravelTimeDatabase.size());
    }
    const std::vector<double> &operator[](const size_t i) const
    {
        return mEstimateTimes[i];
    }
    const std::vector<double> &operator()(const size_t i) const
    {
        return mEstimateTimes.at(i);
    } 
    VectorType::iterator begin() {return mEstimateTimes.begin();}
    VectorType::iterator end()   {return mEstimateTimes.end();}
    VectorType::const_iterator cbegin() const {return mEstimateTimes.cbegin();}
    VectorType::const_iterator cend() const   {return mEstimateTimes.cend();}

    mutable std::map<std::string, std::vector<PrecomputedTime>>
        mQuarryPhaseTravelTimeDatabase;
    std::vector<ULocator::Position::UtahQuarry> mQuarries;
    std::vector<std::vector<double>> mEstimateTimes;
    ULocator::Position::Utah mGeographicRegion;
    std::filesystem::path mCorrectionsFile;
    bool mApplyCorrection{true};
};
*/
}
#endif
