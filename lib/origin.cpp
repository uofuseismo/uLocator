#include <vector>
#include <cmath>
#include <algorithm>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"

using namespace ULocator;

class Origin::OriginImpl
{
public:
    void updateArrivalDistanceAzimuth()
    {
        if (mEpicenter.havePosition())
        {
            for (auto &arrival : mArrivals)
            {
                const auto &station = arrival.getStationReference();
                double azimuth, distance;
                Position::computeDistanceAzimuth(
                    mEpicenter,
                    station.getGeographicPosition(),
                    nullptr,
                    &distance,
                    &azimuth,
                    nullptr);
                arrival.setDistance(distance);
                arrival.setAzimuth(azimuth);
            }
        }
    }
    void computeAzimuthalGapAndNearestStation()
    {
        mAzimuthalGap =-1;
        mNearestStationDistance =-1;
        if (mArrivals.empty()){return;}
        if (!mEpicenter.havePosition()){return;}
        // Compute derivative products
        std::vector<double> azimuths;
        azimuths.reserve(mArrivals.size() + 1);
        double minimumDistance{std::numeric_limits<double>::max()};
        for (const auto &arrival : mArrivals)
        {
           if (arrival.haveAzimuth())
            {
               azimuths.push_back(arrival.getAzimuth());
            }
            if (arrival.haveDistance())
            {
                minimumDistance = std::min(minimumDistance,
                                           arrival.getDistance());
            }
        }
        // Could happen for no arrivals with distances
        if (minimumDistance < std::numeric_limits<double>::max())
        {
            mNearestStationDistance = minimumDistance;
        }
        auto nAzimuths = static_cast<int> (azimuths.size());
        std::sort(azimuths.begin(), azimuths.begin() + nAzimuths);
        azimuths.push_back(azimuths.at(0)); // Trick that makes this work
        double maximumGap{0};
        for (int i = 0; i < nAzimuths; ++i)
        {
            maximumGap = std::max(azimuths[i + 1] - azimuths[i], maximumGap);
        }
        // Could happen for only one station and multiple phases
        if (maximumGap > 0){mAzimuthalGap = maximumGap;}
    }
    void computeWeightedRMSE()
    {
        mWeightedRMSE =-1;
        if (mArrivals.empty()){return;}
        double numerator{0};
        double denominator{0};
        for (const auto &arrival : mArrivals)
        {
            if (arrival.haveStandardError() &&
                arrival.haveTime() &&
                arrival.haveResidual())
            {
#ifndef NDEBUG
                assert(arrival.getStandardError() > 0);
#endif
                double weight = 1./arrival.getStandardError();
                double residual = arrival.getResidual();
                numerator = numerator + weight*(residual*residual);
                denominator = denominator + weight;
            }
        }
        if (denominator > 0){mWeightedRMSE = std::sqrt(numerator/denominator);}
    }
    std::vector<Arrival> mArrivals;
    Position::WGS84 mEpicenter;
    double mOriginTime{0};
    double mDepth{0};
    double mAzimuthalGap{-1};
    double mNearestStationDistance{-1};
    double mWeightedRMSE{-1};
    int64_t mIdentifier{0};
    EventType mEventType{EventType::Unknown};
    bool mDepthFixed{false};
    bool mHaveEpicenter{false};
    bool mHaveDepth{false};
    bool mHaveOriginTime{false};
};

/// Constructor
Origin::Origin() :
    pImpl(std::make_unique<OriginImpl> ())
{
}

/// Copy constructor
Origin::Origin(const Origin &origin)
{
    *this = origin;
}

/// Move constructor
Origin::Origin(Origin &&origin) noexcept
{
    *this = std::move(origin);
}

/// Copy assignment
Origin& Origin::operator=(const Origin &origin)
{
    if (&origin == this){return *this;}
    pImpl = std::make_unique<OriginImpl> (*origin.pImpl);
    return *this;
}

/// Move assignment
Origin& Origin::operator=(Origin &&origin) noexcept
{
    if (&origin == this){return *this;}
    pImpl = std::move(origin.pImpl);
    return *this;
}

/// Reset class
void Origin::clear() noexcept
{
    pImpl = std::make_unique<OriginImpl> ();
}

/// Destructor
Origin::~Origin() = default;

/// Depth
void Origin::setDepth(const double depth, const bool isFixed)
{
    if (depth < -8900 || depth > 6400000)
    {
        throw std::invalid_argument(
           "Depth must be in range [-8900, 6400000] km");
    }
    pImpl->mDepth = depth;
    pImpl->mDepthFixed = isFixed;
    pImpl->mHaveDepth = true;
}

double Origin::getDepth() const
{
    if (!haveDepth()){throw std::runtime_error("Depth not set");}
    return pImpl->mDepth;
}

bool Origin::haveDepth() const noexcept
{
    return pImpl->mHaveDepth;
}

bool Origin::isFixedDepth() const noexcept
{
    return pImpl->mDepthFixed;
}
    
/// Event identifier
void Origin::setIdentifier(int64_t identifier) noexcept
{
    pImpl->mIdentifier = identifier;
}

int64_t Origin::getIdentifier() const noexcept
{
    return pImpl->mIdentifier;
}

/// Arrivals
void Origin::setArrivals(const std::vector<Arrival> &arrivals)
{
    for (const auto &arrival : arrivals)
    {
        if (!arrival.haveTime())
        {
            throw std::invalid_argument("Arrival time not set");
        }
        if (!arrival.haveStation())
        {
            throw std::invalid_argument("Station not set");
        }
    }
    pImpl->mArrivals = arrivals;
    pImpl->updateArrivalDistanceAzimuth();
    pImpl->computeAzimuthalGapAndNearestStation();
    pImpl->computeWeightedRMSE();
}

std::vector<Arrival> Origin::getArrivals() const
{
    return pImpl->mArrivals;
}

const std::vector<Arrival> &Origin::getArrivalsReference() const
{
    return *&pImpl->mArrivals;
}
 
/// Origin time
void Origin::setTime(const double originTime) noexcept
{
    pImpl->mOriginTime = originTime;
    pImpl->mHaveOriginTime = true;
}

double Origin::getTime() const
{
    if (!haveTime()){throw std::runtime_error("Time not set");}
    return pImpl->mOriginTime;
}

bool Origin::haveTime() const noexcept
{
    return pImpl->mHaveOriginTime;
}

/// Location
void Origin::setEpicenter(const Position::WGS84 &epicenter)
{
    pImpl->mEpicenter = epicenter;
    pImpl->mHaveEpicenter = true;
    pImpl->updateArrivalDistanceAzimuth();
    pImpl->computeAzimuthalGapAndNearestStation();
}

Position::WGS84 Origin::getEpicenter() const
{
    if (!haveEpicenter()){throw std::runtime_error("Epicenter not set");}
    return pImpl->mEpicenter;
}

bool Origin::haveEpicenter() const noexcept
{
    return pImpl->mHaveEpicenter;
}

/// Event type
void Origin::setEventType(const EventType eventType) noexcept
{
    pImpl->mEventType = eventType;
}

Origin::EventType Origin::getEventType() const noexcept
{
    return pImpl->mEventType;
}

/// Azimuthal gap
double Origin::getAzimuthalGap() const
{
    if (!haveAzimuthalGap())
    {
        if (!haveEpicenter())
        {
            throw std::runtime_error("Epicenter not set");
        }
        else if (pImpl->mArrivals.empty())
        {
            throw std::runtime_error("No arrivals set");
        }
        throw std::runtime_error("Nearest station distance not computed");
    }
    return pImpl->mAzimuthalGap;
}

bool Origin::haveAzimuthalGap() const noexcept
{
    return (pImpl->mAzimuthalGap >= 0);
}

/// Distance
double Origin::getNearestStationDistance() const
{
    if (!haveNearestStationDistance())
    {
        if (!haveEpicenter())
        {
            throw std::runtime_error("Epicenter not set");
        }
        else if (pImpl->mArrivals.empty())
        {
            throw std::runtime_error("No arrivals set");
        }
        throw std::runtime_error("Nearest station distance not computed");
    }
    return pImpl->mNearestStationDistance;
}

bool Origin::haveNearestStationDistance() const noexcept
{
    return (pImpl->mNearestStationDistance >= 0);
}

/// RMSE
double Origin::getWeightedRootMeanSquaredError() const
{
    if (!haveWeightedRootMeanSquaredError())
    {
        throw std::runtime_error(
           "WeightedRMSE cannot be computed from current arrivals");
    }
    return pImpl->mWeightedRMSE;
}

bool Origin::haveWeightedRootMeanSquaredError() const noexcept
{
    return (pImpl->mWeightedRMSE >= 0);
}
