#include <iostream>
#include <string>
#include <algorithm>
#include <map>
#ifdef USE_TBB
#include <oneapi/tbb/parallel_for.h>
#endif
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/station.hpp"
#include "uLocator/travelTimeCalculator.hpp"

using namespace ULocator;

namespace
{

bool operator<(const std::pair<Station, std::string> &lhs,
               const std::pair<Station, std::string> &rhs)
{
    auto lhsName = lhs.first.getNetwork() + "." + lhs.first.getName()
                 + "." + lhs.second;
    auto rhsName = rhs.first.getNetwork() + "." + rhs.first.getName()
                 + "." + rhs.second;
    return lhsName < rhsName;
}

struct Compare
{
    template<typename T>
    bool operator()(const T &lhs, const T &rhs) const
    {
        return lhs < rhs;
    } 
};

}

class TravelTimeCalculatorMap::TravelTimeCalculatorMapImpl
{
public:
    std::map<std::pair<Station, std::string>, 
             std::unique_ptr<const ITravelTimeCalculator>,
             ::Compare> mTravelTimeCalculatorsMap;
};

/// Constructor
TravelTimeCalculatorMap::TravelTimeCalculatorMap() :
    pImpl(std::make_unique<TravelTimeCalculatorMapImpl> ())
{
}

/// Move constructor
TravelTimeCalculatorMap::TravelTimeCalculatorMap(
    TravelTimeCalculatorMap &&map) noexcept
{
    *this = std::move(map);
}

/// Reset class
void TravelTimeCalculatorMap::clear() noexcept
{
    pImpl->mTravelTimeCalculatorsMap.clear();
}

/// Destructor
TravelTimeCalculatorMap::~TravelTimeCalculatorMap() = default;

/// Number of tables
int TravelTimeCalculatorMap::size() const noexcept
{
    return static_cast<int> (pImpl->mTravelTimeCalculatorsMap.size());
}

/// Evaluate
void TravelTimeCalculatorMap::evaluate(
    const std::vector<std::pair<Station, std::string>> &stationPhases,
    const double originTime,
    const double xSource, const double ySource, const double zSource,
    std::vector<double> *travelTimes,
    const bool applyCorrection) const
{
    if (travelTimes == nullptr)
    {
        throw std::runtime_error("Travel times is NULL");
    }
    if (travelTimes->size() != stationPhases.size())
    {
        travelTimes->resize(stationPhases.size(), 0);
    }
    if (stationPhases.empty()){return;}
#ifdef USE_TBB
    oneapi::tbb::parallel_for
    (
        oneapi::tbb::blocked_range<size_t> (0, stationPhases.size()),
        [&](const oneapi::tbb::blocked_range<size_t> &range)
        {
             for (size_t i = range.begin(); i != range.end(); ++i)
             {
                 travelTimes->at(i)
                     = this->evaluate(stationPhases[i].first,
                                      stationPhases[i].second,
                                      originTime, xSource, ySource, zSource,
                                      applyCorrection);
             }
        }
    );
#else
    for (size_t i = 0; i < stationPhases.size(); ++i)
    {
        travelTimes->at(i)
            = this->evaluate(stationPhases[i].first,
                             stationPhases[i].second,
                             originTime, xSource, ySource, zSource,
                             applyCorrection);
    }
#endif
}

/// Evaluate
void TravelTimeCalculatorMap::evaluate(
    const std::vector<std::pair<Station, std::string>> &stationPhases,
    const double originTime,
    const double xSource, const double ySource, const double zSource,
    std::vector<double> *travelTimes,
    std::vector<double> *dtdt0,
    std::vector<double> *dtdx,
    std::vector<double> *dtdy,
    std::vector<double> *dtdz,
    const bool applyCorrection) const
{
    if (travelTimes == nullptr)
    {   
        throw std::invalid_argument("Travel times is NULL");
    }   
    if (dtdt0 == nullptr){throw std::invalid_argument("dtdt0 is NULL");}
    if (dtdx == nullptr){throw std::invalid_argument("dtdx is NULL");}
    if (dtdy == nullptr){throw std::invalid_argument("dtdy is NULL");}
    if (dtdz == nullptr){throw std::invalid_argument("dtdz is NULL");}
    if (travelTimes->size() != stationPhases.size())
    {   
        travelTimes->resize(stationPhases.size(), 0); 
    }
    if (dtdt0->size() != stationPhases.size())
    {
        dtdt0->resize(stationPhases.size(), 0);
    }
    if (dtdx->size() != stationPhases.size())
    {
        dtdx->resize(stationPhases.size(), 0);
    } 
    if (dtdy->size() != stationPhases.size()) 
    {
        dtdy->resize(stationPhases.size(), 0);
    }
    if (dtdz->size() != stationPhases.size()) 
    {
        dtdz->resize(stationPhases.size(), 0);
    }
    if (stationPhases.empty()){return;}
#ifdef USE_TBB
    oneapi::tbb::parallel_for
    (   
        oneapi::tbb::blocked_range<size_t> (0, stationPhases.size()),
        [&](const oneapi::tbb::blocked_range<size_t> &range)
        {
             for (size_t i = range.begin(); i != range.end(); ++i)
             {
                 double dtdt0i, dtdxi, dtdyi, dtdzi;
                 travelTimes->at(i)
                      = this->evaluate(stationPhases[i].first,
                                       stationPhases[i].second,
                                       originTime, xSource, ySource, zSource,
                                       &dtdt0i, &dtdxi, &dtdyi, &dtdzi,
                                       applyCorrection);
                 dtdt0->at(i) = dtdt0i;
                 dtdx->at(i) = dtdxi;
                 dtdy->at(i) = dtdyi;
                 dtdz->at(i) = dtdzi;
             }   
        }
    );
#else
    for (size_t i = 0; i < stationPhases.size(); ++i)
    {   
        double dtdt0i, dtdxi, dtdyi, dtdzi;
        travelTimes->at(i)
            = this->evaluate(stationPhases[i].first,
                             stationPhases[i].second,
                             originTime, xSource, ySource, zSource,
                             &dtdt0i, &dtdxi, &dtdyi, &dtdzi,
                             applyCorrection);
        dtdt0->at(i) = dtdt0i;
        dtdx->at(i) = dtdxi;
        dtdy->at(i) = dtdyi;
        dtdz->at(i) = dtdzi;
    }
#endif
}

double TravelTimeCalculatorMap::evaluate(
    const Station &station, const std::string &phase,
    const double originTime,
    const double xSource, const double ySource, const double zSource,
    const bool applyCorrection) const
{
    return at(station, phase)->evaluate(originTime, xSource, ySource, zSource,
                                        applyCorrection);
}

double TravelTimeCalculatorMap::evaluate(
    const Station &station, const std::string &phase,
    const double originTime,
    const double xSource, const double ySource, const double zSource,
    double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    return at(station, phase)->evaluate(originTime, xSource, ySource, zSource,
                                        dtdt0, dtdx, dtdy, dtdz,
                                        applyCorrection);
}

/// Insert
void TravelTimeCalculatorMap::insert(
    const Station &station, const std::string &phase,
    std::unique_ptr<const ITravelTimeCalculator> &&travelTimeCalculator)
{
    if (travelTimeCalculator == nullptr)
    {
        throw std::invalid_argument("Travel time calculator is NULL");
    }
    if (!station.haveNetwork())
    {
        throw std::invalid_argument("Station network code not set");
    }
    if (!station.haveName())
    {
        throw std::invalid_argument("Station name not set");
    }
    if (phase.empty()){throw std::invalid_argument("Phase not defined");}
    if (contains(station, phase))
    {
        auto stationName = station.getNetwork() + "." + station.getName();
        throw std::runtime_error(stationName + "." + phase + " already exists");
    }
    pImpl->mTravelTimeCalculatorsMap.insert(
        std::pair{ std::pair{station, phase}, std::move(travelTimeCalculator)});
}

/// Contains?
bool TravelTimeCalculatorMap::contains(
    const Station &station, const std::string &phase) const
{
    if (pImpl->mTravelTimeCalculatorsMap.empty()){return false;}
    auto element = pImpl->mTravelTimeCalculatorsMap.find(std::pair{station, phase});
    return (element != pImpl->mTravelTimeCalculatorsMap.end());
}

const ITravelTimeCalculator * 
    TravelTimeCalculatorMap::at(const Station &station,
                                const std::string &phase) const
{
    if (!contains(station, phase))
    {
        auto stationPhasePair = "(" + station.getNetwork()
                              + "." + station.getName() + "," + phase + ")";
        throw std::invalid_argument("Station/phase pair: " + stationPhasePair
                                  + " does not exist");
    }
    return pImpl->mTravelTimeCalculatorsMap.at(std::pair {station, phase}).get();
}

/// Move assignment 
TravelTimeCalculatorMap& 
TravelTimeCalculatorMap::operator=(TravelTimeCalculatorMap &&map) noexcept
{
    if (&map == this){return *this;}
    pImpl = std::move(map.pImpl);
    return *this;
}
