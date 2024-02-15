#ifndef ULOCATOR_SEARCH_STATIONS_HPP
#define ULOCATOR_SEARCH_STATIONS_HPP
#ifndef NDBEBUG
#include <cassert>
#endif
#include <algorithm>
#include <limits>
#include "uLocator/travelTimeCalculator.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/optimizers/optimizer.hpp"
#include "originTime.hpp"

namespace
{
[[nodiscard]] std::pair<ULocator::Origin, double>
    searchStations(const ULocator::Origin &origin,
                   const ULocator::TravelTimeCalculatorMap &travelTimeCalculatorMap,
                   double defaultDepth,
                   const ULocator::Optimizers::IOptimizer::Norm norm,
                   const double p = 1.5,
                   const double timeWindow = 150,
                   const bool applyCorrection = true)
{
    // Find the earliest arrival and we'll use that station's
    // location as the epicenter
    auto arrivals = origin.getArrivalsReference();
    if (arrivals.empty())
    {
        throw std::invalid_argument("No arrivals to search!");
    }
    std::vector<std::pair<ULocator::Station, std::string>> stationPhases;
    stationPhases.reserve(arrivals.size());
    double earliestArrivalTime{std::numeric_limits<double>::max()};
    int earliestArrivalIndex{-1};
    for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
    {
        const auto &stationReference = arrivals.at(i).getStationReference();
        stationPhases.push_back(std::pair{stationReference,
                                          arrivals[i].getPhase()}); 
        if (arrivals[i].getTime() < earliestArrivalTime)
        {
            earliestArrivalTime = arrivals[i].getTime();
            earliestArrivalIndex = i;
        }
    }
    // Set the closest station as the source position
    auto [xSource, ySource]
        = arrivals.at(earliestArrivalIndex).getStationReference()
                                           .getLocalCoordinates();
    auto zSource = defaultDepth;
    constexpr double zeroOriginTime{0};
    // Tabulate travel times
    std::vector<double> travelTimes;
    travelTimeCalculatorMap.evaluate(stationPhases,
                                     zeroOriginTime, xSource, ySource, zSource,
                                     &travelTimes,
                                     applyCorrection);
    // Optimize for an origin time
    ULocator::Origin candidateOrigin{origin};
    candidateOrigin.setEpicenter(
        arrivals.at(earliestArrivalIndex).getStationReference()
                                         .getGeographicPosition());
    candidateOrigin.setDepth(defaultDepth); 
    // Optimize for origin time
    double originTime{0};
    if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
    {
        originTime = ::optimizeOriginTimeLeastSquares(candidateOrigin,
                                                      travelTimes);
    }
    else if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
    {
        originTime = ::optimizeOriginTimeL1(candidateOrigin, travelTimes);
    }
    else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
    {
        originTime = ::optimizeOriginTimeLp(candidateOrigin, travelTimes, p,
                                            timeWindow);
    }
    candidateOrigin.setTime(originTime);
    // Update the travel times with the origin
    std::transform(travelTimes.begin(), travelTimes.end(), travelTimes.begin(),
                   [&](const double t)
                   {
                       return originTime + t;
                   });
    for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
    {
        arrivals[i].setResidual(arrivals[i].getTime() - travelTimes[i]);
    }
    candidateOrigin.setArrivals(arrivals);
    // Tabulate the objective function
    double objectiveFunction{std::numeric_limits<double>::max()};
    if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
    {
        objectiveFunction = ::leastSquares(candidateOrigin, travelTimes,
                                           ::Measurement::Standard);
    }
    else if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
    {
        objectiveFunction = ::l1(candidateOrigin, travelTimes,
                                 ::Measurement::Standard);
    }
    else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
    {
        objectiveFunction = ::lp(candidateOrigin, travelTimes, p,
                                 ::Measurement::Standard);
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
    return std::pair {candidateOrigin, objectiveFunction};
}
}
#endif
