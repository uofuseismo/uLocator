#include <string>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include <pagmo/problem.hpp>
#include <pagmo/population.hpp>
#include <pagmo/algorithms/pso.hpp>
#include "uLocator/optimizers/pagmo/particleSwarm.hpp"
#include "uLocator/position/wgs84.hpp"
#include "fitnessFunctions.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/topography/topography.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"
#include "uLocator/origin.hpp"

using namespace ULocator::Optimizers::Pagmo;

class ParticleSwarm::ParticleSwarmImpl
{
public:
    double mTimeWindow{50};
    double mInitialDepth{6500};
    int mParticles{500};
    int mGenerations{200};
    bool mHaveOrigin{false};
};

/// Constructor
ParticleSwarm::ParticleSwarm() :
    ULocator::Optimizers::IOptimizer(),
    pImpl(std::make_unique<ParticleSwarmImpl> ())
{
}

/// Destructor
ParticleSwarm::~ParticleSwarm() = default;

/// Locates with 

/// Locates 
void ParticleSwarm::locate(
    const IOptimizer::LocationProblem locationProblem,
    IOptimizer::Norm norm)
{
    pImpl->mHaveOrigin = false;
    // Throws
    auto region = ULocator::Optimizers::IOptimizer::getGeographicRegion();
    const auto &arrivals
        = ULocator::Optimizers::IOptimizer::getArrivalsReference();
    if (arrivals.empty())
    {
        throw std::runtime_error("No arrivals");
    }
    // Figure out the bounds
    double z0 =-ULocator::Optimizers::IOptimizer::getTopography()->
                getMinimumAndMaximumElevation().first;
    double z1 = 65000;
    double t0 =-pImpl->mTimeWindow;
    double t1 = 0; // Reduced arrival time - first arrival is t = 0
    auto [x0, x1] = region->getExtentInX();
    auto [y0, y1] = region->getExtentInY();
    double reductionTime{0};

    // Extract and reduce observations
    std::vector<std::pair<ULocator::Station, std::string>>
        stationPhases(arrivals.size());
    std::vector<double> weights(arrivals.size());
    std::vector<double> observations(arrivals.size());
    std::vector<double> estimateArrivalTimes(arrivals.size(), 0);
    for (int i = 0; i < static_cast<int> (arrivals.size()); ++i)
    {
        stationPhases[i]
            = std::pair{arrivals[i].getStation(), arrivals[i].getPhase()};
        weights[i] = 1./arrivals[i].getStandardError();
        observations[i] = arrivals[i].getTime();
    }

    Origin origin;
    ULocator::Position::WGS84 epicenter;
    double originTime{0};
    double sourceDepth{std::min(std::max(z0, pImpl->mInitialDepth), z1)};

    if (locationProblem == IOptimizer::LocationProblem::ThreeDimensionsAndTime)
    {
        // Instantiate the fitness class 
        ::PagmoProblem3DAndTime fitness;
        fitness.mNorm = norm; 
        fitness.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        fitness.mTopography
            = ULocator::Optimizers::IOptimizer::getTopography();
        // Set bounds 
        fitness.mLowerBounds = pagmo::vector_double {t0, x0, y0, z0};
        fitness.mUpperBounds = pagmo::vector_double {t1, x1, y1, z1};
        // Copy observations/weights/etc.
        fitness.mStationPhases = std::move(stationPhases);
        fitness.mWeights = std::move(weights);
        ::reduceTimes(observations, &fitness.mObservations, &reductionTime);
        // Instantiate the problem
        pagmo::problem problem{std::move(fitness)};
        // Instantiate the particle swarm algorithm
        pagmo::algorithm algorithm{pagmo::pso(pImpl->mGenerations)};
        // Instantiate the population
        auto populationSize = static_cast<size_t> (pImpl->mParticles);
        pagmo::population population{problem, populationSize};
        // Evolve the population
        auto newPopulation = algorithm.evolve(population);
        // Pick a winner and extract the hypocenter and origin time
        auto optimumLocation = newPopulation.champion_x();
        auto optimumObjectiveFunction = newPopulation.champion_f();
        // Extract the origin information
        origin.setTime(reductionTime + optimumLocation.at(0));
        double xSource{optimumLocation.at(1)};
        double ySource{optimumLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region->localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(optimumLocation.at(3));
        // Compute the theoretical arrivals
        auto fitnessPtr
            = reinterpret_cast<const ::PagmoProblem3DAndTime *>
              (problem.get_ptr());
        fitnessPtr->mTravelTimeCalculatorMap->evaluate(
            fitnessPtr->mStationPhases,
            origin.getTime(),
            xSource,
            ySource,
            origin.getDepth(),
            &estimateArrivalTimes,
            fitnessPtr->mApplyCorrection);
    }
    else if (locationProblem ==
             IOptimizer::LocationProblem::FixedToFreeSurfaceAndTime)
    {
        // Instantiate the fitness class 
        ::PagmoProblem2DAndTimeAndDepthAtFreeSurface fitness;
        fitness.mNorm = norm; 
        fitness.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        fitness.mTopography
            = ULocator::Optimizers::IOptimizer::getTopography();
        // Set bounds 
        fitness.mLowerBounds = pagmo::vector_double {t0, x0, y0};
        fitness.mUpperBounds = pagmo::vector_double {t1, x1, y1};
        // Copy observations/weights/etc.
        fitness.mStationPhases = std::move(stationPhases);
        fitness.mWeights = std::move(weights);
        ::reduceTimes(observations, &fitness.mObservations, &reductionTime);
        // Instantiate the problem
        pagmo::problem problem{std::move(fitness)};
        // Instantiate the particle swarm algorithm
        pagmo::algorithm algorithm{pagmo::pso(pImpl->mGenerations)};
        // Instantiate the population
        auto populationSize = static_cast<size_t> (pImpl->mParticles);
        pagmo::population population{problem, populationSize};
        // Evolve the population
        auto newPopulation = algorithm.evolve(population);
        // Pick a winner and extract the hypocenter and origin time
        auto optimumLocation = newPopulation.champion_x();
        auto optimumObjectiveFunction = newPopulation.champion_f();
        // Extract the origin information
        origin.setTime(reductionTime + optimumLocation.at(0));
        double xSource{optimumLocation.at(1)};
        double ySource{optimumLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region->localToGeographicCoordinates(xSource, ySource);
        auto fitnessPtr
            = reinterpret_cast<
                 const ::PagmoProblem2DAndTimeAndDepthAtFreeSurface *>
              (problem.get_ptr());
        auto elevation = fitnessPtr->mTopography->evaluate(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(-elevation);
        // Compute the theoretical arrivals
        fitnessPtr->mTravelTimeCalculatorMap->evaluate(
            fitnessPtr->mStationPhases,
            origin.getTime(),
            xSource,
            ySource,
            origin.getDepth(),
            &estimateArrivalTimes,
            fitnessPtr->mApplyCorrection);
    }
    else if (locationProblem == IOptimizer::LocationProblem::FixedDepthAndTime)
    {
        // Instantiate the fitness class 
        ::PagmoProblem2DAndTimeAndFixedDepth fitness;
        fitness.mNorm = norm;
        fitness.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        fitness.mDepth = sourceDepth;
        // Set bounds 
        fitness.mLowerBounds = pagmo::vector_double {t0, x0, y0};
        fitness.mUpperBounds = pagmo::vector_double {t1, x1, y1};
        // Copy observations/weights/etc.
        fitness.mStationPhases = std::move(stationPhases);
        fitness.mWeights = std::move(weights);
        ::reduceTimes(observations, &fitness.mObservations, &reductionTime);
        // Instantiate the problem
        pagmo::problem problem{std::move(fitness)};
        // Instantiate the particle swarm algorithm
        pagmo::algorithm algorithm{pagmo::pso(pImpl->mGenerations)};
        // Instantiate the population
        auto populationSize = static_cast<size_t> (pImpl->mParticles);
        pagmo::population population{problem, populationSize};
        // Evolve the population
        auto newPopulation = algorithm.evolve(population);
        // Pick a winner and extract the hypocenter and origin time
        auto optimumLocation = newPopulation.champion_x();
        auto optimumObjectiveFunction = newPopulation.champion_f();
        // Extract the origin information
        origin.setTime(reductionTime + optimumLocation.at(0));
        double xSource{optimumLocation.at(1)};
        double ySource{optimumLocation.at(2)};
        auto [sourceLatitude, sourceLongitude]
            = region->localToGeographicCoordinates(xSource, ySource);
        origin.setEpicenter(ULocator::Position::WGS84 {sourceLatitude,
                                                       sourceLongitude});
        origin.setDepth(sourceDepth);
        // Compute the theoretical arrivals
        auto fitnessPtr
            = reinterpret_cast<
                 const ::PagmoProblem2DAndTimeAndFixedDepth *>
              (problem.get_ptr());
        fitnessPtr->mTravelTimeCalculatorMap->evaluate(
            fitnessPtr->mStationPhases,
            origin.getTime(),
            xSource,
            ySource,
            origin.getDepth(),
            &estimateArrivalTimes,
            fitnessPtr->mApplyCorrection);
    }
    else
    {
        throw std::runtime_error("Unhandled location problem");
    }
    // Insert the arrival times
    auto newArrivals = arrivals;
    for (int i = 0; i < static_cast<int> (newArrivals.size()); ++i)
    {
        newArrivals[i].setResidual(observations[i] - estimateArrivalTimes[i]);
    }
    origin.setArrivals(std::move(newArrivals));
    // Sets the origin
    IOptimizer::setOrigin(origin);
    pImpl->mHaveOrigin = true;
}

bool ParticleSwarm::haveOrigin() const noexcept
{
    return pImpl->mHaveOrigin;
}
