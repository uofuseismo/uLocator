#include <string>
#include <limits>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include <pagmo/problem.hpp>
#include <pagmo/population.hpp>
#include <pagmo/algorithms/pso.hpp>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
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
    ParticleSwarmImpl(
        std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    std::pair<double, double> mExtentInX;
    std::pair<double, double> mExtentInY;
    double mTimeWindow{120};
    double mInitialDepth{6500};
    double mMaximumDepth{65000};
    double mOptimalObjectiveFunction{std::numeric_limits<double>::max()};
    int mParticles{50};
    int mGenerations{250};
    bool mHaveOrigin{false};
    bool mHaveExtentInX{false};
    bool mHaveExtentInY{false};
};

/// Constructor
ParticleSwarm::ParticleSwarm() :
    ULocator::Optimizers::IOptimizer(),
    pImpl(std::make_unique<ParticleSwarmImpl> ())
{
}

/// Constructor
ParticleSwarm::ParticleSwarm(std::shared_ptr<UMPS::Logging::ILog> &logger) :
    ULocator::Optimizers::IOptimizer(logger),
    pImpl(std::make_unique<ParticleSwarmImpl> (logger)) 
{
}

/// Destructor
ParticleSwarm::~ParticleSwarm() = default;

void ParticleSwarm::setOriginTimeSearchWindowDuration(
    const double duration)
{
    if (duration <= 0)
    {   
        throw std::invalid_argument("Duration must be positive");
    }
    pImpl->mTimeWindow = duration;
}

double ParticleSwarm::getOriginTimeSearchWindowDuration() const noexcept
{
    return pImpl->mTimeWindow;
}

void ParticleSwarm::setExtentInX(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {
        throw std::invalid_argument("extent.second <= extent.first in x");
    }
    pImpl->mExtentInX = extent;
    pImpl->mHaveExtentInX = true;
}


void ParticleSwarm::setExtentInY(
    const std::pair<double, double> &extent)
{
    if (extent.second <= extent.first)
    {
        throw std::invalid_argument("extent.second <= extent.first in y");
    }
    pImpl->mExtentInY = extent;
    pImpl->mHaveExtentInY = true;
}

std::pair<double, double> ParticleSwarm::getExtentInX() const
{
    if (!pImpl->mHaveExtentInX)
    {
        throw std::runtime_error("Extent in x not set");
    }
    return pImpl->mExtentInX;
}

std::pair<double, double> ParticleSwarm::getExtentInY() const
{
    if (!pImpl->mHaveExtentInY)
    {
        throw std::runtime_error("Extent in y not set");
    }
    return pImpl->mExtentInY;
}

void ParticleSwarm::setNumberOfParticles(const int nParticles)
{
    if (nParticles < 1)
    {
        throw std::runtime_error("Number of particles must be positive");
    }
    pImpl->mParticles = nParticles;
}

int ParticleSwarm::getNumberOfParticles() const noexcept
{
    return pImpl->mParticles;
}

void ParticleSwarm::setNumberOfGenerations(const int nGenerations)
{
    if (nGenerations < 1)
    {
        throw std::runtime_error("Number of generations must be postiive");
    }
    pImpl->mGenerations = nGenerations;
}

int ParticleSwarm::getNumberOfGenerations() const noexcept
{
    return pImpl->mGenerations;
}

/// Locates the event
void ParticleSwarm::locate(
    const ULocator::Origin &initialGuess,
    const IOptimizer::LocationProblem locationProblem,
    const IOptimizer::Norm norm)
{
    constexpr bool reduceTimes{true};
    pImpl->mHaveOrigin = false;
    // Throws
    auto region = ULocator::Optimizers::IOptimizer::getGeographicRegion();
    const auto &arrivals
        = ULocator::Optimizers::IOptimizer::getArrivalsReference();
    if (arrivals.empty())
    {
        throw std::runtime_error("No arrivals");
    }
    // Ensure extent is set
    if (!pImpl->mHaveExtentInX)
    {
        setExtentInX(region->getExtentInX());
    }
    if (!pImpl->mHaveExtentInY)
    {
        setExtentInY(region->getExtentInY());
    }
    // Figure out the bounds
    double z0 =-ULocator::Optimizers::IOptimizer::getTopography()->
                getMinimumAndMaximumElevation().first;
    double z1 = pImpl->mMaximumDepth;
    double t0 =-getOriginTimeSearchWindowDuration(); //pImpl->mTimeWindow;
    double t1 = 0; // Reduced arrival time - first arrival is t = 0
    auto [x0, x1] = getExtentInX();
    auto [y0, y1] = getExtentInY();
/*

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
*/

    // Will use origin information after convergence
    Origin origin;
    ULocator::Position::WGS84 epicenter;
    std::vector<double> estimateArrivalTimes(arrivals.size(), 0);
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
        fitness.setSearchBoundaries(std::vector<double> {t0, x0, y0, z0},
                                    std::vector<double> {t1, x1, y1, z1} );
        // Copy observations/weights/etc.
        fitness.fillObservationsWeightsStationPhasesArraysFromArrivals(
             arrivals, reduceTimes);
        // Instantiate the problem
        pagmo::problem problem{std::move(fitness)};
        // Instantiate the particle swarm algorithm
        pagmo::algorithm algorithm{pagmo::pso(getNumberOfGenerations())};
        // Instantiate the population
        auto populationSize = static_cast<size_t> (getNumberOfParticles());
        pagmo::population population{problem, populationSize};
        // Evolve the population
        pImpl->mLogger->debug("Beginning PSO for 3D and time");
        auto newPopulation = algorithm.evolve(population);
        pImpl->mLogger->debug("PSO finished!");
        // Pick a winner and extract the hypocenter and origin time
        auto optimumLocation = newPopulation.champion_x();
        pImpl->mOptimalObjectiveFunction = newPopulation.champion_f().at(0);
        // Extract origin and Compute the theoretical arrivals
        pImpl->mLogger->debug("Computing estimate travel times");
        auto fitnessPtr
            = reinterpret_cast<const ::PagmoProblem3DAndTime *>
              (problem.get_ptr());
        origin = fitnessPtr->locationToOrigin(optimumLocation, *region);
        fitnessPtr->mTravelTimeCalculatorMap->evaluate(
            fitnessPtr->mStationPhases,
            origin.getTime(),
            optimumLocation.at(1),
            optimumLocation.at(2),
            optimumLocation.at(3),
            &estimateArrivalTimes,
            fitnessPtr->mApplyCorrection);
    }
    else if (locationProblem ==
             IOptimizer::LocationProblem::FixedToFreeSurfaceAndTime)
    {
        // Instantiate the fitness class 
        PagmoProblem2DAndTimeAndDepthAtFreeSurface fitness;
        fitness.mNorm = norm; 
        fitness.mTravelTimeCalculatorMap
            = ULocator::Optimizers::IOptimizer::getTravelTimeCalculatorMap();
        fitness.mTopography
            = ULocator::Optimizers::IOptimizer::getTopography();
        // Set bounds 
        fitness.setSearchBoundaries(std::vector<double> {t0, x0, y0},
                                    std::vector<double> {t1, x1, y1} );
        // Copy observations/weights/etc.
        fitness.fillObservationsWeightsStationPhasesArraysFromArrivals(
             arrivals, reduceTimes);
        // Instantiate the problem
        pagmo::problem problem{std::move(fitness)};
        // Instantiate the particle swarm algorithm
        pagmo::algorithm algorithm{pagmo::pso(getNumberOfGenerations())};
        // Instantiate the population
        auto populationSize = static_cast<size_t> (getNumberOfParticles());
        pagmo::population population{problem, populationSize};
        // Evolve the population
        pImpl->mLogger->debug("Beginning PSO for free-surface problem");
        auto newPopulation = algorithm.evolve(population);
        pImpl->mLogger->debug("PSO finished!");
        // Pick a winner and extract the hypocenter and origin time
        auto optimumLocation = newPopulation.champion_x();
        pImpl->mOptimalObjectiveFunction = newPopulation.champion_f().at(0);
        // Extract the origin information and compute travel times 
        pImpl->mLogger->debug("Computing estimate travel times");
        auto fitnessPtr
            = reinterpret_cast<
                 const ::PagmoProblem2DAndTimeAndDepthAtFreeSurface *>
              (problem.get_ptr());
        origin = fitnessPtr->locationToOrigin(optimumLocation, *region);
        // Compute the theoretical arrivals
        fitnessPtr->mTravelTimeCalculatorMap->evaluate(
            fitnessPtr->mStationPhases,
            origin.getTime(),
            optimumLocation.at(1),
            optimumLocation.at(2),
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
        double sourceDepth = std::min(std::max(z0, pImpl->mInitialDepth), z1);
        if (initialGuess.haveDepth())
        {
            sourceDepth = std::min(std::max(z0, initialGuess.getDepth()), z1);
        }
        fitness.mDepth = sourceDepth;
        // Set bounds 
        fitness.setSearchBoundaries(std::vector<double> {t0, x0, y0},
                                    std::vector<double> {t1, x1, y1} );
        // Copy observations/weights/etc.
        fitness.fillObservationsWeightsStationPhasesArraysFromArrivals(
             arrivals, reduceTimes);
        // Instantiate the problem
        pagmo::problem problem{std::move(fitness)};
        // Instantiate the particle swarm algorithm
        pagmo::algorithm algorithm{pagmo::pso(getNumberOfGenerations())};
        // Instantiate the population
        auto populationSize = static_cast<size_t> (getNumberOfParticles());
        pagmo::population population{problem, populationSize};
        // Evolve the population
        pImpl->mLogger->debug("Beginning PSO for fixed-depth problem");
        auto newPopulation = algorithm.evolve(population);
        pImpl->mLogger->debug("PSO finished!");
        // Pick a winner and extract the hypocenter and origin time
        auto optimumLocation = newPopulation.champion_x();
        pImpl->mOptimalObjectiveFunction = newPopulation.champion_f().at(0);
        // Extract the origin information and compute theoretical times
         pImpl->mLogger->debug("Computing estimate travel times");
        auto fitnessPtr
            = reinterpret_cast<
                 const ::PagmoProblem2DAndTimeAndFixedDepth *>
              (problem.get_ptr());
        origin = fitnessPtr->locationToOrigin(optimumLocation, *region);
        fitnessPtr->mTravelTimeCalculatorMap->evaluate(
            fitnessPtr->mStationPhases,
            origin.getTime(),
            optimumLocation.at(1),
            optimumLocation.at(2),
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
        newArrivals[i].setResidual(arrivals[i].getTime()
                                 - estimateArrivalTimes[i]);
    }
    origin.setArrivals(std::move(newArrivals));
    // Sets the origin
    IOptimizer::setOrigin(origin);
    pImpl->mHaveOrigin = true;
}

/// Optimal objective function value
double ParticleSwarm::getOptimalObjectiveFunction() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return pImpl->mOptimalObjectiveFunction;
}

/// Number of objective function evaluations
int ParticleSwarm::getNumberOfObjectiveFunctionEvaluations() const
{
    if (!haveOrigin()){throw std::runtime_error("Event not located");}
    return getNumberOfParticles()*getNumberOfGenerations();
}

bool ParticleSwarm::haveOrigin() const noexcept
{
    return pImpl->mHaveOrigin;
}
