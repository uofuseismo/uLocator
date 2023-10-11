#include <iostream>
#include <fstream>
#include <iomanip>
#include <functional>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include <nlopt.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/tokenizer.hpp>
#include <umps/logging/standardOut.hpp>
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/firstArrivalRayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/position/wgs84.hpp"

#define UTM_ZONE 12

double qualityToStandardError(const double quality)
{
    if (std::abs(quality - 1) < 1.e-4)
    {   
        return 0.03;
    }   
    else if (std::abs(quality - 0.75) < 1.e-4)
    {   
        return 0.06;
    }   
    else if (std::abs(quality - 0.5) < 1.e-4)
    {   
        return 0.15;
    }   
#ifndef NDEBUG
    assert(false);
#endif
    return 1;
}

std::vector<ULocator::Origin>
    loadCatalog(const std::string &fileName,
                const std::string &trainingFileName)
{
    // Get training events
    std::string line;
    std::vector<int64_t> trainingEventIdentifiers;
    std::ifstream trainingFile(trainingFileName);
    getline(trainingFile, line); // Header
    while (getline(trainingFile, line))
    {   
        std::vector<std::string> splitLine;
        splitLine.reserve(64);
        boost::tokenizer<boost::escaped_list_separator<char>>
           tokenizer(line,
                     boost::escaped_list_separator<char>('\\', ',', '\"'));
        for (auto i = tokenizer.begin(); i != tokenizer.end(); ++i) 
        {
            splitLine.push_back(*i);
        }
        trainingEventIdentifiers.push_back(std::stol(splitLine[0]));
    }   
    // Get entire catalog
    std::vector<ULocator::Origin> origins;
    origins.reserve(10000);
    std::ifstream csvFile(fileName); 
    if (!csvFile.is_open())
    {   
        throw std::runtime_error("Couldn't open: " + fileName);
    }   
    ULocator::Origin origin;
    std::vector<ULocator::Arrival> arrivals;
    int64_t evidOld{-1};
    getline(csvFile, line); // Header
    while (getline(csvFile, line))
    {   
        std::vector<std::string> splitLine;
        splitLine.reserve(64);
        boost::tokenizer<boost::escaped_list_separator<char>>
           tokenizer(line,
                     boost::escaped_list_separator<char>('\\', ',', '\"'));
        for (auto i = tokenizer.begin(); i != tokenizer.end(); ++i) 
        {
            splitLine.push_back(*i);
        }
        auto evid = std::stol(splitLine[0]);
        auto network = splitLine[1];
        auto stationName = splitLine[2];
        auto phase = splitLine[5];
        auto arid = std::stol(splitLine[6]);
        auto time = std::stod(splitLine[7]);
        auto residual = std::stod(splitLine[13]);
        auto uncertainty = qualityToStandardError(std::stod(splitLine[8]));
        auto stationLat = std::stod(splitLine[14]);
        auto stationLon = std::stod(splitLine[15]);
        auto stationElev = std::stod(splitLine[16]);
        auto eventLat = std::stod(splitLine[17]); 
        auto eventLon = std::stod(splitLine[18]);
        auto eventDepth = 1000*std::stod(splitLine[19]);
        auto originTime = std::stod(splitLine[20]);
        auto eventType = splitLine[24];
        if (evid != evidOld)
        {
            // When we get to a new event finish out this origin
            if (evidOld !=-1)
            {
                if (std::find(trainingEventIdentifiers.begin(),
                              trainingEventIdentifiers.end(),
                              origin.getIdentifier())
                    != trainingEventIdentifiers.end())
                {
                    origin.setArrivals(arrivals);
                    origins.push_back(origin);
                }
            }
            // Update
            evidOld = evid;
            arrivals.clear();
            origin.clear();
            // Set initial origin information
            ULocator::Position::WGS84 epicenter{eventLat, eventLon, UTM_ZONE};
            origin.setIdentifier(evid);
            origin.setEpicenter(epicenter);
            origin.setDepth(eventDepth);
            origin.setTime(originTime);
            origin.setEventType(ULocator::Origin::EventType::Earthquake);
            if (eventType == "qb")
            {
                origin.setEventType(ULocator::Origin::EventType::QuarryBlast);
            }
        }
        ULocator::Position::WGS84 stationPosition{stationLat, stationLon, UTM_ZONE};
        ULocator::Arrival arrival;
        ULocator::Station station;
        station.setNetwork(network);
        station.setName(stationName);
        station.setGeographicPosition(stationPosition);
        station.setElevation(stationElev);
        arrival.setTime(time);
        arrival.setStandardError(uncertainty);
        arrival.setResidual(residual);
        arrival.setPhase(phase);
        arrival.setStation(station);
        arrival.setIdentifier(arid);
        arrivals.push_back(std::move(arrival));
    }   
    return origins;    
}

struct Layer
{
    Layer(const double inInterface, const double inVelocity,
          const double inLowerBound, const double inUpperBound) :
        topInterface(inInterface),
        velocity(inVelocity),
        lowerBound(inLowerBound),
        upperBound(inUpperBound)
    {
        if (velocity <= 0)
        {
            throw std::invalid_argument("Velocity must be positive");
        }
        if (velocity < lowerBound)
        {
            throw std::invalid_argument(
                "Velocity cannot be less than lower bound");
        }
        if (velocity > upperBound)
        {
            throw std::invalid_argument(
                "Velocity cannot be greater than upper bound");
        }
    }
    double topInterface{0};
    double velocity{0};
    double lowerBound{0};
    double upperBound{0};
};

struct ObjectiveFunction
{
    void setOrigins(const std::vector<ULocator::Origin> &inOrigins)
    {
        origins.clear();
        arrivalTimes.clear();
        estimates.clear();
        catalogResiduals.clear();
        weights.clear();
        stationPhasePairs.clear();
        origins.reserve(inOrigins.size());
        arrivalTimes.reserve(inOrigins.size()*10);
        weights.reserve(inOrigins.size()*10);
        stationPhasePairs.reserve(inOrigins.size()); 
        double objectiveFunction = 0;
        for (const auto &inOrigin : inOrigins)
        {
            auto origin = inOrigin;
            auto arrivals = origin.getArrivals();
            std::vector<ULocator::Arrival> arrivalsToKeep;
            std::vector<std::pair<ULocator::Station, std::string>> spPair;
            for (const auto &arrival : arrivals)
            {
                if (arrival.getPhase() == phase)
                {
                    arrivalsToKeep.push_back(arrival);
                    arrivalTimes.push_back(arrival.getTime());
                    weights.push_back(1./arrival.getStandardError());
                    catalogResiduals.push_back(arrival.getResidual());
                    auto weightedResidual
                         = arrival.getResidual()/arrival.getStandardError();
                    spPair.push_back(std::pair {arrival.getStation(), phase});
                    objectiveFunction = objectiveFunction
                                      + weightedResidual*weightedResidual;
                }
            }
            if (!arrivalsToKeep.empty())
            {
                origin.setArrivals(arrivalsToKeep);
                origins.push_back(origin);
                stationPhasePairs.push_back(std::move(spPair));
            }
        }
        estimates.resize(arrivalTimes.size());
        logger->info("Catalog obj fn: " + std::to_string(objectiveFunction));
        logger->info("Number of origins: " + std::to_string(origins.size()));
        logger->info("Number of phase arrivals: "
                   + std::to_string(arrivalTimes.size()));
        // Tabulate unique station/phase pairs
        std::vector<std::string> names;
        for (const auto &origin : origins)
        {
            const auto &arrivals = origin.getArrivalsReference();
            for (const auto &arrival : arrivals)
            {
                if (arrival.getPhase() != phase){continue;}
                const auto &station = arrival.getStationReference();
                bool found = false;
                for (const auto &uniqueStation : uniqueStations)
                {
                    if (uniqueStation.getNetwork() == station.getNetwork() &&
                        uniqueStation.getName()  == station.getName())
                    {
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    uniqueStations.push_back(station);
                }
            }
        }
        // Now make a list of station/phase pairs
        logger->info("Number of unique stations: "
                   + std::to_string(uniqueStations.size()));
    }
    void initializeCalculators()
    {
        logger->info("Initializing calculators...");
#ifndef NDEBUG
        assert(interfaces.size() == velocities.size());
#endif
        calculators.clear();
        for (const auto &uniqueStation : uniqueStations)
        {
            auto name = uniqueStation.getNetwork()
                      + "." + uniqueStation.getName() + phase;
            auto calculator
                = std::make_unique<ULocator::FirstArrivalRayTracer> (logger);
            try
            {
                calculator->initialize(
                    uniqueStation, phase, interfaces, velocities);
                calculators.insert(uniqueStation, phase, std::move(calculator));
            }
            catch (const std::exception &e)
            {
                 logger->error(e.what());
            }
        }
    }
    double f(const int n, const double x[])
    {
        nEvaluations = nEvaluations + 1;
        logger->info("Objective function evaluation: "
                   + std::to_string(nEvaluations));
        constexpr bool applyCorrection{false};
        // Reinitialize the travel time map
        if (static_cast<int> (velocities.size()) != n)
        {
            velocities.resize(n);
        }
        //for (int i =0; i < n; ++i){std::cout << x[i] << std::endl;}
        //getchar();
        std::copy(x, x + n, velocities.data());
        initializeCalculators();
        // Tabulate the objective function
        double objectiveFunction = 0;
        int nObservations{0};
        for (int iOrigin = 0;
             iOrigin < static_cast<int> (origins.size()); ++iOrigin)
        {
            auto epicenter = origins[iOrigin].getEpicenter();
            auto depth = origins[iOrigin].getDepth();
            auto time = origins[iOrigin].getTime();
            const auto &arrivals = origins[iOrigin].getArrivalsReference();
            std::vector<double> travelTimes;
            calculators.evaluate(stationPhasePairs[iOrigin], epicenter, depth,
                                 &travelTimes, applyCorrection);
            for (int iArrival = 0;
                 iArrival < static_cast<int> (arrivals.size()); ++iArrival)
            {
                auto observed = arrivalTimes[nObservations]; //arrivals[iArrival].getTime();
                auto weight = weights[nObservations]; //1./arrivals[iArrival].getStandardError();
                double estimate = travelTimes[iArrival] + time;
                auto weightedResidual = weight*(observed - estimate);
                objectiveFunction = objectiveFunction
                                  + weightedResidual*weightedResidual;
                estimates[nObservations] = estimate;
                nObservations = nObservations + 1;
            }
        }
#ifndef NDEBUG
        assert(nObservations == arrivalTimes.size());
#endif
        logger->info("Objective function: "
                   + std::to_string(objectiveFunction));
        return objectiveFunction;
    }
    double operator()(const unsigned int n, const double *x, double *g)
    {
        if (g == nullptr)
        {
            return f(n, x);
        }
        throw std::runtime_error("gradient not done");
    } 
    void writeToCSV(const std::string &filename)
    {
        std::ofstream of(filename);
        of << "observed,estimate,residual,catalogResidual,weight" << std::endl;
        for (int i = 0; i < estimates.size(); ++i)
        {
            of << std::setprecision(16) << arrivalTimes[i]
               << "," << estimates[i] 
               << "," << arrivalTimes[i] - estimates[i] 
               << "," << catalogResiduals[i]
               << "," << weights[i]
               << std::endl;
        }
        of.close();
    }
    std::shared_ptr<UMPS::Logging::ILog> logger{nullptr};
    std::vector<ULocator::Station> uniqueStations;
    std::vector<ULocator::Origin> origins;
    std::vector< std::vector<std::pair<ULocator::Station, std::string>> >
        stationPhasePairs;
    std::vector<double> arrivalTimes;
    std::vector<double> catalogResiduals;
    std::vector<double> estimates;
    std::vector<double> weights;
    std::vector<double> interfaces;
    std::vector<double> velocities;
    ULocator::TravelTimeCalculatorMap calculators;
    std::string phase{"P"};
    int nEvaluations{0};
};

int main(int argc, char *argv[])
{
    auto logger = std::make_shared<UMPS::Logging::StandardOut> ();

    std::string phase{"P"};
    std::string initialResidualsFile{"initialResiduals." + phase + ".csv"};
    std::string finalResidualsFile{"finalResiduals." + phase + ".csv"};
    std::string catalogFile{"../examples/uuss/utah_catalog.csv"};
    std::string trainingFile{"../examples/uuss/utah_events_training.csv"};
    auto origins = loadCatalog(catalogFile, trainingFile);

    // Make some baseline
    std::vector<::Layer> initialModel{
        ::Layer{-4500, 3400, 2000, 4500},
        ::Layer{   40, 5900, 4500, 6000},
        ::Layer{15600, 6400, 6000, 7000},
        ::Layer{26500, 7500, 7000, 7800},
        ::Layer{40500, 7900, 7800, 8400}};
    if (phase == "S")
    {
        initialModel = std::vector<::Layer>
        {
            ::Layer{-4500, 1973, 1000, 2300},
            ::Layer{   40, 3490, 2300, 3500},
            ::Layer{15600, 3804, 3500, 4000},
            ::Layer{26500, 4115, 4000, 4500},
            ::Layer{40500, 5225, 4500, 5500}
        };
    }

    ::ObjectiveFunction objectiveFunction;
    objectiveFunction.logger = logger;
    objectiveFunction.phase = phase;
    objectiveFunction.setOrigins(origins);
    std::vector<double> xInitial;
    for (const auto &layer : initialModel)
    {
        xInitial.push_back(layer.velocity);
        objectiveFunction.interfaces.push_back(layer.topInterface);
    }
    auto f0 = objectiveFunction.f(xInitial.size(), xInitial.data());
    objectiveFunction.writeToCSV(initialResidualsFile);
    logger->info("Initial objective function: " + std::to_string(f0));
    objectiveFunction.nEvaluations = 0;

    std::vector<::Layer> layers{
        ::Layer{-4500, 3695, 2000, 5000},
        ::Layer{   40, 5975, 5000, 6200},
        ::Layer{15600, 6565, 6200, 6800},
        ::Layer{26500, 7010, 6800, 7600},
        ::Layer{40500, 7690, 7600, 8400}};
    if (phase == "S")
    {
        layers = std::vector<::Layer>
        {
            ::Layer{-4500, 2180, 1000, 2300},
            ::Layer{   40, 3435, 2300, 3500},
            ::Layer{15600, 3660, 3500, 4000},
            ::Layer{26500, 4360, 4000, 4500},
            ::Layer{40500, 5000, 4500, 5500}
        };
    }
/*
{ 
        ::Layer{-4500,  2000,  500, 2700},
        ::Layer{-2500,  2800, 2700, 3100},
        ::Layer{-1500,  3300, 3100, 3900},
        ::Layer{ -500,  4000, 3900, 4800},
        ::Layer{ -200,  5000, 4800, 5700},
        ::Layer{   40,  5900, 5700, 6100},
        ::Layer{10000,  6200, 6100, 6300},
        ::Layer{15600,  6400, 6300, 6600},
        ::Layer{20500,  6900, 6600, 7200},
        ::Layer{26500,  7500, 7200, 7600},
        ::Layer{36500,  7700, 7600, 7800},
        ::Layer{40500,  7900, 7800, 8500}};
*/

    auto n = static_cast<int> (layers.size());
    std::vector<double> x;
    std::vector<double> interfaces;
    std::vector<double> lowerBounds;
    std::vector<double> upperBounds;
    for (const auto &layer : layers)
    {
        interfaces.push_back(layer.topInterface);
        x.push_back(layer.velocity);
        lowerBounds.push_back(layer.lowerBound);
        upperBounds.push_back(layer.upperBound);
    } 
    objectiveFunction.interfaces = interfaces;


    //nlopt::opt optimizer(nlopt::GN_DIRECT_L, n);
    nlopt::opt optimizer(nlopt::LN_BOBYQA, n);
    optimizer.set_maxeval(100);
    optimizer.set_lower_bounds(lowerBounds);
    optimizer.set_upper_bounds(upperBounds);
    optimizer.set_xtol_abs(1);
    auto nlOptObjectiveFunction
         = std::bind(&::ObjectiveFunction::operator(),
                     &objectiveFunction,
                     std::placeholders::_1,
                     std::placeholders::_2,
                     std::placeholders::_3);
    optimizer.set_min_objective(nlOptObjectiveFunction);
    double loss{0};
    try
    {
        optimizer.optimize(x, loss);
    }
    catch (const std::exception &e)
    {
        logger->error(e.what());
    }
    std::cout << "Final loss: " << loss << std::endl;
    for (int i = 0; i < objectiveFunction.interfaces.size(); ++i)
    {
        std::cout << objectiveFunction.interfaces[i] << " "
                  << x[i] << std::endl;
    }
    objectiveFunction.nEvaluations = 0;
    auto fFinal = objectiveFunction.f(x.size(), x.data());
    std::cout << "Evaluating final loss: " << fFinal << "; "
              << (f0 - fFinal)/f0*100 <<  " pct reduction" << std::endl;
    objectiveFunction.writeToCSV(finalResidualsFile);
    return EXIT_SUCCESS;
}
