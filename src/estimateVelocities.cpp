#include <iostream>
#include <filesystem>
#include <string>
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
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#endif
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
                const std::string &trainingFileName,
                std::shared_ptr<UMPS::Logging::ILog> logger)
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
        // Assert the arrival does not already exist
        bool exists{false};
        for (const auto &a : arrivals)
        {
            const auto &si = a.getStation();
            if (si.getNetwork() == network &&
                si.getName()    == stationName &&
                a.getPhase()    == phase)
            {
                exists = true;
                logger->warn("Duplicate phase arrival for: "
                            + std::to_string(evid)
                            + " " + network + "." + stationName + "." + phase
                            + "; skipping...");
            }
        }
        if (!exists)
        {
            arrivals.push_back(std::move(arrival));
        }
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
            auto epicenter = origins.at(iOrigin).getEpicenter();
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
#ifdef WITH_UMPS
    std::shared_ptr<UMPS::Logging::ILog> logger{nullptr};
#endif
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

struct ProgramOptions
{
    std::filesystem::path resultsDirectory{"/home/bbaker/Codes/uLocator/examples/uuss/ynpResults"};
    std::filesystem::path catalogFile{"/home/bbaker/Codes/uLocator/examples/uuss/ynpResults/utah_catalog.csv"};
    std::filesystem::path trainingFile{"/home/bbaker/Codes/uLocator/examples/uuss/ynpResults/utah_events_training.csv"};
    std::string phase{"P"};
    bool doUtah{true};
    bool doHelp{false};
};

[[nodiscard]] ::ProgramOptions parseCommandLineOptions(int argc, char *argv[])
{
    ::ProgramOptions options;
    boost::program_options::options_description desc(
R"""(
This utility will help to estimate a new starting model for Utah or Yellowstone.  Example usage:
     estimateVelocities --results_directory=../examples/uuss/ynpResults --catalog_file=../examples/uuss/ynp_catalog.csv --training_file=../examples/uuss/ynp_events_training.csv --phase=P --do_ynp
Allowed options)""");
    desc.add_options()
        ("help",         "Produces this help message")
        ("results_directory", boost::program_options::value<std::string> ()->default_value(options.resultsDirectory),
                          "The directory to which results will be written")
        ("catalog_file",  boost::program_options::value<std::string> ()->default_value(options.catalogFile),
                         "A CSV with the events and travel times")
        ("training_file",  boost::program_options::value<std::string> ()->default_value(options.trainingFile),
                         "A CSV file with the events identifiers in the training set")
        ("phase", boost::program_options::value<std::string> ()->default_value(options.phase),
                  "The phase velocity for which we are inverting e.g., P or S")
        ("do_ynp", "If present then this will use the default Yellowstone models.  Otherwise, this will be for optimizing a Utah model.");
    boost::program_options::variables_map vm;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, desc), vm);
    boost::program_options::notify(vm);
    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        options.doHelp = true;
        return options;
    }
    if (vm.count("catalog_file"))
    {   
        auto catalogFile = vm["catalog_file"].as<std::string> (); 
        if (!std::filesystem::exists(catalogFile))
        {
            throw std::runtime_error(catalogFile + " does not exist");
        }
        options.catalogFile = catalogFile;
    }   
    else 
    {    
        throw std::runtime_error("Catalog file not set");
    }
    if (vm.count("training_file"))
    {   
        auto trainingFile = vm["training_file"].as<std::string> (); 
        if (!std::filesystem::exists(trainingFile))
        {   
            throw std::runtime_error(trainingFile + " does not exist");
        }   
        options.trainingFile = trainingFile;
    }   
    else 
    {       
        throw std::runtime_error("Training file not set");
    }
    if (vm.count("results_directory"))
    {   
        auto resultsDirectory = vm["results_directory"].as<std::string> ();
        if (!std::filesystem::exists(resultsDirectory))
        {
            if (!std::filesystem::create_directories(resultsDirectory) &&
                !resultsDirectory.empty())
            {
                throw std::runtime_error(resultsDirectory + " could not be made");
            }
        }
        options.resultsDirectory = resultsDirectory;
    }
    if (vm.count("phase"))
    {
        auto phase = vm["phase"].as<std::string> ();
        if (phase == "P" or phase == "S")
        {
            options.phase = phase;
        }
        else
        {
            throw std::invalid_argument("Unhandled phase " + phase);
        }
    }
    options.doUtah = true;
    if (vm.count("do_ynp")){options.doUtah = false;}
    return options;
}

int main(int argc, char *argv[])
{
    ::ProgramOptions programOptions;
    try
    {
        programOptions = ::parseCommandLineOptions(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    if (programOptions.doHelp){return EXIT_SUCCESS;}

#ifdef WITH_UMPS
    auto logger = std::make_shared<UMPS::Logging::StandardOut> ();
#endif

    auto phase = programOptions.phase;
    auto initialResidualsFile = programOptions.resultsDirectory / std::filesystem::path{"initialResiduals." + phase + ".csv"};
    auto finalResidualsFile  = programOptions.resultsDirectory / std::filesystem::path{"finalResiduals." + phase + ".csv"};
    auto catalogFile = programOptions.catalogFile;
    auto trainingFile = programOptions.trainingFile;
#ifdef WITH_UMPS
    auto origins = ::loadCatalog(catalogFile, trainingFile, logger);
#else
    auto origins = ::loadCatalog(catalogFile, trainingFile);
#endif

    // Make some baseline
    std::vector<::Layer> initialModel;
    if (programOptions.doUtah)
    {
        if (phase == "P")
        {
#ifdef WITH_UMPS
            logger->info("Creating initial Utah P model...");
#else
            std::cout << "Creating initial Utah P model..." << std::endl;
#endif
            initialModel = std::vector<::Layer>
            {
                ::Layer{-4500, 3400, 2000, 4500},
                ::Layer{   40, 5950, 4500, 6000},
                ::Layer{15600, 6450, 6000, 7000},
                ::Layer{26500, 7550, 7000, 7800},
                ::Layer{40500, 7900, 7800, 8400}
            };
        }
        else if (phase == "S")
        {
#ifdef WITH_UMPS
            logger->info("Creating initial Utah S model...");
#else
            std::cout << "Creating initial Utah S model..." << std::endl;
#endif
            initialModel = std::vector<::Layer>
            {
                ::Layer{-4500, 1950, 1000, 2300},
                ::Layer{   40, 3390, 2300, 3500},
                ::Layer{15600, 3680, 3500, 4000},
                ::Layer{26500, 4310, 4000, 4500},
                ::Layer{40500, 4540, 4500, 5500}
            };
        }
        else
        {
            throw std::runtime_error("Unhandled phase " + phase);
        }
    }
    else //if (!programOptions.doUtah)
    {
        if (phase == "P")
        {
#ifdef WITH_UMPS
            logger->info("Creating initial YNP P model...");
#else
            std::cout << "Creating initial YNP P model..." << std::endl;
#endif
            initialModel = std::vector<::Layer>
            {
                ::Layer{-4500,  2720, 1000, 2750},
                ::Layer{-1000,  2790, 2750, 4000},
                ::Layer{ 2000,  5210, 4000, 5350},
                ::Layer{ 5000,  5560, 5350, 5600},
                ::Layer{ 8000,  5770, 5600, 5850},
                ::Layer{12000,  6070, 5850, 6150},
                ::Layer{16000,  6330, 6150, 6450},
                ::Layer{21000,  6630, 6450, 7200},
                ::Layer{50000,  8000, 7200, 8400}
            };
        }
        else if (phase == "S")
        {
#ifdef WITH_UMPS
            logger->info("Creating initial YNP S model...");
#else
            std::cout << "Creating initial YNP S model..." << std::endl;
#endif
            initialModel = std::vector<::Layer>
            {

                ::Layer{-4500,  1950,  500, 1975},
                ::Layer{-1000,  2000, 1975, 2700},
                ::Layer{ 2000,  3400, 2700, 3410},
                ::Layer{ 5000,  3420, 3410, 3450},
                ::Layer{ 8000,  3490, 3450, 3590},
                ::Layer{12000,  3680, 3590, 3710},
                ::Layer{16000,  3780, 3710, 3890},
                ::Layer{21000,  4000, 3890, 4400},
                ::Layer{50000,  4850, 4400, 5500}
            };
        }
        else
        {
            throw std::runtime_error("Unhandled phase " + phase);
        }
    }
#ifndef NDEBUG
    assert(!initialModel.empty());
#endif

    ::ObjectiveFunction objectiveFunction;
#ifdef WITH_UMPS
    objectiveFunction.logger = logger;
#endif
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
#ifdef WITH_UMPS
    logger->info("Initial objective function: " + std::to_string(f0));
#else
    std::cout << "Initial objective function: " << f0 << std::endl;
#endif
    objectiveFunction.nEvaluations = 0;

    std::vector<::Layer> layers;
    if (programOptions.doUtah)
    {
        if (phase == "P")
        {
#ifdef WITH_UMPS
            logger->info("Creating starting Utah P model...");
#else
            std::cout << "Creating starting Utah P model..." << std::endl;
#endif
            layers = std::vector<::Layer>
            {
                ::Layer{-4500, 3695, 2000, 5000},
                ::Layer{   40, 5975, 5000, 6200},
                ::Layer{15600, 6565, 6200, 6800},
                ::Layer{26500, 7010, 6800, 7600},
                ::Layer{40500, 7690, 7600, 8400}
            };
        }
        else if (phase == "S")
        {
#ifdef WITH_UMPS
            logger->info("Creating starting Utah S model...");
#else
            std::cout << "Creating start Utah S model..." << std::endl;
#endif
            layers = std::vector<::Layer>
            {
                ::Layer{-4500, 2180, 1000, 2300},
                ::Layer{   40, 3435, 2300, 3500},
                ::Layer{15600, 3660, 3500, 4000},
                ::Layer{26500, 4360, 4000, 4500},
                ::Layer{40500, 5000, 4500, 5500}
            };
        }
    }
    else //if (!programOptions.doUtah)
    { 
        if (phase == "P")
        {
#ifdef WITH_UMPS
            logger->info("Creating starting YNP P model...");
#else
            std::cout << "Creating starting YNP P model..." << std::endl;
#endif
            layers = std::vector<::Layer>
            {
                ::Layer{-4500,  2458, 1000, 2750},
                ::Layer{-1000,  3375, 2750, 4000},
                ::Layer{ 2000,  4675, 4000, 5000},
                ::Layer{ 5000,  5475, 5000, 5600},
                ::Layer{ 8000,  5725, 5600, 5900},
                ::Layer{12000,  6133, 5900, 6250},
                ::Layer{16000,  6400, 6250, 6450},
                ::Layer{21000,  6575, 6450, 7200},
                ::Layer{50000,  8200, 7200, 8400}
            };
        }
        else if (phase == "S")
        {
#ifdef WITH_UMPS
            logger->info("Creating starting YNP S model...");
#else
            std::cout << "Creating starting YNP S model..." << std::endl;
#endif
            layers = std::vector<::Layer>
            {
                ::Layer{-4500,  1729,  500, 1975},
                ::Layer{-1000,  2338, 1975, 2700},
                ::Layer{ 2000,  3055, 2700, 3200},
                ::Layer{ 5000,  3430, 3200, 3500},
                ::Layer{ 8000,  3567, 3500, 3650},
                ::Layer{12000,  3690, 3650, 3700},
                ::Layer{16000,  3720, 3700, 3900},
                ::Layer{21000,  3975, 3900, 4400},
                ::Layer{50000,  4950, 4400, 5500}
            };
        }
    }
#ifndef NDEBUG
    assert(!layers.empty());
#endif

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
#ifdef WITH_UMPS
        logger->error(e.what());
#else
        std::cerr << e.what() << std::endl;
#endif
    }
#ifdef WITH_UMPS
    logger->info("Final loss: " + std::to_string(loss));
#else
    std::cout << "Final loss: " << loss << std::endl;
#endif
    for (int i = 0; i < objectiveFunction.interfaces.size(); ++i)
    {
        std::cout << objectiveFunction.interfaces[i] << " "
                  << x[i] << std::endl;
    }
    objectiveFunction.nEvaluations = 0;
    auto fFinal = objectiveFunction.f(x.size(), x.data());
#ifdef WITH_UMPS
    logger->info("Evaluating final loss: " + std::to_string(fFinal)
               + " ; " 
               + std::to_string( (f0 - fFinal)/f0*100 ) + " pct reduction");
#else
    std::cout << "Evaluating final loss: " << fFinal << "; "
              << (f0 - fFinal)/f0*100 <<  " pct reduction" << std::endl;
#endif 
    objectiveFunction.writeToCSV(finalResidualsFile);
    return EXIT_SUCCESS;
}
