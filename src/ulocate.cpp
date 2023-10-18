#include <fstream>
#include <iostream>
#include <string>
#include <filesystem>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/tokenizer.hpp>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
//#include "uLocator/direct.hpp"
//#include "uLocator/directOptions.hpp"
#include "uLocator/uussTravelTimeCalculator.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
//#include "uLocator/hamiltonianMonteCarloOptions.hpp"
//#include "uLocator/hamiltonianMonteCarlo.hpp"
#include "uLocator/nlopt.hpp"
#include "uLocator/nloptOptions.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/topography.hpp"

#define UTM_ZONE 12

using namespace ULocator;

// 
//0    1       2       3        4        5     6          7            8            9             10            11                       12                      13                   14            15          16            17        18        19          20          21        22             23    24 
//evid,network,station,location,channelz,phase,arrival_id,arrival_time,pick_quality,first_motion,take_off_angle,source_receiver_distance,source_receiver_azimuth,travel_time_residual,receiver_lat,receiver_lon,receiver_elev,event_lat,event_lon,event_depth,origin_time,magnitude,magnitude_type,rflag,etype
//60000004,UU,HLJ,01,EHZ,P,228,1349658396.9360406,1.0,1,140,10.8,296.0,0.03,40.6105,-111.40067,1931.0,40.5678333,-111.2855,12.41,1349658393.5900002,0.01,d,F,eq

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

std::vector<Origin> loadCatalog(const std::string &fileName,
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
    std::vector<Origin> origins;
    origins.reserve(10000);
    std::ifstream csvFile(fileName); 
    if (!csvFile.is_open())
    {
        throw std::runtime_error("Couldn't open: " + fileName);
    }
    Origin origin;
    std::vector<Arrival> arrivals;
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
            Position::WGS84 epicenter{eventLat, eventLon, UTM_ZONE};
            origin.setIdentifier(evid);
            origin.setEpicenter(epicenter);
            origin.setDepth(eventDepth);
            origin.setTime(originTime);
            origin.setEventType(Origin::EventType::Earthquake);
            if (eventType == "qb")
            {
                origin.setEventType(Origin::EventType::QuarryBlast);
            }
        }
        Position::WGS84 stationPosition{stationLat, stationLon, UTM_ZONE};
        Arrival arrival;
        Station station;
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

std::vector<Station> uniqueStationsFromOrigins(
    const std::vector<Origin> &origins)
{
    std::vector<Station> uniqueStations;
    std::vector<std::string> uniqueStationNames;
    for (const auto &origin : origins)
    {
        const auto arrivals = origin.getArrivalsReference();
        for (const auto &arrival : arrivals)
        {
            const auto station = arrival.getStationReference();
            auto name = station.getNetwork() + "." + station.getName();
            bool exists{false};
            for (const auto &uniqueStationName : uniqueStationNames)
            {
                if (name == uniqueStationName)
                {
                    exists = true;
                    break;
                }
            }
            if (!exists)
            {
                uniqueStationNames.push_back(name);
                uniqueStations.push_back(station);
            }
        }
    }
    return uniqueStations;
}

/*
std::vector<ULocator::Arrival> 
    loadJSON(const std::filesystem::path &jsonFileName)
{
    std::ifstream jsonFile(jsonFileName);
    auto obj = nlohmann::json::parse(jsonFile);
    jsonFile.close();
    std::vector<ULocator::Arrival> arrivals;
    for (const auto &item : obj)
    {
        auto stationLatitude  = item["station_latitude"].get<double> ();
        auto stationLongitude = item["station_longitude"].get<double> ();
        ULocator::Position::WGS84 position{stationLatitude, stationLongitude, UTM_ZONE};
        ULocator::Station station;
        station.setNetwork(item["network"].get<std::string> ());
        station.setName(item["station"].get<std::string> ());
        station.setGeographicPosition(position);
        ULocator::Arrival arrival;
        arrival.setStation(station);
        arrival.setIdentifier(item["arrival_id"].get<int64_t> ());
        arrival.setTime(item["arrival_time"].get<double> ()); 
        arrival.setStandardError(item["standard_error"].get<double> ());
        auto phase = item["phase"].get<std::string> ();
        if (phase == "P")
        {
            arrival.setPhase(ULocator::Arrival::PhaseType::P);
        }
        else if (phase == "S")
        {
            arrival.setPhase(ULocator::Arrival::PhaseType::S);
        }
        else
        {
            assert(false);
            std::cout << "Unhandled phase type: " << phase << std::endl;
            continue;
        }
        arrivals.push_back(arrival);
    }
    return arrivals;
}
*/

std::shared_ptr<UMPS::Logging::ILog> makeLogger()
{
    return std::make_shared<UMPS::Logging::StandardOut> ();
}

/// @brief Parses the command line options.
struct ProgramOptions
{
    std::string catalogFile;
    std::string eventsFile;
    std::string outputFile;
    std::string travelTimeFile;
    std::string topographyFile;
    std::string correctionsFile;
    std::string region{"utah"};
    bool doSourceSpecificStationCorrections{true};
    bool doStaticCorrections{true};
    bool locate{true};
};

[[nodiscard]] ::ProgramOptions parseCommandLineOptions(int argc, char *argv[])
{
    ::ProgramOptions options;
    boost::program_options::options_description desc(
R"""(
The ulocate utility locates all evens in  catalog.  Example usage:
    ulocate --catalog_file=utah_catalog.csv --events_file=utah_events_training.csv  --output_file=relocatedUtahCatalog.csv
Allowed options)""");
    desc.add_options()
        ("help",         "Produces this help message")
        ("catalog_file", boost::program_options::value<std::string> ()->default_value("../examples/utah/utah_catalog.csv"),
                         "A CSV file with the picks for each event")
        ("events_file",  boost::program_options::value<std::string> ()->default_value("../examples/utah/utah_events_training.csv"),
                         "A CSV with the events")
        ("output_file",  boost::program_options::value<std::string> ()->default_value("relocatedUtahCatalog.csv"),
                         "The output catalog file")
        ("traveltime_file", boost::program_options::value<std::string> ()->default_value("../examples/utah/utahTravelTimes.h5"),
                         "The archive with the travel time grids")
        ("corrections_file", boost::program_options::value<std::string> ()->default_value(""), //correctionsArchive.h5"),
                         "The travel time corrections file")
        ("topography_file",  boost::program_options::value<std::string> ()->default_value("utahTopo.h5"),
                         "The topography file")
        ("region", boost::program_options::value<std::string> ()->default_value("utah"),
                   "The region - e.g., utah or ynp")
        ("disable_static_corrections", "If present and the corrections file is specified then the static corrections will not be applied")
        ("disable_source_specific_station_corrections", "If present and the corrections file is specified then the source specific station corrections will not be applied")
        ("predict", "If present then we will only predict travel times.");
    boost::program_options::variables_map vm;
    boost::program_options::store(
        boost::program_options::parse_command_line(argc, argv, desc), vm); 
    boost::program_options::notify(vm);
    if (vm.count("help"))
    {    
        std::cout << desc << std::endl;
        return options;
    }    
    if (vm.count("catalog_file"))
    {
        options.catalogFile = vm["catalog_file"].as<std::string> ();
        if (!std::filesystem::exists(options.catalogFile))
        {
            throw std::runtime_error("Catalog file: " + options.catalogFile
                                   + " does not exist");
        }
    }
    else 
    {    
        throw std::runtime_error("Catalog file was not set");
    }
    if (vm.count("events_file"))
    {
        options.eventsFile = vm["events_file"].as<std::string> ();
        if (!std::filesystem::exists(options.eventsFile))
        {
            throw std::runtime_error("Events file: " + options.eventsFile
                                   + " does not exist");
        }
    }
    else
    {
        throw std::runtime_error("Events file was not set");
    }
    if (vm.count("output_file"))
    {
        options.outputFile = vm["output_file"].as<std::string> ();
        if (std::filesystem::exists(options.outputFile))
        {
            std::cerr << "Output file: " << options.outputFile
                      << " will be overwritten" << std::endl;
        }
    }
    else
    {
        throw std::runtime_error("Output file was not set");
    }
    if (vm.count("traveltime_file"))
    {
        options.travelTimeFile = vm["traveltime_file"].as<std::string> ();
        if (!std::filesystem::exists(options.travelTimeFile))
        {
            throw std::runtime_error("Travel time file: " + 
                                     options.travelTimeFile
                                   + " does not exist");
        } 
    }
    else
    {
        throw std::runtime_error("Travel time file was not set");
    }
    if (vm.count("topography_file"))
    {
        auto topographyFile = vm["topography_file"].as<std::string> ();
        if (std::filesystem::exists(topographyFile))
        {
            options.topographyFile = topographyFile;
        }
    }
    else
    {
        std::cout << "Will not use topography" << std::endl;
    }
    if (vm.count("corrections_file"))
    {
        auto correctionsFile = vm["corrections_file"].as<std::string> ();
        if (std::filesystem::exists(correctionsFile))
        {
            options.correctionsFile = correctionsFile;
        }
        else
        {
            if (!correctionsFile.empty())
            {
                throw std::invalid_argument(correctionsFile + " does not exist");
            }
        }
    }
    else
    {
        std::cout << "Will not use corrections" << std::endl;
    }
    if (vm.count("region"))
    {
        auto region = vm["region"].as<std::string> ();
        if (region != "utah" && region != "ynp")
        {
            std::cerr << "Unhandled region: " << region << std::endl;
        }
        options.region = region;
    }
    else
    {
        throw std::invalid_argument("Region not set");
    }
    options.locate = true;
    if (vm.count("predict"))
    {
        options.locate = false;
    }
    options.doStaticCorrections = true;
    if (vm.count("disable_static_corrections"))
    {
        options.doStaticCorrections = false;
    }
    options.doSourceSpecificStationCorrections = true;
    if (vm.count("disable_source_specific_station_corrections"))
    {
        options.doSourceSpecificStationCorrections = false;
    }
    return options;
}

int main(int argc, char *argv[])
{
    ::ProgramOptions programOptions;
    try
    {
        programOptions = parseCommandLineOptions(argc, argv);
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    if (programOptions.catalogFile.empty()){return EXIT_SUCCESS;}
 
    auto logger = makeLogger();
//60086357,41.9576667,-112.8155,8.3,1411640115.6399994,0.2047364678717393,1.99,l,eq
//    std::filesystem::path jsonFileName{"60086357.json"}; 
    //std::filesystem::path jsonFileName{"../examples/utah/60000071.json"};
    //std::filesystem::path jsonFileName{"../examples/utah/60528932.json"};
    std::filesystem::path travelTimeFile = programOptions.travelTimeFile; //{"../examples/utah/utahTravelTimes.h5"};
    std::filesystem::path topographyFile = programOptions.topographyFile; // = {"utahTopo.h5"};
    std::filesystem::path correctionsFile = programOptions.correctionsFile; //{"correctionsArchive.h5"};

    logger->info("Loading origins...");
    auto origins = ::loadCatalog(programOptions.catalogFile,
                                 programOptions.eventsFile,
                                 logger);
    logger->info("Tabulating unique stations...");
    auto uniqueStations = ::uniqueStationsFromOrigins(origins); 
    logger->info("Found: " + std::to_string(uniqueStations.size()) + " stations");
    // Load the topography
    logger->info("Loading topography...");
    auto topography = std::make_unique<ULocator::Topography> (); 
    try 
    {
        topography->load(topographyFile);
    }
    catch (const std::exception &e) 
    {
        logger->error(e.what());
        return EXIT_FAILURE;
    }
    // Load the travel time tables
    std::vector<Station> stationShitList;
    auto travelTimeCalculators = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    for (const auto &uniqueStation : uniqueStations)
    {
        auto pCalculator = std::make_unique<ULocator::UUSSTravelTimeCalculator> (logger);
        auto sCalculator = std::make_unique<ULocator::UUSSTravelTimeCalculator> (logger);
        try
        {
            pCalculator->load(travelTimeFile,
                              uniqueStation,
                              "P",
                              programOptions.region);
            if (!correctionsFile.empty())
            {
                if (programOptions.doStaticCorrections)
                {
                    logger->info("Loading P static corrections from: "
                               + std::string {correctionsFile});
                    pCalculator->loadStaticCorrection(correctionsFile);
                }
                if (programOptions.doSourceSpecificStationCorrections)
                {
                    logger->info(
                        "Loading P source specific station corrections from: "
                       + std::string {correctionsFile});
                    pCalculator->loadSourceSpecificStationCorrections(
                       correctionsFile);
                }
            }
            sCalculator->load(travelTimeFile,
                              uniqueStation,
                              "S",
                              programOptions.region);
            if (!correctionsFile.empty())
            {
                if (programOptions.doStaticCorrections)
                {
                    logger->info("Loading S static corrections from: "
                               + std::string {correctionsFile});
                    sCalculator->loadStaticCorrection(correctionsFile);
                }
                if (programOptions.doSourceSpecificStationCorrections)
                {
                    logger->info(
                        "Loading S source specific station corrections from: "
                       + std::string {correctionsFile});
                    sCalculator->loadSourceSpecificStationCorrections(
                       correctionsFile);
                }
            }
        }
        catch (const std::exception &e)
        {
            bool exists{false};
            for (const auto &station : stationShitList)
            {
                if (station.getNetwork() == uniqueStation.getNetwork() && 
                    station.getName() == uniqueStation.getName())
                {
                    exists = true;
                }
            }
            if (!exists){stationShitList.push_back(uniqueStation);}
            std::cerr << e.what() << std::endl;
            continue;
        }
        travelTimeCalculators->insert(uniqueStation, "P", std::move(pCalculator));
        travelTimeCalculators->insert(uniqueStation, "S", std::move(sCalculator));
    }
    auto region = NLOptOptions::Region::Utah;
    if (programOptions.region == "ynp")
    {
        region = NLOptOptions::Region::Yellowstone;
        logger->info("Using the Yellowstone region"); 
    }
    else
    {
        logger->info("Using the Utah region");
    }
    bool locate = true;
    if (programOptions.locate)
    {
        locate = true;
        logger->info("Will perform location");
    }
    else
    {
        locate = false;
        logger->info("Will compute predictions for given events");
    }
/*
    ULocator::DirectOptions options(region);
    options.setObjectiveFunction(DirectOptions::ObjectiveFunction::L1);
    options.setAbsoluteModelTolerance(1.e-6);
*/
    ULocator::NLOptOptions options(region);
    options.setObjectiveFunction(NLOptOptions::ObjectiveFunction::L1);
    ULocator::NLOpt nloptSolver(logger);
    nloptSolver.setOptions(options);
    nloptSolver.setTopography(std::move(topography));
    nloptSolver.setTravelTimeCalculatorMap(std::move(travelTimeCalculators));
/*
    ULocator::Direct directSolver(logger);
    directSolver.setOptions(options);
    directSolver.setTopography(std::move(topography));
    directSolver.setTravelTimeCalculatorMap(std::move(travelTimeCalculators));
*/
    std::ofstream outCatalog(programOptions.outputFile);
    std::ofstream locationFailures("locationFailures.csv");
    outCatalog << "event_identifier,latitude,longitude,depth,origin_time,network,station,phase,arrival_time,standard_error,residual,uncorrected_travel_time,source_receiver_distance,event_type,n_objective_function_evaluations,weightedRMS" << std::endl;
    locationFailures << "event_identifier,latitude,longitude,depth,event_type,reason" << std::endl; 
    for (int iOrigin = 0; iOrigin < static_cast<int> (origins.size()); ++iOrigin)
    {
        logger->info("Processing " + std::to_string(iOrigin + 1) + " out of "
                   + std::to_string(origins.size()));
        auto origin = origins[iOrigin];
//if (origin.getIdentifier() != 60086357){continue;}
//if (origin.getIdentifier() != 60029937){continue;}
        auto arrivals = origin.getArrivals();
        arrivals.erase(
            std::remove_if(arrivals.begin(), arrivals.end(),
                           [&](const Arrival &arrival)
                      {
                          auto arrivalStation = arrival.getStationReference();
                          for (const auto &station : stationShitList)
                          {
                              if (arrivalStation.getNetwork() ==
                                  station.getNetwork() &&
                                  arrivalStation.getName() == station.getName())
                              {
                                  logger->warn("Purging: " + station.getNetwork()
                                             + "." + station.getName());
                                  return true;
                              }
                          }
                          return false;
                      }),
             arrivals.end());
        if (arrivals.size() < 4 && locate)
        {
            logger->error("Cannot locate event: "
                        + std::to_string(origin.getIdentifier())
                        + " because too few observations");
            continue;
        }
        try
        {
            nloptSolver.setArrivals(arrivals);
        }
        catch (const std::exception &e)
        {
            logger->error("Could not set arrivals for: "
                        + std::to_string(origin.getIdentifier())
                        + ".  Failed with: " + std::string {e.what()});
            continue;
        }
        try
        {
            bool isQuarryBlast = false;
            if (origin.getEventType() == Origin::EventType::QuarryBlast)
            {
                isQuarryBlast = true;
            }
            int nObjectiveFunctionEvaluations = 1;
            if (locate)
            {
                if (isQuarryBlast)
                {
                    logger->info("Locating quarry blast: "
                               + std::to_string(origin.getIdentifier()));
                    nloptSolver.locateQuarryBlast();
                    //(Direct::SourceDepthConstraint::FixedToFreeSurface);
                }
                else
                {
                    logger->info("Locating event: "
                               + std::to_string(origin.getIdentifier()));
                    nloptSolver.locateEarthquake(); //directSolver.locate();
                }
                nObjectiveFunctionEvaluations
                    = nloptSolver.getNumberOfObjectiveFunctionEvaluations()
                    + nloptSolver.getNumberOfGradientEvaluations();

            }
            else
            {
                if (isQuarryBlast)
                {
                    logger->info("Computing travel times for quarry blast: "
                               + std::to_string(origin.getIdentifier()));
                }
                else
                {
                    logger->info("Computing travel times for event: "
                               + std::to_string(origin.getIdentifier()));
                }
            }
            // Tabulate the raw travel times
            Origin newOrigin;
            if (locate)
            {
                newOrigin = nloptSolver.getOrigin();
            }
            else
            {
                newOrigin = origin;
            }
            constexpr bool applyCorrection{false};
            newOrigin = nloptSolver.predict(newOrigin, applyCorrection);
            auto newArrivals = newOrigin.getArrivals();
            auto newEpicenter = newOrigin.getEpicenter();
            auto newDepth = newOrigin.getDepth();
            double numerator = 0;
            double denominator = 0;
            std::vector<double> residuals;
            for (const auto &newArrival : newArrivals)
            {
                auto weight = 1./std::max(1.e-14, newArrival.getStandardError());
                numerator = numerator + weight*std::pow(newArrival.getResidual(), 2);
                denominator = denominator + weight;
                residuals.push_back(newArrival.getResidual());
            }
            double weightedRMS = std::sqrt(numerator/std::max(1.e-14, denominator));
            
/*
            std::vector<double> travelTimes;
            travelTimes.resize(newArrivals.size(), 0);
            auto travelTimeMap = directSolver.releaseTravelTimeCalculatorMap();
            int iArrival = 0;
            for (const auto &newArrival : newArrivals)
            {
                constexpr bool applyCorrection{false};
                travelTimes.at(iArrival) = 
                    travelTimeMap->evaluate(newArrival.getStation(), newArrival.getPhase(),
                                            newEpicenter, newDepth,
                                            applyCorrection); 
                iArrival = iArrival + 1;
            }
            directSolver.setTravelTimeCalculatorMap(std::move(travelTimeMap));
*/
            // Write it all out
            std::string eventType = "eq";
            if (isQuarryBlast){eventType = "qb";}
            int iArrival = 0;
            for (const auto &newArrival : newArrivals)
            {
                // r = obs - est = obs - (tt + ot) = obs - tt - ot
                // tt = obs - r - ot
                double travelTime = newArrival.getTime()
                                  - newArrival.getResidual()
                                  - newOrigin.getTime();
                outCatalog << std::setprecision(16) 
                           << origin.getIdentifier() << "," 
                           << newEpicenter.getLatitude() << ","
                           << newEpicenter.getLongitude() << ","
                           << newDepth << ","
                           << newOrigin.getTime() << ","
                           << newArrival.getStationReference().getNetwork() << ","
                           << newArrival.getStationReference().getName() << ","
                           << newArrival.getPhase() << ","
                           << newArrival.getTime() << ","
                           << newArrival.getStandardError() << ","
                           << residuals[iArrival] << "," //newArrival.getResidual() << ","
                           << travelTime << ","
                           << newArrival.getDistance() << ","
                           << eventType << ","
                           << nObjectiveFunctionEvaluations << ","
                           << weightedRMS << std::endl;
                iArrival = iArrival + 1;
            }
            //for (const auto &a : newArrivals){std::cout << a.getResidual() << std::endl;}
        }
        catch (const std::exception &e)
        {
            logger->error("Failed to locate: "
                        + std::to_string(origin.getIdentifier())
                        + " Failed with: " + std::string {e.what()});
            auto eventType = "eq";
            if (origin.getEventType() == Origin::EventType::QuarryBlast)
            {
                eventType = "qb";
            }
            locationFailures << origin.getIdentifier() << ","
                             << origin.getEpicenter().getLatitude() << ","
                             << origin.getEpicenter().getLongitude() << ","
                             << origin.getDepth() << ","
                             << eventType << ","
                             << "\"" << e.what() << "\""
                             << std::endl;
            continue;
        }
    } 
    outCatalog.close();
    locationFailures.close();
    return EXIT_SUCCESS;
}
