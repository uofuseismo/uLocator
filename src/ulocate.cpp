#include <fstream>
#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include <filesystem>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/tokenizer.hpp>
#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/LocalCartesian.hpp>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/optimizers/nlopt/dividedRectangles.hpp"
#include "uLocator/optimizers/pagmo/particleSwarm.hpp"
#include "uLocator/optimizers/nlopt/boundOptimizationByQuadraticApproximation.hpp"
#include "uLocator/optimizers/prima/boundOptimizationByQuadraticApproximation.hpp"
//#include "uLocator/direct.hpp"
//#include "uLocator/directOptions.hpp"
//#include "uLocator/uussTravelTimeCalculator.hpp"
#include "uLocator/uussRayTracer.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/topography/constant.hpp"
#include "uLocator/topography/gridded.hpp"
#include "uLocator/corrections/static.hpp"
#include "uLocator/corrections/sourceSpecific.hpp"
#include "loadCatalog.hpp"
#include "utahQuarries.hpp"
#include "utahQuarryBlastTravelTimeDatabase.hpp"
#include "originTime.hpp"
#include "searchStations.hpp"

using namespace ULocator;

// 
//0    1       2       3        4        5     6          7            8            9             10            11                       12                      13                   14            15          16            17        18        19          20          21        22             23    24 
//evid,network,station,location,channelz,phase,arrival_id,arrival_time,pick_quality,first_motion,take_off_angle,source_receiver_distance,source_receiver_azimuth,travel_time_residual,receiver_lat,receiver_lon,receiver_elev,event_lat,event_lon,event_depth,origin_time,magnitude,magnitude_type,rflag,etype
//60000004,UU,HLJ,01,EHZ,P,228,1349658396.9360406,1.0,1,140,10.8,296.0,0.03,40.6105,-111.40067,1931.0,40.5678333,-111.2855,12.41,1349658393.5900002,0.01,d,F,eq
//0,              1,       2,     3,                4             5
//event_identifier,network,station,vertical_channel,location_code,arrival_time,residual,phase,station_latitude,station_longitude,station_elevation,standard_error,arrival_identifier,arrival_evaluation_type,event_latitude,event_longitude,event_depth,origin_time,source_receiver_distance_km,source_to_receiver_azimuth,receiver_to_source_azimuth,event_type

/// @result An application logger.
std::shared_ptr<UMPS::Logging::ILog> makeLogger()
{
    return std::make_shared<UMPS::Logging::StandardOut> (UMPS::Logging::Level::Info);
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
    std::unique_ptr<Position::IGeographicRegion> geographicRegion{nullptr};
    int catalogVersion{3};
    double utahDefaultDepth{4780};
    double ynpDefaultDepth{6600};
    double utahTimeWindow{140};  // Longest seen is 80 s (50 pct extra)
    double ynpTimeWindow{61}; // Longest seen is 42 seconds (50 pct extra)
    double pNorm{1.5};
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
        ("catalog_file", boost::program_options::value<std::string> ()->default_value("../examples/uuss/utah_catalog.csv"),
                         "A CSV file with the picks for each event")
        ("catalog_version", boost::program_options::value<int> ()->default_value(options.catalogVersion),
                         "The catalog file version number")
        ("events_file",  boost::program_options::value<std::string> ()->default_value("../examples/uuss/utah_events_training.csv"),
                         "A CSV with the events")
        ("output_file",  boost::program_options::value<std::string> ()->default_value("relocatedUtahCatalog.csv"),
                         "The output catalog file")
        ("traveltime_file", boost::program_options::value<std::string> ()->default_value("../examples/uuss/utahTravelTimes.h5"),
                         "The archive with the travel time grids")
        ("corrections_file", boost::program_options::value<std::string> ()->default_value(""), //correctionsArchive.h5"),
                         "The travel time corrections file")
        ("topography_file",  boost::program_options::value<std::string> ()->default_value("utahTopography.h5"),
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
        std::cout << "Will use constant topography" << std::endl;
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
            throw std::invalid_argument("Unhandled region: " + region);
        }
        else
        {
            if (region == "utah")
            {
                options.geographicRegion = std::make_unique<Position::UtahRegion> ();
            }
            else if (region == "ynp")
            {
                options.geographicRegion = std::make_unique<Position::YNPRegion> ();
            }
        } 
        options.region = region;
    }
    else
    {
        throw std::invalid_argument("Region not set");
    }
    if (vm.count("catalog_version"))
    {
        auto version = vm["catalog_version"].as<int> ();
        if (version < 1 || version > 3)
        {
            throw std::invalid_argument("catalog_version must be 1, 2, or 3");
        }
        options.catalogVersion = version;
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
    if (options.doStaticCorrections &&
        !std::filesystem::exists(options.correctionsFile))
    {
        options.doStaticCorrections = false;
    }
    options.doSourceSpecificStationCorrections = true;
    if (vm.count("disable_source_specific_station_corrections"))
    {
        options.doSourceSpecificStationCorrections = false;
    }
    if (options.doSourceSpecificStationCorrections &&
        !std::filesystem::exists(options.correctionsFile))
    {
        options.doSourceSpecificStationCorrections = false;
    }
    return options;
}

/*
class RefinedRegion : public ULocator::Position::IGeographicRegion
{
public:
    ~RefinedRegion() = default;
    RefinedRegion(const ULocator::Position::WGS84 &epicenter,
                  const ULocator::Position::IGeographicRegion &region,
                  const double refinedSearchInX = 50000,
                  const double refinedSearchInY = 50000) :
        IGeographicRegion()
    {
        auto latitude = epicenter.getLatitude();
        auto longitude = epicenter.getLongitude();
        GeographicLib::Geocentric mEarth{GeographicLib::Constants::WGS84_a(),
                                         GeographicLib::Constants::WGS84_f()};
        GeographicLib::LocalCartesian mProjection{latitude,
                                                  longitude,
                                                  0.0, // Height above ellipsoid at origin (meters)
                                                  mEarth};
        mMinimumX =-std::abs(refinedSearchInX);
        mMaximumX = std::abs(refinedSearchInX);
        mMinimumY =-std::abs(refinedSearchInY);
        mMaximumY = std::abs(refinedSearchInY);
    }
    /// Forward transformation
    std::pair<double, double> 
        geographicToLocalCoordinates(const double latitude,
                                     const double longitude) const override
    {
        if (latitude < -90 || latitude > 90) 
        {   
            throw std::invalid_argument("Latitude " + std::to_string(latitude)
                                      + " must be in range [-90,90]");
        }   
        double x, y, z;
        mProjection.Forward(latitude, longitude, 0.0, x, y, z); 
        return std::pair {x, y}; 
    }
    /// Reverse transformation
    std::pair<double, double> 
        localToGeographicCoordinates(const double x, const double y) const override
    {
        double latitude, longitude, h;
        mProjection.Reverse(x, y, 0.0, latitude, longitude, h); 
        return std::pair {latitude, longitude}; 
    }
    /// X-extent
    [[nodiscard]] std::pair<double, double> getExtentInX() const noexcept override
    {   
        return std::pair {mMinimumX, mMaximumX};
    }
    /// Y-extent
    std::pair<double, double> getExtentInY() const noexcept override
    {
        return std::pair {mMinimumY, mMaximumY};
    }
    GeographicLib::LocalCartesian mProjection;
    double mMinimumX{0};
    double mMaximumX{0};
    double mMinimumY{0};
    double mMaximumY{0};
};
*/

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
 
    auto logger = ::makeLogger();
    
/*
//60086357,41.9576667,-112.8155,8.3,1411640115.6399994,0.2047364678717393,1.99,l,eq
//    std::filesystem::path jsonFileName{"60086357.json"}; 
    //std::filesystem::path jsonFileName{"../examples/utah/60000071.json"};
    //std::filesystem::path jsonFileName{"../examples/utah/60528932.json"};
    std::filesystem::path travelTimeFile = programOptions.travelTimeFile; //{"../examples/utah/utahTravelTimes.h5"};
    std::filesystem::path topographyFile = programOptions.topographyFile; // = {"utahTopo.h5"};
    std::filesystem::path correctionsFile = programOptions.correctionsFile; //{"correctionsArchive.h5"};
*/

    // Load up the catalog to (re)locate, station information, etc. 
    logger->info("Loading origins...");
    constexpr int utmZone{12};
    auto origins = ::loadCatalog(programOptions.catalogFile,
                                 programOptions.eventsFile,
                                 *programOptions.geographicRegion,
                                 logger,
                                 programOptions.catalogVersion,
                                 utmZone);
    logger->info("Tabulating unique stations...");
    auto uniqueStations = ::uniqueStationsFromOrigins(origins); 
    logger->info("Found: "
               + std::to_string(uniqueStations.size()) + " stations");

    // Load the topography
    std::unique_ptr<ULocator::Topography::ITopography> topography{nullptr};
    if (std::filesystem::exists(programOptions.topographyFile))
    {
        logger->info("Loading topography from "
                   + programOptions.topographyFile);
        auto griddedTopography
            = std::make_unique<ULocator::Topography::Gridded> ();
        try 
        {
            griddedTopography->load(programOptions.topographyFile,
                                    *programOptions.geographicRegion);
        }
        catch (const std::exception &e) 
        {
            logger->error(e.what());
            return EXIT_FAILURE;
        }
        topography = std::move(griddedTopography);
    }
    else
    {
        logger->info("Using constant topography of 2 km");
        auto constantTopography
            = std::make_unique<ULocator::Topography::Constant> ();
        constantTopography->set(2000);
        topography = std::move(constantTopography);
    }

    // Load the quarries
    std::unique_ptr<::UtahQuarryBlastTravelTimeDatabase> quarryLocator{nullptr};
    if (programOptions.locate && programOptions.region == "utah")
    {
        quarryLocator = std::make_unique<::UtahQuarryBlastTravelTimeDatabase> ();
        if (programOptions.doStaticCorrections)
        {
            quarryLocator->setStaticCorrectionsFile(
                programOptions.correctionsFile);
        }
        if (programOptions.doSourceSpecificStationCorrections)
        {
            quarryLocator->setSourceSpecificCorrectionsFile(
                programOptions.correctionsFile);
        }
    }
    // Load the known locations
    std::unique_ptr<::UtahEventTravelTimeDatabase> utahEventDatabase{nullptr};
    std::unique_ptr<::YNPEventTravelTimeDatabase> ynpEventDatabase{nullptr};
    if (programOptions.locate)
    {
        if (programOptions.region == "utah")
        {
            logger->info("Making Utah event travel time database...");
            utahEventDatabase = std::make_unique<::UtahEventTravelTimeDatabase> ();
            if (programOptions.doStaticCorrections)
            {
                utahEventDatabase->setStaticCorrectionsFile(
                    programOptions.correctionsFile);
            }
            if (programOptions.doSourceSpecificStationCorrections)
            {
                utahEventDatabase->setSourceSpecificCorrectionsFile(
                    programOptions.correctionsFile);
            }
        }
        else
        {
            logger->info("Making YNP event travel time database...");
            ynpEventDatabase = std::make_unique<::YNPEventTravelTimeDatabase> ();
            if (programOptions.doStaticCorrections)
            {
                ynpEventDatabase->setStaticCorrectionsFile(
                    programOptions.correctionsFile);
            }
            if (programOptions.doSourceSpecificStationCorrections)
            {
                ynpEventDatabase->setSourceSpecificCorrectionsFile(
                    programOptions.correctionsFile);
            }
        }
    }

    // Generate the ray-tracer travel time calculators
    logger->info("Initializing travel time calculators...");
    std::vector<Station> stationBlackList;
    auto travelTimeCalculators = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    for (const auto &uniqueStation : uniqueStations)
    {
        auto rayTracerRegion = UUSSRayTracer::Region::Utah;
        if (programOptions.region == "ynp")
        {
            rayTracerRegion = UUSSRayTracer::Region::YNP;
        }
        try
        {
            ULocator::Corrections::Static pStatic;
            ULocator::Corrections::Static sStatic;
            ULocator::Corrections::SourceSpecific pSSSC;
            ULocator::Corrections::SourceSpecific sSSSC;
            if (programOptions.doStaticCorrections)
            {
                pStatic.setStationNameAndPhase(uniqueStation.getNetwork(),
                                               uniqueStation.getName(),
                                               "P");
                sStatic.setStationNameAndPhase(uniqueStation.getNetwork(),
                                               uniqueStation.getName(),
                                               "S");
                try
                {
                    pStatic.load(programOptions.correctionsFile);
                    logger->info("Loaded P static for "
                               + uniqueStation.getNetwork()
                               + " " + uniqueStation.getName());
                }
                catch (const std::exception &e)
                {
                    logger->warn(e.what());
                    pStatic.clear();
                }
                try
                {
                    sStatic.load(programOptions.correctionsFile);
                    logger->info("Loaded S static for "
                               + uniqueStation.getNetwork()
                               + " " + uniqueStation.getName());
                }
                catch (const std::exception &e)
                {
                    logger->warn(e.what());
                    sStatic.clear();
                }
            }
            if (programOptions.doSourceSpecificStationCorrections)
            {
                pSSSC.setStationNameAndPhase(uniqueStation.getNetwork(),
                                             uniqueStation.getName(),
                                             "P");
                sSSSC.setStationNameAndPhase(uniqueStation.getNetwork(),
                                             uniqueStation.getName(),
                                             "S");
                try
                {
                    pSSSC.load(programOptions.correctionsFile);
                    logger->info("Loaded P SSSC for "
                               + uniqueStation.getNetwork() 
                               + " " + uniqueStation.getName());
                }
                catch (const std::exception &e)
                {
                    logger->warn(e.what());
                    pSSSC.clear();
                }
                try
                {
                    sSSSC.load(programOptions.correctionsFile);
                    logger->info("Loaded S SSSC for " 
                               + uniqueStation.getNetwork()
                               + " " + uniqueStation.getName());
                }
                catch (const std::exception &e)
                {
                    logger->warn(e.what());
                    sSSSC.clear();
                }
            }
            auto pCalculator
                = std::make_unique<UUSSRayTracer> (uniqueStation,
                                                   UUSSRayTracer::Phase::P,
                                                   rayTracerRegion,
                                                   std::move(pStatic),
                                                   std::move(pSSSC),
                                                   logger);
            auto sCalculator
                = std::make_unique<UUSSRayTracer> (uniqueStation,
                                                   UUSSRayTracer::Phase::S,
                                                   rayTracerRegion,
                                                   std::move(sStatic),
                                                   std::move(sSSSC),
                                                   logger);
            travelTimeCalculators->insert(uniqueStation, "P",
                                          std::move(pCalculator));
            travelTimeCalculators->insert(uniqueStation, "S",
                                          std::move(sCalculator));
        }
        catch (const std::exception &e)
        {
            bool exists{false};
            for (const auto &station : stationBlackList)
            {
                if (station.getNetwork() == uniqueStation.getNetwork() && 
                    station.getName() == uniqueStation.getName())
                {
                    exists = true;
                }
            }
            if (!exists)
            {
                logger->warn("Blacklisting "
                           + uniqueStation.getNetwork() + "."
                           + uniqueStation.getName());
                stationBlackList.push_back(uniqueStation);
            }
            logger->error(e.what());
            continue;
        }
    }

    std::ofstream outCatalog(programOptions.outputFile);
    //std::ofstream locationFailures("locationFailures.csv");
    outCatalog << "event_identifier,latitude,longitude,depth,origin_time,network,station,phase,arrival_time,standard_error,residual,uncorrected_travel_time,source_receiver_distance,event_type,n_objective_function_evaluations,weightedRMSE" << std::endl;

    for (int iOrigin = 0; iOrigin < static_cast<int>(origins.size()); ++iOrigin)
    {
//if (iOrigin + 1 != 68){continue;}
//if (iOrigin > 20){break;}
        if ((iOrigin + 1)%50 == 0 || true)
        {
            double percentDone = (iOrigin*100.)/(origins.size() - 1);
            logger->info("Processing " + std::to_string(iOrigin + 1)
                       + " out of "
                       + std::to_string(origins.size())
                       + " (" + std::to_string(percentDone) + " %)");
        }
        auto origin = origins[iOrigin];
        auto arrivals = origin.getArrivals();
        bool isQuarryBlast = false;
        if (origin.getEventType() == Origin::EventType::QuarryBlast)
        {
            isQuarryBlast = true;
        }
        // If any arrivals are on a black-listed station then purge it 
        arrivals.erase(
            std::remove_if(arrivals.begin(), arrivals.end(),
                           [&](const Arrival &arrival)
                      {
                          auto arrivalStation = arrival.getStationReference();
                          for (const auto &station : stationBlackList)
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
        // Time to work
        int nObjectiveFunctionEvaluations{0};
        if (programOptions.locate)
        {
            // Define the origin time search window and fixed depth for the
            // initial optimization 
            double timeWindow = programOptions.utahTimeWindow;
            double initialDepth = programOptions.utahDefaultDepth;
            if (programOptions.region == "ynp")
            {
                timeWindow = programOptions.ynpTimeWindow; 
                initialDepth = programOptions.ynpDefaultDepth;
            }
            // First we build an initial solution.  Effectively, we will
            // always check quarries and sometimes check earthquake locations.
            double bestQuarryObjectiveFunction
                = std::numeric_limits<double>::max();
            ULocator::Origin bestQuarryOrigin;
            std::string bestQuarryName{"Undefined"};
            if (quarryLocator)
            {
                logger->debug("Performing quarry optimization...");
                //auto [bestQuarryIndex, bestQuarryOrigin,
                //      bestQuarryObjectiveFunction]
                auto bestQuarryLocatorResult
                    = quarryLocator->findBestQuarry(
                         origin,
                         ULocator::Optimizers::IOptimizer::Norm::Lp,
                         programOptions.pNorm,
                         timeWindow);
                auto bestQuarryIndex = std::get<0> (bestQuarryLocatorResult);
                bestQuarryOrigin = std::get<1> (bestQuarryLocatorResult);
                bestQuarryObjectiveFunction
                    = std::get<2> (bestQuarryLocatorResult);
                bestQuarryName = quarryLocator->at(bestQuarryIndex).getName();
            }
            ULocator::Origin bestInitialOrigin;
            double bestInitialObjectiveFunction{std::numeric_limits<double>::max()};
            ULocator::Optimizers::IOptimizer::LocationProblem problem;
            //std::vector<double> trialDepths;
            // For quarry blasts we fix the depth and call it a day
            if (origin.getEventType() ==
                ULocator::Origin::EventType::QuarryBlast)
            {
                problem = ULocator::Optimizers::IOptimizer::LocationProblem::FixedToFreeSurfaceAndTime;
                // Use a known quarry as an initial guess for a quarry solution
                ULocator::Origin initialGuess{bestQuarryOrigin};
                // Now optimize with DIRECT
                logger->info("Performing initial quarry blast search...");
                ULocator::Optimizers::NLOpt::DividedRectangles direct(logger);
                direct.setMaximumNumberOfObjectiveFunctionEvaluations(1200);
                direct.setLocationTolerance(1000);
                direct.setOriginTimeTolerance(1);
                direct.setOriginTimeSearchWindowDuration(timeWindow);
                direct.setTravelTimeCalculatorMap(std::move(travelTimeCalculators));
                direct.setTopography(std::move(topography));
                direct.setGeographicRegion(*programOptions.geographicRegion);
                direct.setArrivals(arrivals);
                direct.locate(initialGuess,
                              problem,
                              ULocator::Optimizers::IOptimizer::Norm::Lp);
                nObjectiveFunctionEvaluations = direct.getNumberOfObjectiveFunctionEvaluations();
                bestInitialOrigin = direct.getOrigin();
                bestInitialObjectiveFunction = direct.getOptimalObjectiveFunction();
                travelTimeCalculators = direct.releaseTravelTimeCalculatorMap();
                topography = direct.releaseTopography();
                if (bestQuarryObjectiveFunction < direct.getOptimalObjectiveFunction())
                {
                    logger->debug("Using quarry location "
                                + bestQuarryName
                                + " instead of DIRECT solution"); 
                    bestInitialOrigin = bestQuarryOrigin;
                    bestInitialObjectiveFunction = bestQuarryObjectiveFunction;
                }
                else
                {
                    logger->debug("Using DIRECT quarry solution");
                }
                // Locate for real
                logger->info("Fine-tuning quarry blast location...");
                //trialDepths.push_back(bestInitialOrigin.getDepth());
            }
            else
            {
                problem = ULocator::Optimizers::IOptimizer::LocationProblem::FixedDepthAndTime;
                // Assume the closest station (smallest arrival) time is where
                // the event is and locate based on that
                logger->debug("Performing station location optimization...");
                auto [bestStationOrigin, stationObjectiveFunction]
                    = ::searchStations(origin,
                                       *travelTimeCalculators,
                                       initialDepth,
                                       ULocator::Optimizers::IOptimizer::Norm::Lp,
                                       programOptions.pNorm,
                                       timeWindow,
                                       true);
                // Search the cluster centers 
                logger->debug("Searching earthquake cluster centroids...");
                double bestKnownEventObjectiveFunction{std::numeric_limits<double>::max()};
                ULocator::Origin bestKnownEventOrigin; 
                int bestKnownEventIndex;
                if (programOptions.region == "utah")
                {
                    auto bestEventInDatabase
                        = utahEventDatabase->findBestEvent(
                            origin,
                            ULocator::Optimizers::IOptimizer::Norm::Lp,
                            programOptions.pNorm,
                            timeWindow,
                            true);
                    bestKnownEventOrigin = std::get<1> (bestEventInDatabase);
                    bestKnownEventObjectiveFunction = std::get<2> (bestEventInDatabase);
                }
                else
                {
                    auto bestEventInDatabase
                        = ynpEventDatabase->findBestEvent(
                            origin,
                            ULocator::Optimizers::IOptimizer::Norm::Lp,
                            programOptions.pNorm,
                            timeWindow,
                            true);
                    bestKnownEventOrigin = std::get<1> (bestEventInDatabase);
                    bestKnownEventObjectiveFunction = std::get<2> (bestEventInDatabase);
                }

                // Do a crude (x,y,t) search with a fixed depth
                logger->info("Performing initial earthquake search...");
                ULocator::Optimizers::NLOpt::DividedRectangles direct(logger);
                direct.setMaximumNumberOfObjectiveFunctionEvaluations(1200);
                direct.setLocationTolerance(1000); // Doesn't matter
                direct.setOriginTimeTolerance(1); // Doesn't matter
                direct.setOriginTimeSearchWindowDuration(timeWindow);
                direct.setTravelTimeCalculatorMap(std::move(travelTimeCalculators));
                direct.setGeographicRegion(*programOptions.geographicRegion);
                direct.setTopography(std::move(topography));
                direct.setArrivals(arrivals);
                direct.locateEventWithFixedDepth(
                    initialDepth,
                    ULocator::Optimizers::IOptimizer::Norm::Lp);
                nObjectiveFunctionEvaluations
                    = direct.getNumberOfObjectiveFunctionEvaluations();
                bestInitialOrigin = direct.getOrigin();
                auto directOptimalObjectiveFunction = direct.getOptimalObjectiveFunction();
                travelTimeCalculators = direct.releaseTravelTimeCalculatorMap();
                topography = direct.releaseTopography();
                // Choose the best starting guess
                bestInitialObjectiveFunction = directOptimalObjectiveFunction;
                if (bestQuarryObjectiveFunction < directOptimalObjectiveFunction &&
                    bestQuarryObjectiveFunction < stationObjectiveFunction &&
                    bestQuarryObjectiveFunction < bestKnownEventObjectiveFunction)
                {
                    logger->debug("Using quarry location "
                                + bestQuarryName
                                + " for earthquake instead of DIRECT solution");
                    bestInitialOrigin = bestQuarryOrigin;
                    bestInitialObjectiveFunction = bestQuarryObjectiveFunction;
                }
                else if (stationObjectiveFunction < directOptimalObjectiveFunction &&
                         stationObjectiveFunction < bestQuarryObjectiveFunction &&
                         stationObjectiveFunction < bestKnownEventObjectiveFunction)
                {
                    logger->debug("Using closest station location for earthquake location instead of DIRECT solution");
                    bestInitialOrigin = bestStationOrigin;
                    bestInitialObjectiveFunction = stationObjectiveFunction;
                }
                else if (bestKnownEventObjectiveFunction < directOptimalObjectiveFunction &&
                         bestKnownEventObjectiveFunction < bestQuarryObjectiveFunction &&
                         bestKnownEventObjectiveFunction < stationObjectiveFunction)
                {
                    logger->debug("Using saved event location for earthquake instead of DIRECT solution");
                    bestInitialOrigin = bestKnownEventOrigin;
                    bestInitialObjectiveFunction = bestKnownEventObjectiveFunction;
                }
                // Locate for real
                logger->info("Fine-tuning earthquake location...");
/*
                //trialDepths.push_back(bestInitialOrigin.getDepth());
                if (programOptions.region == "utah")
                {
                    auto interfaces = ULocator::UUSSRayTracer::getInterfaces(
                        ULocator::UUSSRayTracer::Region::Utah);
                    for (int layer = 0;
                         layer < static_cast<int> (interfaces.size()) - 1;
                         ++layer)
                    {
                        auto averageDepth = 0.5*(interfaces[layer]
                                               + interfaces[layer + 1]);
                        trialDepths.push_back( std::max(-800., averageDepth) );
                    }
                }
                else
                {
                    auto interfaces = ULocator::UUSSRayTracer::getInterfaces(
                        ULocator::UUSSRayTracer::Region::YNP);
                    for (int layer = 0;
                         layer < static_cast<int> (interfaces.size()) - 1;
                         ++layer)
                    {
                        auto averageDepth = 0.5*(interfaces[layer]
                                               + interfaces[layer + 1]);
                        trialDepths.push_back( std::max(-1200., averageDepth) );
                    }
                }
*/
                // Define search depths
                problem = ULocator::Optimizers::IOptimizer::LocationProblem::ThreeDimensionsAndTime;
            }
            // Set some extra stuff 
            bestInitialOrigin.setIdentifier(origins[iOrigin].getIdentifier());
            bestInitialOrigin.setEventType(origins[iOrigin].getEventType());
            // Now locate for real 
            auto initialLatitude = bestInitialOrigin.getEpicenter().getLatitude();
            auto initialLongitude = bestInitialOrigin.getEpicenter().getLongitude();
            auto [initialX, initialY]
                = programOptions.geographicRegion->geographicToLocalCoordinates(
                      initialLatitude, initialLongitude);
            auto [x0, x1] = programOptions.geographicRegion->getExtentInX();
            auto [y0, y1] = programOptions.geographicRegion->getExtentInY();
            double refineX{50000};
            double refineY{50000};
            int nParticles{20};
            int nGenerations{140};
            if (programOptions.region == "ynp")
            {
                refineX = 35000;
                refineY = 35000;
                nParticles = 15;
                nGenerations = 120;
            }
            std::pair<double, double> newExtentInX {std::max(x0, initialX - refineX),
                                                    std::min(x1, initialX + refineX)};
            std::pair<double, double> newExtentInY {std::max(y0, initialY - refineY),
                                                    std::min(y1, initialY + refineY)};
            ULocator::Optimizers::Pagmo::ParticleSwarm pso(logger);
            pso.setTravelTimeCalculatorMap(std::move(travelTimeCalculators));
            pso.setTopography(std::move(topography));
            pso.setNumberOfParticles(nParticles);
            pso.setNumberOfGenerations(nGenerations);
            pso.setOriginTimeSearchWindowDuration(timeWindow);
            pso.setGeographicRegion(*programOptions.geographicRegion);
            pso.setExtentInX(newExtentInX);
            pso.setExtentInY(newExtentInY);
            pso.setArrivals(arrivals);
            auto optimalOrigin = bestInitialOrigin;
            double optimalDepthObjectiveFunction{bestInitialObjectiveFunction};
            try
            {
                pso.locate(bestInitialOrigin,
                           problem,
                           ULocator::Optimizers::IOptimizer::Norm::Lp);
                if (pso.getOptimalObjectiveFunction() < bestInitialObjectiveFunction)
                {
                    logger->debug("Using PSO solution");
                    optimalOrigin = pso.getOrigin();
                }
                else
                {
                    logger->warn("Using initial solution for "
                            + std::to_string(origins[iOrigin].getIdentifier()));
                }
            }
            catch (const std::exception &e)
            {
                 logger->error("PSO failed " + std::string {e.what()} + " for "
                            + std::to_string(origins[iOrigin].getIdentifier()));
            }
            travelTimeCalculators = pso.releaseTravelTimeCalculatorMap();
            topography = pso.releaseTopography();

/*
            ULocator::Optimizers::Prima::
                BoundOptimizationByQuadraticApproximation bobyqa(logger);
            bobyqa.setMaximumNumberOfObjectiveFunctionEvaluations(2500);
            bobyqa.setLocationTolerance(1.e-1); ////bobyqa.setLocationTolerance(1);
            bobyqa.setOriginTimeTolerance(1.e-3);
            bobyqa.setOriginTimeTolerance(0.001);
            bobyqa.setOriginTimeSearchWindowDuration(timeWindow);
            bobyqa.setTravelTimeCalculatorMap(std::move(travelTimeCalculators));
            bobyqa.setGeographicRegion(*programOptions.geographicRegion);
            bobyqa.setExtentInX(newExtentInX);
            bobyqa.setExtentInY(newExtentInY);
            bobyqa.setTopography(std::move(topography));
            bobyqa.setArrivals(arrivals);
            // Iterate through the velocity model and locate at different depths
            bool reducedObjectiveFunction{false};
            auto optimalOrigin = bestInitialOrigin;
            double optimalDepthObjectiveFunction{bestInitialObjectiveFunction};
std::cout << bestInitialObjectiveFunction << std::endl;
            for (const auto &trialDepth : trialDepths)
            {
//std::cout << trialDepth << std::endl;
                ULocator::Origin trialOrigin{optimalOrigin};
                trialOrigin.setDepth(trialDepth);
                try
                {
                    bobyqa.locate(trialOrigin,
                                  problem,
                                  ULocator::Optimizers::IOptimizer::Norm::Lp);
                }
                catch (const std::exception &e)
                {
                    logger->warn("BOBYQA problem detected: "
                               + std::string {e.what()} + " at depth "
                               + std::to_string(trialDepth) + "; skipping...");
                    continue;
                } 
                double objectiveFunctionThisDepth
                    = bobyqa.getOptimalObjectiveFunction();
std::cout<<objectiveFunctionThisDepth<<std::endl;
                nObjectiveFunctionEvaluations 
                    = nObjectiveFunctionEvaluations 
                    + bobyqa.getNumberOfObjectiveFunctionEvaluations();
                if (objectiveFunctionThisDepth < optimalDepthObjectiveFunction)
                {
                    optimalOrigin = bobyqa.getOrigin();
                    optimalDepthObjectiveFunction = objectiveFunctionThisDepth;
                    reducedObjectiveFunction = true;
                }
            }
            travelTimeCalculators = bobyqa.releaseTravelTimeCalculatorMap();
            topography = bobyqa.releaseTopography();
            if (!reducedObjectiveFunction)
            {
                logger->warn("Using initial solution for event "
                           + std::to_string(origins[iOrigin].getIdentifier()));
            }
            else
            {
                logger->debug("Using BOBYQA solution for "
                            + std::to_string(origins[iOrigin].getIdentifier()));
            }
*/
//std::cout << "Done" << std::endl;
            // Note something for my own edification
/*
//std::cout << origin.getEpicenter().getLatitude() << " " << origin.getEpicenter().getLongitude() << std::endl;
std::cout << "initial solution: " << bestInitialOrigin.getEpicenter().getLatitude() << " " << bestInitialOrigin.getEpicenter().getLongitude() << " " << bestInitialOrigin.getDepth() << std::endl;
std::cout << "refined solution: " << optimalOrigin.getEpicenter().getLatitude() << " " << optimalOrigin.getEpicenter().getLongitude() << " " << optimalOrigin.getDepth() << std::endl;
std::cout << "real solution: " << origins[iOrigin].getEpicenter().getLatitude() << " " << origins[iOrigin].getEpicenter().getLongitude() << " " << origins[iOrigin].getDepth() << std::endl;
*/
            // Over-write origin 
            origin = optimalOrigin;
            origin.setIdentifier(origins.at(iOrigin).getIdentifier());
            origin.setEventType(origins.at(iOrigin).getEventType());
            // Now compute the uncorrected travel times
            arrivals = origin.getArrivals();
            std::vector<std::pair<ULocator::Station, std::string>>
                stationPhases;
            for (const auto &arrival : arrivals)
            {
                stationPhases.push_back(std::pair {arrival.getStation(),
                                                   arrival.getPhase()});
            }
            // Compute uncorrected times
            std::vector<double> estimateTimes;
            try
            {
                constexpr bool applyCorrection{false};
                auto [xSource, ySource]
                    = programOptions.geographicRegion->geographicToLocalCoordinates(
                         origin.getEpicenter().getLatitude(),
                         origin.getEpicenter().getLongitude());
                travelTimeCalculators->evaluate(stationPhases,
                                                0, //origin.getTime(), 
                                                xSource, ySource,
                                                origin.getDepth(),
                                                &estimateTimes,
                                                applyCorrection);
            }
            catch (const std::exception &e) 
            {
                logger->error(e.what());
            }
            // Write the catalog
            int iArrival = 0;
            for (const auto &arrival : origin.getArrivals())
            {
                double uncorrectedTravelTime = estimateTimes.at(iArrival);
                outCatalog << std::setprecision(16)  
                           << origin.getIdentifier() << "," 
                           << origin.getEpicenter().getLatitude() << ","
                           << origin.getEpicenter().getLongitude() << ","
                           << origin.getDepth() << ","
                           << origin.getTime() << ","
                           << arrival.getStationReference().getNetwork() << ","
                           << arrival.getStationReference().getName() << ","
                           << arrival.getPhase() << ","
                           << arrival.getTime() << ","
                           << arrival.getStandardError() << ","
                           << arrival.getResidual() << ","
                           << uncorrectedTravelTime << ","
                           << arrival.getDistance() << ","
                           << ::eventTypeToString(origin.getEventType()) << ","
                           << nObjectiveFunctionEvaluations << ","
                           << origin.getWeightedRootMeanSquaredError() << std::endl;
                iArrival = iArrival + 1;
            }
        }
        else // Job is forward problem only
        {
            nObjectiveFunctionEvaluations = 1;
            auto epicenter = origin.getEpicenter();
            auto zSource = origin.getDepth();
            auto [xSource, ySource]
                = programOptions.geographicRegion->geographicToLocalCoordinates(
                     epicenter.getLatitude(),
                     epicenter.getLongitude());
            std::vector<std::pair<ULocator::Station, std::string>>
                stationPhases;
            for (const auto &arrival : arrivals)
            {
                stationPhases.push_back(std::pair {arrival.getStation(),
                                                   arrival.getPhase()});
            }
            // Compute residuals for corrected times and uncorrected travel times
            std::vector<double> estimateTimes;
            try
            {
                bool applyCorrection{true};
                travelTimeCalculators->evaluate(stationPhases,
                                                origin.getTime(),
                                                xSource, ySource, zSource,
                                                &estimateTimes,
                                                applyCorrection);
                for (int i = 0; i < static_cast<int>(arrivals.size()); ++i)
                {
                    auto residual = arrivals.at(i).getTime()
                                  - estimateTimes.at(i);
                    arrivals.at(i).setResidual(residual);
                }
                origin.setArrivals(arrivals);
                // Uncorrected estimate times
                applyCorrection = false;
                travelTimeCalculators->evaluate(stationPhases,
                                                0, //origin.getTime(),
                                                xSource, ySource, zSource,
                                                &estimateTimes,
                                                applyCorrection);
            }
            catch (const std::exception &e)
            {
                logger->error(e.what());
            }
            // Write the catalog
            int iArrival = 0;
            for (const auto &arrival : origin.getArrivals())
            {
                // r = obs - est = obs - (tt + ot) = obs - tt - ot
                // tt = obs - r - ot
                double uncorrectedTravelTime = estimateTimes.at(iArrival);
                outCatalog << std::setprecision(16) 
                           << origin.getIdentifier() << "," 
                           << origin.getEpicenter().getLatitude() << ","
                           << origin.getEpicenter().getLongitude() << ","
                           << origin.getDepth() << ","
                           << origin.getTime() << ","
                           << arrival.getStationReference().getNetwork() << ","
                           << arrival.getStationReference().getName() << ","
                           << arrival.getPhase() << ","
                           << arrival.getTime() << ","
                           << arrival.getStandardError() << ","
                           << arrival.getResidual() << ","
                           << uncorrectedTravelTime << ","
                           << arrival.getDistance() << ","
                           << ::eventTypeToString(origin.getEventType()) << ","
                           << nObjectiveFunctionEvaluations << ","
                           << origin.getWeightedRootMeanSquaredError() << std::endl;
                  iArrival = iArrival + 1;
            }
        }
    } 
    outCatalog.close();
    return EXIT_SUCCESS;
}
