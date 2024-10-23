#include <iostream>
#include <filesystem>
#include <string>
#include <fstream>
#include <iomanip>
#include <functional>
#include <vector>
#include <prima/prima.h>
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
#else
#include "logging/standardOut.hpp"
#endif
#include "uLocator/travelTimeCalculatorMap.hpp"
#include "uLocator/rayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "optimizers/objectiveFunctions.hpp"
#include "loadCatalog.hpp"

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
                "Velocity " + std::to_string(velocity)
              + " cannot be less than lower bound " 
              + std::to_string(lowerBound));
        }
        if (velocity > upperBound)
        {
            throw std::invalid_argument(
                "Velocity " + std::to_string(velocity)
              + " cannot be greater than upper bound "
              + std::to_string(upperBound));
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
        estimates.reserve(inOrigins.size()*10);
        stationPhasePairs.reserve(inOrigins.size()); 
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
                    estimates.push_back(
                        arrival.getTime() - arrival.getResidual());
                    catalogResiduals.push_back(arrival.getResidual());
                    spPair.push_back(std::pair {arrival.getStation(), phase});
                }
            }
            if (!arrivalsToKeep.empty())
            {
                origin.setArrivals(arrivalsToKeep);
                origins.push_back(origin);
                stationPhasePairs.push_back(std::move(spPair));
            }
        }
        // Tabulate objective function
        double objectiveFunction = evaluateObjectiveFunction();
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
    void initializeCalculators() const
    {
        logger->debug("Initializing calculators...");
#ifndef NDEBUG
        assert(interfaces.size() == velocities.size());
#endif
        calculators.clear();
        for (const auto &uniqueStation : uniqueStations)
        {
            auto name = uniqueStation.getNetwork()
                      + "." + uniqueStation.getName() + "." + phase;
            try
            {
                logger->debug("Initializing: " + name);
                auto calculator
                    = std::make_unique<ULocator::RayTracer> (uniqueStation,
                                                             phase,
                                                             interfaces,
                                                             velocities,
                                                             logger);
                calculators.insert(uniqueStation, phase, std::move(calculator));
            }
            catch (const std::exception &e)
            {
                 logger->error(e.what());
            }
        }
    }
    double f(const int n, const double x[]) const
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
        int nOrigins = static_cast<int> (origins.size());
        int iObservation{0};
        for (int iOrigin = 0; iOrigin < nOrigins; ++iOrigin)
        {
            auto epicenter = origins.at(iOrigin).getEpicenter();
            auto [xSource, ySource]
                = region->geographicToLocalCoordinates(
                    epicenter.getLatitude(), epicenter.getLongitude());
            auto zSource = origins[iOrigin].getDepth();
            double originTime = origins[iOrigin].getTime();
            const auto &arrivals = origins[iOrigin].getArrivalsReference();
            std::vector<double> travelTimes;
            constexpr bool applyCorrection{false};
            calculators.evaluate(stationPhasePairs.at(iOrigin),
                                 originTime, xSource, ySource, zSource,
                                 &travelTimes,
                                 applyCorrection);
            std::copy(travelTimes.begin(), travelTimes.end(),
                      estimates.begin() + iObservation);
            iObservation = iObservation + travelTimes.size();
        }
#ifndef NDEBUG
        assert(iObservation == estimates.size());
#endif
        auto objectiveFunction = evaluateObjectiveFunction();
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
    double evaluateObjectiveFunction() const
    {
        // Tabulate objective function
        double objectiveFunction{0};
        if (norm == ULocator::Optimizers::IOptimizer::Norm::L1)
        {
            objectiveFunction
                = ::l1(weights, arrivalTimes, estimates, ::Measurement::Standard);
            logger->debug("L1 objective function: "
                        + std::to_string(objectiveFunction));
        }
        else if (norm == ULocator::Optimizers::IOptimizer::Norm::LeastSquares)
        {
            objectiveFunction
                = ::leastSquares(weights, arrivalTimes, estimates,
                                 ::Measurement::Standard);
            logger->debug("Catalog least-squares objective function: "
                        + std::to_string(objectiveFunction));
        }
        else if (norm == ULocator::Optimizers::IOptimizer::Norm::Lp)
        {
            objectiveFunction
                = ::lp(weights, arrivalTimes, estimates, pNorm,
                      ::Measurement::Standard);
            logger->debug("Catalog lp objective function: "
                        + std::to_string(objectiveFunction));
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        return objectiveFunction;
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
    mutable std::vector<double> estimates;
    std::vector<double> weights;
    std::vector<double> interfaces;
    mutable std::vector<double> velocities;
    mutable ULocator::TravelTimeCalculatorMap calculators;
    std::unique_ptr<ULocator::Position::IGeographicRegion> region{nullptr};
    std::string phase{"P"};
    double pNorm{1.5};
    ULocator::Optimizers::IOptimizer::Norm
        norm{ULocator::Optimizers::IOptimizer::Norm::LeastSquares};
    mutable int nEvaluations{0};
};

static void primaCallbackFunction(const double x[],
                                  double *f,
                                  const void *data)
{
     auto objectiveFunction
         = reinterpret_cast<const ObjectiveFunction *> (data);
     auto n = static_cast<int> (objectiveFunction->velocities.size());
     *f = objectiveFunction->f(n, x);
}


struct ProgramOptions
{
    std::filesystem::path resultsDirectory{"/home/bbaker/Codes/uLocator/examples/uuss/utahResults"};
    std::filesystem::path catalogFile{"/home/bbaker/Codes/uLocator/examples/uuss/utahReviewedCatalog.csv"};
    std::filesystem::path trainingFile{"/home/bbaker/Codes/uLocator/examples/uuss/utahReviewedTrainingEvents.csv"};
    std::string phase{"P"};
    std::unique_ptr<ULocator::Position::IGeographicRegion> region{nullptr};
    double pNorm{1.5};
    ULocator::Optimizers::IOptimizer::Norm 
        norm{ULocator::Optimizers::IOptimizer::Norm::LeastSquares};
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
        ("norm", boost::program_options::value<std::string> ()->default_value("l2"),
                "This defines the norm in which to minimize.  This can be l1, l2 (least-squares), lp")
        ("p", boost::program_options::value<double> ()->default_value(options.pNorm),
              "If using norm = lp then this is the exponent.  This must be greater than or equal to 1")
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
    if (vm.count("norm"))
    {
        auto snorm = vm["norm"].as<std::string> ();
        if (snorm == "l2")
        {
            options.norm = ULocator::Optimizers::IOptimizer::Norm::LeastSquares;
        }
        else if (snorm == "l1")
        {
            options.norm = ULocator::Optimizers::IOptimizer::Norm::L1;
        }
        else if (snorm == "lp")
        {
           options.norm = ULocator::Optimizers::IOptimizer::Norm::Lp;
           if (vm.count("p"))
           {
               auto p = vm["p"].as<double> ();
               if (p < 1)
               {
                   throw std::invalid_argument("p must be >= 1");
               }
               if (p == 1)
               {
                   options.pNorm = 1;
                   options.norm = ULocator::Optimizers::IOptimizer::Norm::L1;
               }
               else if (p == 2) 
               {
                   options.pNorm = 2;
                   options.norm
                       = ULocator::Optimizers::IOptimizer::Norm::LeastSquares;
               }
            }
        }
    }
    options.doUtah = true;
    options.region = std::make_unique<ULocator::Position::UtahRegion> ();
    if (vm.count("do_ynp"))
    {
        options.region = std::make_unique<ULocator::Position::YNPRegion> ();
        options.doUtah = false;
    }
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

    auto logger = std::make_shared<UMPS::Logging::StandardOut> ();

    auto phase = programOptions.phase;
    auto initialResidualsFile
        = programOptions.resultsDirectory /
          std::filesystem::path{"initialResiduals." + phase + ".csv"};
    auto finalResidualsFile 
        = programOptions.resultsDirectory /
          std::filesystem::path{"finalResiduals." + phase + ".csv"};
    auto catalogFile = programOptions.catalogFile;
    auto trainingFile = programOptions.trainingFile;
    constexpr int catalogVersion{3};
    auto origins = ::loadCatalog(catalogFile, trainingFile,
                                 *programOptions.region,
                                 logger, catalogVersion);
    // Make some baseline
    std::vector<::Layer> initialModel;
    if (programOptions.doUtah)
    {
        if (phase == "P")
        {
            logger->info("Creating initial Utah P model...");
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
            logger->info("Creating initial Utah S model...");
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
            logger->info("Creating initial YNP P model...");
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
            logger->info("Creating initial YNP S model...");
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
    objectiveFunction.logger = logger;
    objectiveFunction.phase = phase;
    objectiveFunction.region = programOptions.region->clone();
    objectiveFunction.norm = programOptions.norm;
    objectiveFunction.pNorm = programOptions.pNorm;
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

    std::vector<::Layer> layers;
    if (programOptions.doUtah)
    {
        if (phase == "P")
        {
            logger->info("Creating starting Utah P model...");
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
            logger->info("Creating starting Utah S model...");
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
            logger->info("Creating starting YNP P model...");
            layers = std::vector<::Layer>
            {
/*
                ::Layer{-4500,  2512, 1000, 2750 - 1.e-1},
                ::Layer{-1000,  3398, 2750, 3800 - 1.e-1},
                ::Layer{ 2000,  4689, 3800, 5000 - 1.e-1},
                ::Layer{ 5000,  5456, 5000, 5600 - 1.e-1},
                ::Layer{ 8000,  5674, 5600, 5900 - 1.e-1},
                ::Layer{12000,  6250, 5900, 6300 - 1.e-1},
                ::Layer{16000,  6398, 6300, 6450 - 1.e-1},
                ::Layer{21000,  6575, 6450, 7200 - 1.e-1},
                ::Layer{50000,  8200, 7200, 8400}
*/
                ::Layer{-4500,  2720, 1000, 2750 - 1.e-1},
                ::Layer{-1000,  2790, 2750, 4000 - 1.e-1},
                ::Layer{ 2000,  5210, 4000, 5350 - 1.e-1},
                ::Layer{ 5000,  5560, 5350, 5600 - 1.e-1},
                ::Layer{ 8000,  5770, 5600, 5850 - 1.e-1},
                ::Layer{12000,  6070, 5850, 6150 - 1.e-1},
                ::Layer{16000,  6330, 6150, 6450 - 1.e-1},
                ::Layer{21000,  6630, 6450, 7200 - 1.e-1},
                ::Layer{50000,  8000, 7200, 8400}
            };
        }
        else if (phase == "S")
        {
            logger->info("Creating starting YNP S model...");
            layers = std::vector<::Layer>
            {
                ::Layer{-4500,  1729,  500, 1975 - 1.e-1},
                ::Layer{-1000,  2338, 1975, 2700 - 1.e-1},
                ::Layer{ 2000,  3055, 2700, 3200 - 1.e-1},
                ::Layer{ 5000,  3430, 3200, 3500 - 1.e-1},
                ::Layer{ 8000,  3567, 3500, 3650 - 1.e-1},
                ::Layer{12000,  3690, 3650, 3700 - 1.e-1},
                ::Layer{16000,  3750, 3700, 3900 - 1.e-1},
                ::Layer{21000,  3975, 3900, 4400 - 1.e-1},
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

    std::vector<double> xWork(x);
    prima_problem_t primaProblem;
    prima_init_problem(&primaProblem, n);
    primaProblem.x0 = xWork.data();
    primaProblem.calfun = &primaCallbackFunction;
    primaProblem.xl = lowerBounds.data();
    primaProblem.xu = upperBounds.data();

    prima_options_t primaOptions;
    prima_init_options(&primaOptions);
    primaOptions.rhobeg = 25;  // reasonable initial change in velocities
    primaOptions.rhoend = 0.1; // reasonable final change in velocities
    primaOptions.data = &objectiveFunction;
    primaOptions.maxfun = std::max(static_cast<int> (x.size()) + 3, 50);

    prima_result_t primaResult;
    prima_minimize(PRIMA_BOBYQA, primaProblem, primaOptions, &primaResult);
    
    std::copy(primaResult.x, primaResult.x + x.size(), x.begin());
    if (primaResult.f < f0)
    {
        logger->info("Final loss: " + std::to_string(primaResult.f));
    }
    else
    {
        logger->warn("Optimization failed.  Using initial model");
        
        for (size_t i = 0; i < initialModel.size(); ++i)
        {   
            x[i] = initialModel[i].velocity;
        }
    }
    
    prima_free_result(&primaResult); 
    //prima_free_problem(&primaProblem);
 
    for (int i = 0; i < objectiveFunction.interfaces.size(); ++i)
    {
        std::cout << objectiveFunction.interfaces[i] << " "
                  << x[i] << std::endl;
    }
/*
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
    logger->info("Final loss: " + std::to_string(loss));
    for (int i = 0; i < objectiveFunction.interfaces.size(); ++i)
    {
        std::cout << objectiveFunction.interfaces[i] << " "
                  << x[i] << std::endl;
    }
    objectiveFunction.nEvaluations = 0;
*/
    auto fFinal = objectiveFunction.f(x.size(), x.data());
    logger->info("Evaluating final loss: " + std::to_string(fFinal)
               + " ; " 
               + std::to_string( (f0 - fFinal)/f0*100 ) + " pct reduction");
    objectiveFunction.writeToCSV(finalResidualsFile);
    return EXIT_SUCCESS;
}
