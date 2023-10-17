#include <iostream>
#include <filesystem>
#include <iomanip>
#include <numeric>
#include <random>
#include <fstream>
#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
//#include <oneapi/dal.hpp>
//#include <oneapi/dal/algo/knn.hpp>
//#include <oneapi/dal/io/csv.hpp>
//#include <oneapi/dal/table/homogen.hpp>
//#include <oneapi/dal/exceptions.hpp>
#include <oneapi/dal/algo/knn.hpp>
#include <oneapi/dal/table/homogen.hpp>
#include <oneapi/dal/table/row_accessor.hpp>
#include <daal.h>
#include <boost/tokenizer.hpp>
#include "uLocator/position/wgs84.hpp"
#include "uLocator/staticCorrection.hpp"
#include "uLocator/sourceSpecificStationCorrection.hpp"
#include "weightedMedian.hpp"
#include "weightedMean.hpp"

struct Row
{
    Row(const std::vector<std::string> &line)
    {
#ifndef NDEBUG
        assert(line.size() == 16);
#endif
        eventIdentifier = std::stol(line[0]);
        latitude = std::stod(line[1]);
        longitude = std::stod(line[2]);
        depth = std::stod(line[3]);
        originTime = std::stod(line[4]);
        network = line[5];
        station = line[6];
        phase = line[7];
        arrivalTime = std::stod(line[8]);
        standardError = std::stod(line[9]);
        residual = std::stod(line[10]);
        uncorrectedTravelTime = std::stod(line[11]); 
        distance = std::stod(line[12]);
        eventType = line[13];
        weight = 1./standardError;
        weightedRMS = std::stod(line[15]);
    }
    [[nodiscard]] std::string getName() const
    {
        return network + "." + station + "." + phase;
    }
    std::string network;
    std::string station;
    std::string phase;
    std::string eventType;
    double latitude;
    double longitude;
    double distance;
    double depth;
    double originTime;
    double arrivalTime;
    double standardError;
    double residual;
    double uncorrectedTravelTime;
    double weight;
    double weightedRMS;
    int64_t eventIdentifier;
};

std::pair<std::vector<std::pair<std::string, int>>,
          std::vector<::Row> >
loadTable(const std::string &fileName,
          const double maximumResidualForP = 0.7,
          const double maximumResidualForS = 1.4,
          const double maximumWeightedRMS = 0.8)
{
    //event_identifier,latitude,longitude,depth,origin_time,network,station,phase,arrival_time,standard_error,residual,uncorrected_travel_time,event_type
    // Get training events
    std::string line;
    std::vector<int64_t> trainingEventIdentifiers;
    std::ifstream trainingFile(fileName);
    getline(trainingFile, line); // Header
    std::vector<::Row> rows;
    std::vector<std::pair<std::string, int>> workList;
    rows.reserve(15000);
    workList.reserve(1000);
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
        ::Row row(splitLine);
        auto name = row.getName();
        auto residual = row.arrivalTime
                      - (row.uncorrectedTravelTime + row.originTime);
        if (row.weightedRMS > maximumWeightedRMS)
        {
            continue;
        }
        if (row.phase == "P")
        {
            if (std::abs(residual) > maximumResidualForP ||
                std::abs(row.residual) > maximumResidualForP)
            {
                //std::cerr << "Residual execeeded for: " << name << std::endl;
                continue; 
            }
        }
        else if (row.phase == "S")
        {
            if (std::abs(residual) > maximumResidualForS ||
                std::abs(row.residual) > maximumResidualForS)
            {
                //std::cerr << "Residual execeeded for: " << name << std::endl;
                continue; 
            }
        }
#ifndef NDEBUG
        else
        {
            assert(false);
        }
#endif
        rows.push_back(std::move(row));
        bool exists = false;
        for (auto &item : workList)
        {
            if (name == item.first)
            {
                item.second = item.second + 1;
                exists = true;
                break;
            }
        }
        if (!exists)
        {
            workList.push_back( std::pair{name, 1} ); 
        }
    }
    std::sort(workList.begin(), workList.end(),
              [=](const std::pair<std::string, int> &lhs,
                  const std::pair<std::string, int> &rhs)
             {
                  return lhs.second > rhs.second;
             });
    int stationNumber = 1;
    for (const auto &item : workList)
    {
        std::cout << stationNumber << "," << item.first << "," << item.second << std::endl;
        stationNumber++;
    }
    return std::pair {workList, rows}; 
}

void getArrivalTimesAndWeightsForStaticCorrection(
    const std::string &name,
    const std::vector<::Row> &table,
    std::vector<double> *observedArrivalTimes,
    std::vector<double> *predictedArrivalTimes,
    std::vector<double> *weights)
{
    observedArrivalTimes->clear();
    predictedArrivalTimes->clear();
    weights->clear();
    observedArrivalTimes->reserve(table.size());
    predictedArrivalTimes->reserve(table.size());
    weights->reserve(table.size());
    for (const auto &row : table)
    {
        if (row.getName() == name)
        {
            observedArrivalTimes->push_back(row.arrivalTime);
            predictedArrivalTimes->push_back(
               row.originTime + row.uncorrectedTravelTime);
            weights->push_back(row.weight);
        }
    }
}

/// Gets weights and residuals for static corrections 
void getWeightsAndResidualsForStaticCorrection(
    const std::string &name,
    const std::vector<::Row> &table,
    std::vector<double> *residuals,
    std::vector<double> *weights)
{
    residuals->clear();
    weights->clear();
    residuals->reserve(table.size());
    weights->reserve(table.size());
    for (const auto &row : table)
    {
        if (row.getName() == name)
        {
            auto tObs = row.arrivalTime;
            auto tEst = row.originTime + row.uncorrectedTravelTime;
            residuals->push_back(tObs - tEst);
            weights->push_back(row.weight);
        }
    }
}

void getFeaturesAndTargets(
    const std::string &name,
    const std::vector<::Row> &table,
    std::vector<double> *latitudes,
    std::vector<double> *longitudes,
    std::vector<double> *depths,
    std::vector<double> *observedArrivalTimes,
    std::vector<double> *predictedArrivalTimes,
    std::vector<double> *weights,
    const ULocator::StaticCorrection &correction)
{
    latitudes->clear();
    longitudes->clear();
    depths->clear();
    observedArrivalTimes->clear();
    predictedArrivalTimes->clear();
    weights->clear();
 
    latitudes->reserve(table.size());
    longitudes->reserve(table.size());
    depths->reserve(table.size());
    observedArrivalTimes->reserve(table.size());
    predictedArrivalTimes->reserve(table.size());
    weights->reserve(table.size());
    for (const auto &row : table)
    {
        if (row.getName() == name)
        {
            auto tObs = row.arrivalTime;
            auto tEst = row.originTime + row.uncorrectedTravelTime;
            latitudes->push_back(row.latitude);
            longitudes->push_back(row.longitude);
            depths->push_back(row.depth);
            observedArrivalTimes->push_back(tObs);
            predictedArrivalTimes->push_back(correction.evaluate(tEst));
            weights->push_back(row.weight);
        }
    }

}

double computeRMSE(const std::vector<double> &observed,
                   const std::vector<double> &estimated)
{
    auto n = std::min(static_cast<int> (observed.size()),
                      static_cast<int> (estimated.size()));
    double sumResidual2 = 0;
    for (int i = 0; i < n; ++i)
    {
        auto residual = observed[i] - estimated[i];
        sumResidual2 = sumResidual2 + residual*residual;
    }
    return std::sqrt(sumResidual2/std::max(0, n)); 
}

double computeRMSE(const std::vector<double> &observed,
                   const std::vector<double> &estimated,
                   const std::vector<double> &weights)
{
#ifndef NDEBUG
    assert(observed.size() == estimated.size());
    assert(observed.size() == weights.size());
#endif
    double sumOfWeights = std::accumulate(weights.begin(), weights.end(), 0.0);
    if (sumOfWeights <= 0){throw std::invalid_argument("Weight problem");}
    double numerator = 0;
    for (int i = 0; i < static_cast<int> (observed.size()); ++i)
    {
        auto residual = observed[i] - estimated[i]; 
        numerator = numerator + weights[i]*(residual*residual);
    }
    return std::sqrt(numerator/sumOfWeights);
}

/*
void getFeaturesAndTargets(
    const std::string &name,
    const std::vector<::Row> &table,
    std::vector<double> *features,
    std::vector<double> *residuals,
    const double correction = 0,
    const int nFeatures = 3)
{
    if (nFeatures < 2 || nFeatures > 3)
    {
        throw std::invalid_argument("Numer of features must be 2 or 3");
    }
    constexpr int zone{12};
    ULocator::Position::WGS84 position;
    features->clear();
    residuals->clear();
    features->reserve(3*table.size());
    residuals->reserve(table.size());
    for (const auto &row : table)
    {
        if (row.getName() == name)
        {
            auto tObs = row.arrivalTime;
            auto tEst = row.originTime + row.uncorrectedTravelTime;
            ULocator::Position::WGS84 position{
                row.latitude, row.longitude, zone};
            features->push_back(position.getEasting()*1.e-3);
            features->push_back(position.getNorthing()*1.e-3);
            if (nFeatures == 3){features->push_back(row.depth*1.e-3);}
            residuals->push_back(tObs - (tEst + correction));
        }
    }
}
*/

///
std::vector<std::vector<int>> 
    createFolds(const int nObservations, const int nFolds = 10,
                const int seed =-1)
{
    if (nFolds < 1)
    {
        throw std::invalid_argument("Number of folds must be positive");
    }
    if (nObservations < nFolds)
    {
        throw std::invalid_argument("Fewer observations than folds");
    }
    std::mt19937 randomGenerator{std::random_device{}()};
    if (seed >= 0){randomGenerator.seed(seed);}
    std::vector<int> rows(nObservations);
    std::iota(rows.begin(), rows.end(), 0);
    std::vector<std::vector<int>> result(nFolds);
    if (nFolds == 1)
    {
        result[0] = std::move(rows);
        return result;
    }
    auto dFold = static_cast<double> (nObservations)/static_cast<double> (nFolds);
    int i0 = 0;
    int i1 = dFold*1; 
    for (int iFold = 0; iFold < nFolds; ++iFold)
    {
        i0 = static_cast<int> ( dFold*iFold );
        if (iFold == 0){i0 = 0;}
        i1 = static_cast<int> ( dFold*(iFold + 1) );
        if (iFold == nFolds - 1){i1 = nObservations;}
        int nExamples = i1 - i0;
        // Sample rows
        std::vector<int> sampledRows;
        std::sample(rows.begin(), rows.end(), std::back_inserter(sampledRows),
                    nExamples, randomGenerator); 
        std::sort(sampledRows.begin(), sampledRows.end());
        // Remove
        for (const auto &row : sampledRows)
        {
            rows.erase(std::remove(rows.begin(), rows.end(), row),
                       rows.end());
        }
        // Save sampled rows
        result[iFold] = sampledRows;
    }
#ifndef NDEBUG
    assert(rows.empty());
    int nRows = 0;
    for (const auto &fold : result)
    {
        nRows = nRows + fold.size();
        //for (const auto &r : fold){std::cout << r << " ";}
        //std::cout << std::endl;
    }
    assert(nRows == nObservations);
#endif
    return result;
}

/// Trains the KNN Model
void trainKNN(const std::string &name,
              const std::string &resultsFile,
              const std::vector<double> &latitudes,
              const std::vector<double> &longitudes,
              const std::vector<double> &depths,
              const std::vector<double> &observedArrivalTimes,
              const std::vector<double> &predictedArrivalTimes,
              const std::vector<double> &weights,
              const int nFolds = 10,
              const std::vector<int> nNeighborsSearch = std::vector<int>{1, 2, 3, 5, 7, 10, 13, 15, 17, 20, 25, 30, 40, 50, 100},
              const std::vector<double> maximumDistancesSearch = std::vector<double>{5000,
                                                                                     10000,
                                                                                     15000,
                                                                                     20000,
                                                                                     25000,
                                                                                     30000,
                                                                                     35000,
                                                                                     40000,
                                                                                     50000},
              auto method = ULocator::SourceSpecificStationCorrection::EvaluationMethod::InverseDistanceWeighted,
              bool verbose = true,
              const int seed = 8434)
{
    std::filesystem::path pathToResult{resultsFile};
    std::filesystem::path gridSearchPath{"knnGridSearch"};
    if (!pathToResult.parent_path().empty())
    {
        gridSearchPath = pathToResult.parent_path()/gridSearchPath; 
    }
    if (!std::filesystem::exists(gridSearchPath))
    {
        std::filesystem::create_directories(gridSearchPath);
    } 
    std::filesystem::path outputGridSearchFileName
       = gridSearchPath/std::filesystem::path{name + ".csv"};
    std::ofstream outputGridSearchFile{outputGridSearchFileName};
    std::vector<std::string> splitName;
    boost::split(splitName, name, boost::is_any_of("."));
    auto nObservations = static_cast<int> (observedArrivalTimes.size());
    auto dataFolds = ::createFolds(nObservations, nFolds, seed);
    double minimumRMSE = std::numeric_limits<double>::max();
    double baselineRMSE = 0;
    double optimumDistance = 35000;
    int optimumNeighbors = 10;
    bool improvementExists{false};
    // Grid-search
    for (const auto &nNeighbors : nNeighborsSearch)
    {
        for (const auto &maximumDistance : maximumDistancesSearch)
        {
            std::vector<double> allZeros(nObservations, 0.0);
            std::vector<double> observedResidualsAllFolds;
            std::vector<double> maxDistanceResidualsAllFolds;
            std::vector<double> idwResidualsAllFolds;
            std::vector<double> weightedResidualsAllFolds;
            std::vector<double> weightsAllFolds;
            // Loop on folds
            for (const auto &fold : dataFolds)
            {
                std::vector<int> trainingRows(nObservations);
                std::iota(trainingRows.begin(), trainingRows.end(), 0);
                auto nValidate = static_cast<int> (fold.size());
                auto nTrain = nObservations - nValidate;
                // Extract validation features and remove validation rows from
                // training 
                std::vector<double> latitudesValidate(nValidate);
                std::vector<double> longitudesValidate(nValidate);
                std::vector<double> depthsValidate(nValidate);
                std::vector<double> residualsValidate(nValidate);
                std::vector<double> weightsValidate(nValidate);
                std::vector<double> zeros(nValidate, 0.0);
                for (int i = 0; i < nValidate; ++i)
                {
                    auto row = fold[i];
                    latitudesValidate[i]  = latitudes[row];
                    longitudesValidate[i] = longitudes[row];
                    depthsValidate[i]     = depths[row];
                    residualsValidate[i]  = observedArrivalTimes[row]
                                          - predictedArrivalTimes[row];
                    weightsValidate[i] = weights[row];
                    trainingRows.erase(std::remove(trainingRows.begin(),
                                                   trainingRows.end(), row),
                                       trainingRows.end());
                }
#ifndef NDEBUG
                assert(static_cast<int> (trainingRows.size()) == nTrain);
#endif
                // Extract the testing features  
                std::vector<double> latitudesTrain(nTrain);
                std::vector<double> longitudesTrain(nTrain);
                std::vector<double> depthsTrain(nTrain);
                std::vector<double> observedArrivalTimesTrain(nTrain);
                std::vector<double> predictedArrivalTimesTrain(nTrain);
                std::vector<double> weightsTrain(nTrain);
                for (int i = 0; i < nTrain; ++i)
                {
                    auto row = trainingRows[i];
                    latitudesTrain[i]  = latitudes[row];
                    longitudesTrain[i] = longitudes[row];
                    depthsTrain[i]     = depths[row];
                    observedArrivalTimesTrain[i] = observedArrivalTimes[row];
                    predictedArrivalTimesTrain[i] = predictedArrivalTimes[row];
                    weightsTrain[i] = weights[row];
                }
                // Train the model
                ULocator::SourceSpecificStationCorrection model;
                model.initialize(splitName[0], splitName[1], splitName[2]);
                model.train(latitudesTrain, longitudesTrain, depthsTrain,
                            observedArrivalTimesTrain, 
                            predictedArrivalTimesTrain,
                            weightsTrain,
                            nNeighbors,
                            maximumDistance);
                // Apply the model
                auto weightedResiduals
                    = model.evaluate(latitudesValidate,
                                     longitudesValidate,
                                     depthsValidate,
                                     ULocator::SourceSpecificStationCorrection::
                                     EvaluationMethod::WeightedAverage);
                auto maxDistanceResiduals
                    = model.evaluate(latitudesValidate,
                                     longitudesValidate,
                                     depthsValidate,
                                     ULocator::SourceSpecificStationCorrection::
                                     EvaluationMethod::MaximumDistanceWeighted);
                auto idwResiduals
                    = model.evaluate(latitudesValidate,
                                     longitudesValidate,
                                     depthsValidate,
                                     ULocator::SourceSpecificStationCorrection::
                                     EvaluationMethod::InverseDistanceWeighted);
                // Tabulate the RMSE
                baselineRMSE
                    = ::computeRMSE(residualsValidate, zeros,
                                    weightsValidate);
                auto weightedRMSE
                    = ::computeRMSE(residualsValidate, weightedResiduals,
                                    weightsValidate);
                auto maxDistanceRMSE
                    = ::computeRMSE(residualsValidate, maxDistanceResiduals,
                                    weightsValidate);
                auto idwRMSE
                    = ::computeRMSE(residualsValidate, idwResiduals,
                                    weightsValidate);
                observedResidualsAllFolds.insert(
                    observedResidualsAllFolds.end(),
                    residualsValidate.begin(),
                    residualsValidate.end());
                weightedResidualsAllFolds.insert(
                    weightedResidualsAllFolds.end(),
                    weightedResiduals.begin(),
                    weightedResiduals.end());
                maxDistanceResidualsAllFolds.insert(
                    maxDistanceResidualsAllFolds.end(),
                    maxDistanceResiduals.begin(),
                    maxDistanceResiduals.end());
                idwResidualsAllFolds.insert(
                    idwResidualsAllFolds.end(),
                    idwResiduals.begin(),
                    idwResiduals.end());
                weightsAllFolds.insert(
                    weightsAllFolds.end(),
                    weightsValidate.begin(),
                    weightsValidate.end()); 
           //std::cout << baselineRMSE  << " " << weightedRMSE << " " << maxDistanceRMSE << " " << idwRMSE << std::endl;
            } // Loop on folds
            baselineRMSE
                = ::computeRMSE(observedResidualsAllFolds,
                                allZeros, weightsAllFolds);
            auto weightedRMSE
                = ::computeRMSE(observedResidualsAllFolds,
                                weightedResidualsAllFolds, weightsAllFolds);
            auto maxDistanceRMSE
                = ::computeRMSE(observedResidualsAllFolds,
                                maxDistanceResidualsAllFolds, weightsAllFolds);
            auto idwRMSE
                = ::computeRMSE(observedResidualsAllFolds,
                                idwResidualsAllFolds, weightsAllFolds);
//std::cout << baselineRMSE << " " << weightedRMSE << " " << maxDistanceRMSE << " " << idwRMSE << std::endl;
            if (baselineRMSE > 2)
            {
                std::cerr << "Large baseline RMSE for" << name << "!"
                          << std::endl;
/*
    for (int i = 0; i < observedResidualsAllFolds.size();++i)
 {
 std::cout << observedResidualsAllFolds[i] << " " << weightedResidualsAllFolds[i] << std::endl;
}
*/
            }
            if (verbose)
            {
                std::cout << nFolds
                          << " Cross-validated Baseline, Average RMS, Clipped RMS, IDW RMS "
                          << baselineRMSE << " " << weightedRMSE << " "
                          << maxDistanceRMSE << " " << idwRMSE << std::endl;
            }
            double preferredRMSE = weightedRMSE;
            if (method == ULocator::SourceSpecificStationCorrection::EvaluationMethod::WeightedAverage)
            {
                preferredRMSE = weightedRMSE;
            }
            else if (method == ULocator::SourceSpecificStationCorrection::EvaluationMethod::MaximumDistanceWeighted)
            {
                preferredRMSE = maxDistanceRMSE;
            }
            else if (method == ULocator::SourceSpecificStationCorrection::EvaluationMethod::InverseDistanceWeighted)
            {
                preferredRMSE = idwRMSE;
            }
#ifndef NDEBUG
            else
            {
               assert(false);
            }
#endif
            // Update optimum in grid search
            if (preferredRMSE < minimumRMSE) 
            {
                minimumRMSE = preferredRMSE;
                optimumDistance = maximumDistance;
                optimumNeighbors = nNeighbors;
                if (preferredRMSE < baselineRMSE)
                {
                    improvementExists = true;
                } 
            }
            outputGridSearchFile << nNeighbors << "," 
                                 << maximumDistance << "," 
                                 << baselineRMSE << ","
                                 << weightedRMSE << ","
                                 << maxDistanceRMSE << ","
                                 << idwRMSE << std::endl;
        } // Loop on maximum distance
        outputGridSearchFile << std::endl;
    } // Loop on neighbors 
    outputGridSearchFile.close();
    // Fit to everything
    if (improvementExists)
    {
        ULocator::SourceSpecificStationCorrection model;
        model.initialize(splitName[0], splitName[1], splitName[2]);
        model.train(latitudes, longitudes, depths,
                    observedArrivalTimes,
                    predictedArrivalTimes,
                    weights,
                    optimumNeighbors);
        try
        {
            model.save(resultsFile); 
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            return;
        }
        double rmseReduction = (baselineRMSE - minimumRMSE)/baselineRMSE*100;
        std::cout << "Optimum {distance,neighbors} for " << name
                  << " is "
                  << "{" << optimumDistance << "," << optimumNeighbors << "}"
                  << ".  The baseline and optimum RMSE is " 
                  << baselineRMSE << " " << minimumRMSE
                  << ".  RMSE reduction " << rmseReduction
                  << " (pct)" << std::endl;
    }
    else
    {
        std::cerr << "No improvement seen for " << name << std::endl;
    }
}
                
/*
/// Trains the KNN model
void trainKNN(std::vector<double> &features,
              std::vector<double> &targets,
              const int nNeighbors,
              const int nFeatures = 3,
              const int nFolds = 10,
              const double maxDistance = 35,
              const int seed = 8434)
{
    // Split into nFolds
    auto nObservations = static_cast<int> (targets.size());
    auto dataFolds = ::createFolds(nObservations, nFolds, seed);
    namespace dknn = oneapi::dal::knn;
    daal::services::Status status;
    double averageBaselineRMS{0};
    double averageFoldRMS{0};
    double averageIDWFoldRMS{0};
    double averageClippedFoldRMS{0};
    for (const auto &fold : dataFolds)
    {
        auto nValidationObservations = static_cast<int> (fold.size());
        // Define the training rows and validation data 
        std::vector<int> trainingRows(nObservations);
        std::iota(trainingRows.begin(), trainingRows.end(), 0);
        auto nTrainingObservations = nObservations - nValidationObservations;
        std::vector<double> validationFeatures(nValidationObservations*nFeatures);
        std::vector<double> validationTargets(nValidationObservations);
        // Extract validation features and remove validation rows from training 
        for (int i = 0; i < nValidationObservations; ++i)
        {
            auto row = fold[i];
            validationFeatures[nFeatures*i + 0] = features[nFeatures*row + 0];
            validationFeatures[nFeatures*i + 1] = features[nFeatures*row + 1];
            if (nFeatures == 3)
            {
                validationFeatures[nFeatures*i + 2] = features[nFeatures*row + 2];
            }
            validationTargets[i] = targets[row];

            trainingRows.erase(std::remove(trainingRows.begin(),
                                           trainingRows.end(), row),
                               trainingRows.end());
        }
#ifndef NDEBUG
        assert(trainingRows.size() == nTrainingObservations);
#endif
        // Populate the training data
        std::vector<double> trainingFeatures(nTrainingObservations*nFeatures);
        std::vector<double> trainingTargets(nTrainingObservations);
        // Extract validation features and remove validation rows from training 
        for (int i = 0; i < nTrainingObservations; ++i)
        {
            auto row = trainingRows[i];
            trainingFeatures[nFeatures*i + 0] = features[nFeatures*row + 0]; 
            trainingFeatures[nFeatures*i + 1] = features[nFeatures*row + 1]; 
            if (nFeatures == 3)
            {
                trainingFeatures[nFeatures*i + 2] = features[nFeatures*row + 2]; 
            }
            trainingTargets[i] = targets[row];
        }
        // Populate tables for DAL
        auto XTrain
            = oneapi::dal::homogen_table::wrap<double>
              (trainingFeatures.data(), nTrainingObservations, nFeatures);
        //auto yTrain
        //    = oneapi::dal::homogen_table::wrap<double>
        //      (trainingTargets.data(),  nTrainingObservations, 1);
        const auto knnDescriptor
            = dknn::descriptor<double,
                               dknn::method::kd_tree,
                               oneapi::dal::knn::task::search,
                               oneapi::dal::minkowski_distance::descriptor<double>>
              (nNeighbors);
        // Train the model
        auto trainingResult = oneapi::dal::train(knnDescriptor, XTrain);//, yTrain);
        // Apply the model 
        auto XValidate
            = oneapi::dal::homogen_table::wrap<double>
              (validationFeatures.data(), nValidationObservations, nFeatures);
        //auto yValidate
        //    = oneapi::dal::homogen_table::wrap<double>
        //     (validationTargets.data(),  nValidationObservations, 1); 
        const auto testResult
            = oneapi::dal::infer(knnDescriptor, XValidate, trainingResult.get_model());
        // Tabulate the RMS 
        const auto &trainingIndices = testResult.get_indices();
#ifndef NDEBUG
        assert(trainingIndices.get_column_count() == nNeighbors); 
        assert(trainingIndices.get_row_count() == nValidationObservations);
#endif
        std::vector<int> nearestNeighbors(nNeighbors);
        oneapi::dal::row_accessor<const double> accessor{trainingIndices};
        double baselineRMS = 0;
        double foldRMS = 0;
        double idwFoldRMS = 0;
        double clippedFoldRMS = 0;
        for (int64_t row = 0; row < nValidationObservations; ++row)
        {
            double xi = validationFeatures[nFeatures*row + 0]; 
            double yi = validationFeatures[nFeatures*row + 1]; 
            double zi = 0;
            //if (nFeatures == 3){zi = validationFeatures[nFeatures*row + 2];}
            double observedResidual = validationTargets[row]; 
            const auto rowValues = accessor.pull({row, row + 1});
#ifndef NDEBUG
            assert(rowValues.get_count() == nNeighbors);
#endif
            double idwResidual = 0;
            double idwSumOfWeights = 0;
            double averageResidual = 0;
            double clippedSumOfWeights = 0;
            double clippedResidual = 0;
            for (int i = 0; i < nNeighbors; ++i)
            {
                auto trainingRow = static_cast<int> (std::round(rowValues[i]));
                double xk = trainingFeatures[nFeatures*trainingRow + 0];
                double yk = trainingFeatures[nFeatures*trainingRow + 1];
                double zk = 0;
                if (nFeatures == 3)
                {
                    zk = trainingFeatures[nFeatures*trainingRow + 2];
                }
                double residual = trainingTargets[trainingRow];
                // Inverse distance weighting
                auto dx = xi - xk;
                auto dy = yi - yk;
                auto dz = zi - zk;
                // Can't be closer than a meter
                double distance
                    = std::max(1.e-3, std::sqrt(dx*dx + dy*dy + dz*dz));
                double idwWeight = 1/distance;
                // Can't be further than x km
                idwResidual = idwWeight*residual;
                idwSumOfWeights = idwSumOfWeights + idwWeight; 
                // Average weights
                averageResidual = averageResidual + residual;
                // Clipped
                double clipWeight = 1;
                if (distance > maxDistance){clipWeight = 0;}
                clippedSumOfWeights = clippedSumOfWeights + clipWeight;
                clippedResidual = clippedResidual + clipWeight*residual;
                //std::cout << distance << " " << residual << std::endl;
            }
            averageResidual = averageResidual/nNeighbors;
            idwResidual = idwResidual/idwSumOfWeights;
            clippedResidual = clippedResidual/std::max(1.e-10, clippedSumOfWeights);
            baselineRMS = baselineRMS
                        + std::pow(observedResidual, 2);
            foldRMS = foldRMS
                    + std::pow(observedResidual - averageResidual, 2);
            idwFoldRMS = idwFoldRMS
                       + std::pow(observedResidual -  idwResidual, 2);
            clippedFoldRMS = clippedFoldRMS
                           + std::pow(observedResidual - clippedResidual, 2);
            //std::cout << observedResidual << "," << averageResidual << "," << idwResidual << std::endl;
        } // Loop on validation set
        baselineRMS = std::sqrt(baselineRMS/nValidationObservations);
        foldRMS = std::sqrt(foldRMS/nValidationObservations);
        idwFoldRMS = std::sqrt(idwFoldRMS/nValidationObservations); 
        clippedFoldRMS = std::sqrt(clippedFoldRMS/nValidationObservations);
std::cout << foldRMS << std::endl;
        averageBaselineRMS = averageBaselineRMS + baselineRMS/nFolds;
        averageFoldRMS = averageFoldRMS + foldRMS/nFolds;
        averageIDWFoldRMS = averageIDWFoldRMS + idwFoldRMS/nFolds;
        averageClippedFoldRMS = averageClippedFoldRMS + clippedFoldRMS/nFolds;
        //std::cout << "Fold RMS, IDW Fold RMS: {"
        //          << foldRMS << "," << idwFoldRMS << "}" << std::endl;
break;
    }
    std::cout << std::setprecision(8) << nFolds << " Cross-validated Baseline, Average RMS, Clipped RMS, IDW RMS: "
              << averageBaselineRMS << "," << averageFoldRMS << ","
              << averageClippedFoldRMS << "," << averageIDWFoldRMS
              << " for " << nObservations << " observations " << std::endl;
}
*/

struct ProgramOptions
{
    std::string correctionsFile{"correctionsArchive.h5"};
    std::string eventsFile; 
    double maximumResidualForP{0.7};
    double maximumResidualForS{1.4};
    double maximumWeightedRMS{1};
    int minimumObservationsForStaticCorrection{20};
    int minimumObservationsForSSSC{250}; // For KNN
    int nFolds{5}; // Number of folds
    bool doHelp{false};
    bool doMedian{true};
};

[[nodiscard]] ::ProgramOptions parseCommandLineOptions(int argc, char *argv[])
{
    ::ProgramOptions options;
    boost::program_options::options_description desc(
R"""(
This utility will compute the static and source specific station corrections.  Example usage:
     computeCorrections --catalog_file=relocatedUtahCatalog.csv --corrections_file=correctionsArchive.h5
Allowed options)""");
    desc.add_options()
        ("help",         "Produces this help message")
        ("events_file",  boost::program_options::value<std::string> ()->default_value("relocatedUtahCatalog.csv"),
                         "A CSV with the events located by the location utility")
        ("corrections_file",  boost::program_options::value<std::string> ()->default_value(options.correctionsFile),
                         "The output HDF5 archive with the static corrections and KNN information")
        ("maximum_weighted_rms", boost::program_options::value<double> ()->default_value(options.maximumWeightedRMS),
                         "If the event's weighted RMS exceeds this time in seconds then all phase arrivals will be skipped")
        ("maximum_residual_for_p", boost::program_options::value<double> ()->default_value(options.maximumResidualForP),
                         "If the absolute value of any P residual exceeds this time in seconds then it will not be used for either correction")
        ("maximum_residual_for_s", boost::program_options::value<double> ()->default_value(options.maximumResidualForS),
                         "If the absolute value of any S residual exceeds this time in seconds then it will not be used for either correction")
        ("min_observations_for_static", boost::program_options::value<int> ()->default_value(options.minimumObservationsForStaticCorrection),
                         "The minimum number of observations required to compute a static correction")
        ("min_observations_for_ssst", boost::program_options::value<int> ()->default_value(options.minimumObservationsForSSSC),
                         "The minimum number of observations required to compute a source-specific station correction")
        ("n_folds", boost::program_options::value<int> ()->default_value(options.nFolds),
                    "The number of folds in the KNN hyper-parameter optimization.")
        ("do_mean_static", "If present then the static corrections will be computed with the mean residual.  Otherwise, the static correcions will be computed with the median residual.");
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
    if (vm.count("events_file"))
    {
        auto eventsFile = vm["events_file"].as<std::string> ();
        if (!std::filesystem::exists(eventsFile))
        {
            throw std::invalid_argument("Events file does not exist");
        }
        options.eventsFile = eventsFile;
    }
    else
    {
        throw std::runtime_error("Events file not set");
    }
    if (vm.count("corrections_file"))
    {
        auto correctionsFile = vm["corrections_file"].as<std::string> ();  
        options.correctionsFile = correctionsFile;
    }
    else
    {
        throw std::runtime_error("Corrections file not set");
    }
    if (vm.count("maximum_weighted_rms"))
    {
        auto maxWeightedRMS = vm["maximum_weighted_rms"].as<double> ();
        if (maxWeightedRMS < 0)
        {   
            throw std::invalid_argument("Max weightedRMS must be non-negative");
        }
        options.maximumWeightedRMS = maxWeightedRMS;
    }
    if (vm.count("maximum_residual_for_p"))
    {
        auto maxResidual = vm["maximum_residual_for_p"].as<double> (); 
        if (maxResidual < 0)
        {
            throw std::invalid_argument("Max p residual must be non-negative");
        }
        options.maximumResidualForP = maxResidual; 
    }
    if (vm.count("maximum_residual_for_s"))
    {
        auto maxResidual = vm["maximum_residual_for_s"].as<double> ();  
        if (maxResidual < 0) 
        {
            throw std::invalid_argument("Max s residual must be non-negative");
        }
        options.maximumResidualForS = maxResidual; 
    }
    if (vm.count("n_folds"))
    {    
        auto nFolds = vm["n_folds"].as<int> ();  
        if (nFolds < 1)
        {
            throw std::invalid_argument("Number of folds must be positive");
        }
        options.nFolds = nFolds;
    }
    if (vm.count("min_observations_for_static"))
    {    
        auto minObservations = vm["min_observations_for_static"].as<int> ();  
        if (minObservations < 1) 
        {
            throw std::invalid_argument("Minimum number of observations for static correction must be positive");
        }
        options.minimumObservationsForStaticCorrection = minObservations;
    }
    if (vm.count("min_observations_for_sssc"))
    {
        auto minObservations = vm["min_observations_for_sssc"].as<int> ();  
        if (minObservations < 1)
        {
            throw std::invalid_argument("Minimum number of observations for SSSC must be positive");
        }
        options.minimumObservationsForSSSC = minObservations;
    }
    options.doMedian = true;
    if (vm.count("do_mean_static"))
    {   
        options.doMedian = false;
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

    std::string resultsFile = programOptions.correctionsFile; //{"correctionsArchive.h5"};
    std::string locatorOutputFile = programOptions.eventsFile; //{"relocatedUtahCatalog.csv"};
    // Parameters for weighted median/mean
    bool doMedian = programOptions.doMedian;
    double maximumResidualForP = programOptions.maximumResidualForP;// 0.7;
    double maximumResidualForS = programOptions.maximumResidualForS;//1.4;
    double maximumWeightedRMS = programOptions.maximumWeightedRMS;
    int minimumObservationsForStaticCorrection = programOptions.minimumObservationsForStaticCorrection; //{20};
    // Parameters for KNN and SSST
    const int64_t nFolds = programOptions.nFolds; //{5}; // n-fold cross-validation for KNN fitting
    const int64_t minObservationsForSSSC = programOptions.minimumObservationsForSSSC;// {250}; // 250 observations for KNN fitting
    int seed = 86332;
    constexpr int64_t nFeatures{3}; // UTM X, UTM Y, Depth
    const std::vector<int> nNeighborsSearch = std::vector<int>{//100,  
                                                               75, 50, 40, 30, 25, 20, 17, 15, 13, 10, 7, 5, 3, 2, 1};
    const std::vector<double> maximumDistancesSearch = std::vector<double>{//2500,
                                                                           //5000,
                                                                           //10000,
                                                                           15000,
                                                                           20000,
                                                                           25000,
                                                                           30000,
                                                                           50000};
    // Load the data
    auto [workList, totalDataTable] = ::loadTable(locatorOutputFile, //"newTrainingUtahCatalog.0.csv",
                                                  maximumResidualForP,
                                                  maximumResidualForS,
                                                  maximumWeightedRMS); 
    std::cout << "Tabulating median static corrections..." << std::endl;
    //std::ofstream correctionsFile("medianUtahCorrections.csv");
    //correctionsFile << "network,station,phase,correction" << std::endl;
    //std::vector<std::pair<std::string, double>> corrections;
    std::vector<std::pair<std::string, ULocator::StaticCorrection>> corrections;
    for (const auto &workItem : workList)
    {
        std::vector<double> residuals, observedArrivalTimes, predictedArrivalTimes, weights;
        std::vector<std::pair<double, int>> workSpace;
        if (workItem.second >= minimumObservationsForStaticCorrection)
        {
            std::vector<std::string> splitName;
            boost::split(splitName, workItem.first, boost::is_any_of("."));
            ULocator::StaticCorrection staticCorrection;
            staticCorrection.initialize(splitName[0],
                                        splitName[1],
                                        splitName[2],
                                        0);
            getArrivalTimesAndWeightsForStaticCorrection(workItem.first,
                                                         totalDataTable,
                                                         &observedArrivalTimes,
                                                         &predictedArrivalTimes,
                                                         &weights);
            if (doMedian)
            {
                staticCorrection.train(observedArrivalTimes,
                                       predictedArrivalTimes,
                                       weights,
                                       ULocator::StaticCorrection::Method::Median);
            }
            else
            {
                staticCorrection.train(observedArrivalTimes,
                                       predictedArrivalTimes,
                                       weights,
                                       ULocator::StaticCorrection::Method::Mean);
            }
            staticCorrection.save(resultsFile);
            ::getWeightsAndResidualsForStaticCorrection(
                workItem.first, totalDataTable,
                &residuals, &weights);
            auto medianCorrection
                = ::weightedMedian(residuals, weights, workSpace);
             auto meanCorrection
                = ::weightedMean(residuals, weights);
            // Currently we have computed
            //   obs - est = bias 
            // and know the systematic bias
            // Say, this bias is 1, i.e.,
            //   obs - est = 2 - 1 = 1.  
            // Now, we wish to remove this bias, i.e.,
            //   obs - est - bias = obs - (est + bias) = 0
            // Hence, our correction is simply the bias.
            std::cout << "Median correction for " << workItem.first << " is "
                      << medianCorrection << " (s)" << " " << staticCorrection.evaluate(0) << std::endl;
            //corrections.push_back(std::pair {workItem.first, medianCorrection});
            corrections.push_back(std::pair {workItem.first, staticCorrection});
        }
    }
    //correctionsFile.close();
    // Fit the KNN models
    for (const auto &workItem : workList)
    {
        if (workItem.second >= minObservationsForSSSC)
        {
            // Apply first model
            bool found{false};
            ULocator::StaticCorrection correction;
            for (const auto &c : corrections)
            {
                if (workItem.first == c.first)
                {
                    correction = c.second;
                    found = true;
                    break;
                }
            }
#ifndef NDEBUG
            assert(found);
#endif
            std::vector<double> latitudes, longitudes, depths;
            std::vector<double> observedArrivalTimes, predictedArrivalTimes, residuals, weights;
            ::getFeaturesAndTargets(workItem.first, totalDataTable,
                                    &latitudes,
                                    &longitudes,
                                    &depths,
                                    &observedArrivalTimes,
                                    &predictedArrivalTimes,
                                    &weights,
                                    correction);

//std::vector<int> nNeighborsSearch{15};
//std::vector<double> maxDistances{35000};
//int nFeatures = 3;
            ::trainKNN(workItem.first,
                       resultsFile,
                       latitudes, longitudes, depths,
                       observedArrivalTimes, predictedArrivalTimes, weights,
                       nFolds,
                       nNeighborsSearch,
                       maximumDistancesSearch,
                       ULocator::SourceSpecificStationCorrection::EvaluationMethod::InverseDistanceWeighted,
                       false,
                       seed);
/*
            std::vector<double> trainingFeatures;
            std::vector<double> trainingTargets;
            ::getFeaturesAndTargets(workItem.first, totalDataTable,
                                    &trainingFeatures,
                                    &trainingTargets,
                                    medianCorrection,
                                    nFeatures);
            // Fit a model
            ::trainKNN(trainingFeatures, trainingTargets, nNeighbors, nFeatures, nFolds, maxDistance);
break;
*/
        }
    }
    
    //const auto dataSource = oneapi::dal::csv::data_source("newUtahCatalog.csv");
    //const auto allDataTable = dal::read<dal::table> (dataSource);
    //dal::column_accessor<const int64_t> eventIdentifierAccessor{allDataTable};
    //dal::array<int64_t> evids = eventIdentifierAccessor.pull(0);


//std::cout << allData.get_row_count() << std::endl;
    //auto readOptions = dataSource.get_read_options();

    return EXIT_SUCCESS;
}
