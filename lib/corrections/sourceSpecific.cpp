#include <iomanip>
#include <iostream>
#include <filesystem>
#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <oneapi/dal/algo/knn.hpp>
#include <oneapi/dal/table/homogen.hpp>
#include <oneapi/dal/table/row_accessor.hpp>
#include <hdf5.h>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uLocator/corrections/sourceSpecific.hpp"
#include "uLocator/position/geographicRegion.hpp"
#include "uLocator/position/wgs84.hpp"

#define KNN_GROUP_NAME "KNN"
#define FEATURES_NAME "KNNTrainingFeaturesXYZ"
#define TARGETS_NAME  "KNNTargetsResidualsWeights"
#define NEIGHBORS_NAME "NumberOfNeighbors"
#define MAXIMUM_DISTANCE_NAME "MaximumDistance"
#define EVALUATION_METHOD_NAME "EvaluationMethod"
#define OPTIMIZATION_ALGORITHM_NAME "OptimizationAlgorithm"

using namespace ULocator::Corrections;

class SourceSpecific::SourceSpecificImpl
{
public:
    void fillFeatureVector(
        const std::vector<double> &xSources,
        const std::vector<double> &ySources,
        const std::vector<double> &zSources,
        const bool isTraining)
    {
        bool doDepth = true;
        if (mFeatures != 3){doDepth = false;}
#ifndef NDEBUG
        assert(xSources.size() == ySources.size());
        if (doDepth){assert(xSources.size() == zSources.size());}
#endif
        auto nObservations = static_cast<int> (xSources.size());
        std::vector<double> featuresMatrix(mFeatures*nObservations, 0);
        for (int i = 0; i < nObservations; ++i)
        {
            auto i0 = mFeatures*i;
            featuresMatrix[i0] = xSources[i];
            featuresMatrix[i0 + 1] = ySources[i];
            if (doDepth){featuresMatrix[i0 + 2] = zSources[i];}
        }
        if (isTraining)
        {
            mTrainingFeaturesMatrix = std::move(featuresMatrix);
        }
        else
        {
            mFeaturesMatrix = std::move(featuresMatrix);
        }
    }
    /// Fills the target vector
    void fillTargetVector(
        const std::vector<double> &observedArrivalTimes,
        const std::vector<double> &predictedArrivalTimes,
        const std::vector<double> &weights)
    {
#ifndef NDEBUG
        assert(observedArrivalTimes.size()  == weights.size());
	assert(predictedArrivalTimes.size() == weights.size());
#endif
        auto nObservations = static_cast<int> (observedArrivalTimes.size());
        mTrainingTargetsMatrix.resize(2*nObservations, 0.0);
        for (int i = 0; i < nObservations; ++i)
        {
            double residual = observedArrivalTimes[i]
                            - predictedArrivalTimes[i];
            mTrainingTargetsMatrix[2*i]     = residual;
            mTrainingTargetsMatrix[2*i + 1] = weights[i]; 
        } 
    }
    /// Make a descriptor
    void makeDescriptor(const int nNeighbors)
    {
        mKDTreeDescriptor
            = std::make_unique<oneapi::dal::knn::descriptor
              <
                 double,
                 oneapi::dal::knn::method::kd_tree,
                 oneapi::dal::knn::task::search,
                 oneapi::dal::minkowski_distance::descriptor<double>>
              > (nNeighbors);
        mBruteForceDescriptor
            = std::make_unique<oneapi::dal::knn::descriptor
              <
                 double,
                 oneapi::dal::knn::method::brute_force,
                 oneapi::dal::knn::task::search,
                 oneapi::dal::minkowski_distance::descriptor<double>>
              > (nNeighbors);
    }
/*
    /// Load the model
    void load(const std::string &fileName)
    {
        if (!std::filesystem::exists(fileName))
        {
            throw std::invalid_argument("HDF5 file " + fileName
                                      + " does not exist");
        }
        auto file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        auto groupName = mNetwork + "." + mStation + "." + mPhase; 
        herr_t status;
        if (!H5Lexists(file, groupName.c_str(), H5P_DEFAULT))
        {
            H5Fclose(file);
            throw std::runtime_error(groupName + " station doesn't exist");
        }
        auto stationGroup = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
        if (!H5Lexists(stationGroup, KNN_GROUP_NAME, H5P_DEFAULT))
        {
            H5Gclose(stationGroup);
            H5Fclose(file);
            throw std::runtime_error(std::string{KNN_GROUP_NAME}
                                 + " KNN group doesn't exist");
        }
        auto knnGroup = H5Gopen(stationGroup, KNN_GROUP_NAME, H5P_DEFAULT);
        // Load the training features
        auto dataSet = H5Dopen2(knnGroup, FEATURES_NAME, H5P_DEFAULT);
        auto dataSpace = H5Dget_space(dataSet);
#ifndef NDEBUG
        assert(H5Sget_simple_extent_ndims(dataSpace) == 2);
#endif
        std::array<hsize_t, 2> dimensions;
        H5Sget_simple_extent_dims(dataSpace, dimensions.data(), nullptr);
        hsize_t length = 1;
        for (int i = 0; i < static_cast<int> (dimensions.size()); ++i)
        {
            length = length*dimensions.at(i);
        }
        mTrainingFeaturesMatrix.resize(length);
        auto memorySpace = H5Screate_simple(dimensions.size(),
                                            dimensions.data(), nullptr);
        status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memorySpace, dataSpace,
                         H5P_DEFAULT, mTrainingFeaturesMatrix.data());
        H5Sclose(memorySpace);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            throw std::runtime_error("Failed to read "
                                   + std::string{FEATURES_NAME});
        }
        // Load the training features
        dataSet = H5Dopen2(knnGroup, TARGETS_NAME, H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
#ifndef NDEBUG
        assert(H5Sget_simple_extent_ndims(dataSpace) == 2);
#endif
        H5Sget_simple_extent_dims(dataSpace, dimensions.data(), nullptr);
        length = 1;
        for (int i = 0; i < static_cast<int> (dimensions.size()); ++i)
        {
            length = length*dimensions.at(i);
        }
        mTrainingTargetsMatrix.resize(length);
        memorySpace = H5Screate_simple(dimensions.size(),
                                       dimensions.data(), nullptr);
        status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memorySpace, dataSpace,
                         H5P_DEFAULT, mTrainingTargetsMatrix.data());
        H5Sclose(memorySpace);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            throw std::runtime_error("Failed to read "
                                   + std::string{TARGETS_NAME});
        }
        // Load the scalars
        const std::array<size_t, 1> onesDimension{1};
        auto scalarSpace = H5Screate_simple(1, onesDimension.data(), nullptr);

        dataSet = H5Dopen2(knnGroup, NEIGHBORS_NAME, H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        int nNeighbors;
        status = H5Dread(dataSet, H5T_NATIVE_INT, scalarSpace, dataSpace,
                         H5P_DEFAULT, &nNeighbors);
        mNeighbors = nNeighbors;
        H5Sclose(dataSpace);
        H5Dclose(dataSet);

        dataSet = H5Dopen2(knnGroup, MAXIMUM_DISTANCE_NAME, H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, scalarSpace, dataSpace,
                         H5P_DEFAULT, &mMaximumDistance);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        H5Sclose(scalarSpace); 
        // Close the rest up
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        // Train the model
        train();
    }
*/
    /// Train the model
    void train()
    {
        auto nTrainingObservations
            = static_cast<int64_t> (mTrainingFeaturesMatrix.size()/mFeatures);
        if (nTrainingObservations < mNeighbors)
        {
            throw std::runtime_error("nObservations < nNeighbors");
        }
        auto nFeatures = static_cast<int64_t> (mFeatures);
        const auto XTrain
            = oneapi::dal::homogen_table::wrap<double>
              (mTrainingFeaturesMatrix.data(), nTrainingObservations, nFeatures);
        makeDescriptor(mNeighbors);
        if (mOptimizationAlgorithm == OptimizationAlgorithm::KDTree)
        {
            mKNNModel
               = oneapi::dal::train(*mKDTreeDescriptor, XTrain).get_model();
        }
        else
        {
            mKNNModel
               = oneapi::dal::train(*mBruteForceDescriptor, XTrain).get_model();
        }
        mHaveModel = true;
    }
    /// Predict
    void predict(double *corrections, double *dtdx, double *dtdy, double *dtdz)
    {
        if (!haveModel())
        {
            throw std::runtime_error("Model not yet set");
        } 
        // Initialize result
        auto nObservations
            = static_cast<int64_t> (mFeaturesMatrix.size()/mFeatures);
        // Doing depth?
        bool doDepth = true;
        if (mFeatures != 3){doDepth = false;}
        // Perform inference
#ifndef NDEBUG
        assert(mFeaturesMatrix.size() ==
               static_cast<size_t> (nObservations*mFeatures));
#endif
        const auto X
            = oneapi::dal::homogen_table::wrap<double>
              (mFeaturesMatrix.data(), nObservations, mFeatures);
        oneapi::dal::table trainingIndices;
        if (mOptimizationAlgorithm == OptimizationAlgorithm::KDTree)
        {
            const auto testResult
                = oneapi::dal::infer(*mKDTreeDescriptor, X, mKNNModel);
            trainingIndices = testResult.get_indices(); 
        }
        else
        {
            const auto testResult
                = oneapi::dal::infer(*mBruteForceDescriptor, X, mKNNModel);
            trainingIndices = testResult.get_indices();
        }
        //const auto testResult = oneapi::dal::infer(*mKDTreeDescriptor, X, mKNNModel);
        // Extract result
        //const auto &trainingIndices = testResult.get_indices();
#ifndef NDEBUG
        assert(trainingIndices.get_row_count() == nObservations);
        assert(trainingIndices.get_column_count() == mNeighbors);
#endif
        oneapi::dal::row_accessor<const double> accessor{trainingIndices};
        for (int64_t row = 0; row < nObservations; ++row)
        {
            double xi = mFeaturesMatrix.at(mFeatures*row + 0); 
            double yi = mFeaturesMatrix.at(mFeatures*row + 1); 
            double zi = 0;
            if (doDepth){zi = mFeaturesMatrix.at(mFeatures*row + 2);}
            const auto rowValues = accessor.pull({row, row + 1});
            // Extract
            bool allZeroWeights{true};
            double correction{0};
            double numerator{0};
            double denominator{0};
            double dsdxi{0};
            double dsdyi{0};
            double dsdzi{0};
            if (mEvaluationMethod == EvaluationMethod::WeightedAverage)
            {
                for (int j = 0; j < static_cast<int> (mNeighbors); ++j)
                {
                    auto trainingRow
                         = static_cast<int> (std::round(rowValues[j]));
                    auto iSrcTargets = 2*trainingRow;
                    double residual = mTrainingTargetsMatrix[iSrcTargets];
                    double weight = mTrainingTargetsMatrix[iSrcTargets + 1];
                    numerator = numerator + weight*residual;
                    denominator = denominator + weight;
                    if (weight > 0){allZeroWeights = false;}
                }
                if (!allZeroWeights){correction = numerator/denominator;}
            }
            else if (mEvaluationMethod == 
                     EvaluationMethod::MaximumDistanceWeighted)
            {
                for (int j = 0; j < static_cast<int> (mNeighbors); ++j)
                {
                    auto trainingRow
                         = static_cast<int> (std::round(rowValues[j]));
                    auto iSrcFeatures = mFeatures*trainingRow;
                    auto iSrcTargets = 2*trainingRow;
                    double residual = mTrainingTargetsMatrix[iSrcTargets];
                    double weight = mTrainingTargetsMatrix[iSrcTargets + 1];
                    double dx = xi - mTrainingFeaturesMatrix[iSrcFeatures];
                    double dy = yi - mTrainingFeaturesMatrix[iSrcFeatures + 1];
                    double dz = 0;
                    if (doDepth)
                    {
                        dz = zi - mTrainingFeaturesMatrix[iSrcFeatures + 2];
                    }
                    double distance = std::sqrt( dx*dx + dy*dy + dz*dz );
                    if (distance > mMaximumDistance)
                    {
                        weight = 0;
                    }
                    numerator = numerator + weight*residual;
                    denominator = denominator + weight;
                    if (weight > 0){allZeroWeights = false;}
                }
                if (!allZeroWeights){correction = numerator/denominator;}
            }
            else
            {
                double dSumRWDdxi{0};
                double dSumRWDdyi{0}; 
                double dSumRWDdzi{0};
                double dSumWd2dDxi{0};
                double dSumWd2dDyi{0};
                double dSumWd2dDzi{0};
                // Canononical form of derivative is:
                // d/dx [ \sum_{i} r_i w_i (x) / \sum_{i} w_i ]
                // = \sum_{i} r_i dw_i/dx / \sum_{i} w_i
                // - (\sum_{i} r_i w_i / (\sum_i{i} w_i)^2)*(\sum_{i} dw_i/dx)
                for (int j = 0; j < static_cast<int> (mNeighbors); ++j)
                {
                    auto trainingRow
                         = static_cast<int> (std::round(rowValues[j]));
                    auto iSrcFeatures = mFeatures*trainingRow;
                    auto iSrcTargets = 2*trainingRow;
                    double residual = mTrainingTargetsMatrix.at(iSrcTargets);
                    double weight = mTrainingTargetsMatrix.at(iSrcTargets + 1); 
                    double dx = xi - mTrainingFeaturesMatrix.at(iSrcFeatures);
                    double dy = yi - mTrainingFeaturesMatrix.at(iSrcFeatures + 1); 
                    double dz = 0;
                    if (doDepth)
                    {
                        dz = zi - mTrainingFeaturesMatrix[iSrcFeatures + 2]; 
                    }
                    // Don't let anything be closer than a micrometer for
                    // numerical reasons
                    double distance
                        = std::max(1.e-6, std::sqrt( dx*dx + dy*dy + dz*dz ));
                    double inverseDistance{0};
                    double inverseDistance3{0};
                    if (distance > mMaximumDistance)
                    {
                        weight = 0;
                        continue;
                    }
                    else
                    {
                        constexpr double one{1};
                        inverseDistance = one/distance;
                        inverseDistance3 = one/(distance*distance*distance);
                        allZeroWeights = false;
                    }
                    numerator = numerator + (weight*inverseDistance)*residual;
                    denominator = denominator + weight*inverseDistance;
                    double wdxid3 = weight*(dx*inverseDistance3);
                    double wdyid3 = weight*(dy*inverseDistance3);
                    double wdzid3 = weight*(dz*inverseDistance3);
                    dSumRWDdxi = dSumRWDdxi + residual*wdxid3;
                    dSumRWDdyi = dSumRWDdyi + residual*wdyid3;
                    dSumRWDdzi = dSumRWDdzi + residual*wdzid3;
                    dSumWd2dDxi = dSumWd2dDxi + wdxid3;
                    dSumWd2dDyi = dSumWd2dDyi + wdyid3;
                    dSumWd2dDzi = dSumWd2dDzi + wdzid3;
                } // Loop
                if (!allZeroWeights)
                {
                    correction = numerator/denominator;
                    double denominator2 = denominator*denominator;
                    dsdxi =-dSumRWDdxi/denominator
                          + (numerator/denominator2)*dSumWd2dDxi;
                    dsdyi =-dSumRWDdyi/denominator
                          + (numerator/denominator2)*dSumWd2dDyi;
                    dsdzi =-dSumRWDdzi/denominator
                          + (numerator/denominator2)*dSumWd2dDzi;
                }
            } // End picking job
            // Save results
            if (corrections != nullptr){corrections[row] = correction;}
            if (dtdx != nullptr){dtdx[row] = dsdxi;}
            if (dtdy != nullptr){dtdy[row] = dsdyi;}
            if (dtdz != nullptr){dtdz[row] = dsdzi;}
        } // Loop on observations
    } // End function
    /// Have a model?
    [[nodiscard]] bool haveModel() const noexcept
    {
        return mHaveModel;
    }
//private:
    std::unique_ptr
    <
        oneapi::dal::knn::descriptor
        <
            double,
            oneapi::dal::knn::method::kd_tree,
            oneapi::dal::knn::task::search,
            oneapi::dal::minkowski_distance::descriptor<double>
        >
    > mKDTreeDescriptor{nullptr};
    std::unique_ptr
    <    
        oneapi::dal::knn::descriptor
        <
            double,
            oneapi::dal::knn::method::brute_force,
            oneapi::dal::knn::task::search,
            oneapi::dal::minkowski_distance::descriptor<double>
        >
    > mBruteForceDescriptor{nullptr};

    std::string mNetwork;
    std::string mStation;
    std::string mPhase;
    std::vector<double> mTrainingFeaturesMatrix;
    std::vector<double> mFeaturesMatrix;
    std::vector<double> mTrainingTargetsMatrix;
    oneapi::dal::knn::model<oneapi::dal::knn::task::search> mKNNModel;
    double mMaximumDistance{35000}; // Contribution 0 beyond this distance (m)
    int64_t mNeighbors{5}; // Number of nearest neighbors
    int mFeatures{3};  // X, Y, Z
    EvaluationMethod mEvaluationMethod{EvaluationMethod::InverseDistanceWeighted};
    OptimizationAlgorithm mOptimizationAlgorithm{OptimizationAlgorithm::KDTree};
    bool mInitialized{false};
    bool mHaveModel{false};
};

/// Constructor
SourceSpecific::SourceSpecific() :
    pImpl(std::make_unique<SourceSpecificImpl> ())
{
}

/// Copy constructor
SourceSpecific::SourceSpecific(const SourceSpecific &correction)
{
    *this = correction;
}

/// Move constructor
SourceSpecific::SourceSpecific(SourceSpecific &&correction) noexcept
{
    *this = std::move(correction);
}

/// Move assignment
SourceSpecific& SourceSpecific::operator=(SourceSpecific &&correction) noexcept
{
    if (&correction == this){return *this;}
    pImpl = std::move(correction.pImpl);
    return *this;
}

/// Copy assignment
SourceSpecific& SourceSpecific::operator=(const SourceSpecific &correction)
{
    if (&correction == this){return *this;}
    pImpl = std::make_unique<SourceSpecificImpl> ();
    if (correction.pImpl->mInitialized)
    {
        pImpl->mNetwork = correction.pImpl->mNetwork;
        pImpl->mStation = correction.pImpl->mStation;
        pImpl->mPhase   = correction.pImpl->mPhase;
        pImpl->mTrainingFeaturesMatrix = correction.pImpl->mTrainingFeaturesMatrix;
        pImpl->mFeaturesMatrix = correction.pImpl->mFeaturesMatrix;;
        pImpl->mTrainingTargetsMatrix = correction.pImpl->mTrainingTargetsMatrix;
        pImpl->mMaximumDistance = correction.pImpl->mMaximumDistance;
        pImpl->mFeatures = correction.pImpl->mFeatures;
        pImpl->mNeighbors = correction.pImpl->mNeighbors;
        pImpl->mEvaluationMethod = correction.pImpl->mEvaluationMethod;
        pImpl->mInitialized = correction.pImpl->mInitialized;
        pImpl->mHaveModel = correction.pImpl->mHaveModel;
        pImpl->mOptimizationAlgorithm
             = correction.pImpl->mOptimizationAlgorithm;
        if (pImpl->mOptimizationAlgorithm == OptimizationAlgorithm::KDTree)
        {
            if (correction.pImpl->mKDTreeDescriptor != nullptr)
            {
                pImpl->makeDescriptor(pImpl->mNeighbors);
            }
        }
        else
        {
            if (correction.pImpl->mBruteForceDescriptor != nullptr)
            {
                pImpl->makeDescriptor(pImpl->mNeighbors);
            }
        }
        if (pImpl->mInitialized)
        {
            pImpl->mKNNModel = correction.pImpl->mKNNModel; 
        }
    }
    return *this;
}

/// Reset class
void SourceSpecific::clear() noexcept
{
    pImpl = std::make_unique<SourceSpecificImpl> ();
}

/// Destructor
SourceSpecific::~SourceSpecific() = default;

/// Identification
void SourceSpecific::setStationNameAndPhase(const std::string &network,
                                            const std::string &station,
                                            const std::string &phase)
{
    if (network.empty()){throw std::invalid_argument("Network cannot be empty");}
    if (station.empty()){throw std::invalid_argument("Station cannot be empty");}
    if (phase.empty()){throw std::invalid_argument("Phase cannot be empty");}
    pImpl->mNetwork = network;
    pImpl->mStation = station;
    pImpl->mPhase = phase;
}

std::string SourceSpecific::getNetwork() const
{
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    return pImpl->mNetwork;
}

std::string SourceSpecific::getStation() const
{
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    return pImpl->mStation;
}

std::string SourceSpecific::getPhase() const
{
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    return pImpl->mPhase;
}

bool SourceSpecific::haveStationNameAndPhase() const noexcept
{
    return !pImpl->mNetwork.empty();
}

/// Save the model
void SourceSpecific::save(const std::string &fileName) const
{
    if (!compiledWithHDF5()){throw std::runtime_error("Recompile with HDF5");}
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name/phase not set");
    }
    if (!haveModel())
    {
        throw std::runtime_error("Model not yet created");
    }
    if (fileName.empty()){throw std::invalid_argument("File name is empty");}
    auto path = std::filesystem::path(fileName).parent_path();
    if (!path.empty() && !std::filesystem::exists(path))
    {
        std::filesystem::create_directories(path);
        if (!std::filesystem::exists(path))
        {
            throw std::invalid_argument("Failed to create path: "
                                      + std::string(path));
        }
    }
#ifndef NDEBUG
    assert(pImpl->mTrainingFeaturesMatrix.size()%3 == 0); 
    assert(pImpl->mTrainingTargetsMatrix.size()%2 == 0); 
#endif
#ifdef WITH_HDF5
    hid_t file;
    if (std::filesystem::exists(fileName))
    {
        file = H5Fopen(fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else
    {
        file = H5Fcreate(fileName.c_str(), H5F_ACC_TRUNC,
                         H5P_DEFAULT, H5P_DEFAULT);
    }
    if (file == H5I_INVALID_HID)
    {
        throw std::runtime_error("Failed to open file: " + fileName);
    }
    // Create a station group
    auto groupName = getNetwork() + "." + getStation() + "." + getPhase();
    herr_t status;
    hid_t stationGroup;
    if (!H5Lexists(file, groupName.c_str(), H5P_DEFAULT))
    {
        stationGroup = H5Gcreate2(file, groupName.c_str(),
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (stationGroup == H5I_INVALID_HID)
        {
            H5Fclose(file);
            throw std::runtime_error("Could not station group " + groupName);
        }   
    }
    else
    {
        stationGroup = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
    }
    // Create a KNN group in the station group
    hid_t knnGroup;
    if (!H5Lexists(stationGroup, KNN_GROUP_NAME, H5P_DEFAULT))
    {
        knnGroup = H5Gcreate2(stationGroup, KNN_GROUP_NAME,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (knnGroup == H5I_INVALID_HID)
        {
            H5Gclose(stationGroup);
            H5Fclose(file);
            throw std::runtime_error("Could not KNN group "
                                   + std::string(KNN_GROUP_NAME));
        }
    }
    else
    { 
        knnGroup = H5Gopen(stationGroup, KNN_GROUP_NAME, H5P_DEFAULT);
    }
    if (knnGroup == H5I_INVALID_HID)
    {
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Could not open KNN group");
    }
    // Data sizes
    hsize_t nObservations = pImpl->mTrainingFeaturesMatrix.size()/3;
    constexpr hsize_t nFeatures{3};
    std::array<hsize_t, 2> featuresDimensions{nObservations, nFeatures};
    std::array<hsize_t, 2> targetsDimensions{nObservations, 2};
#ifndef NDEBUG
    assert(nObservations*2 == pImpl->mTrainingTargetsMatrix.size());
    assert(nObservations*nFeatures == pImpl->mTrainingFeaturesMatrix.size());
#endif
    // Ensure datasets do not exist
    if (H5Lexists(knnGroup, FEATURES_NAME, H5P_DEFAULT))
    {
        status = H5Ldelete(knnGroup, FEATURES_NAME, H5P_DEFAULT);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            std::runtime_error("Failed to delete "
                             + std::string {FEATURES_NAME});
        }
    }
    if (H5Lexists(knnGroup, TARGETS_NAME, H5P_DEFAULT))
    {
        status = H5Ldelete(knnGroup, TARGETS_NAME, H5P_DEFAULT);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            std::runtime_error("Failed to delete "
                             + std::string {TARGETS_NAME});
        }
    }
    // Write features
    auto memorySpace = H5Screate_simple(2, featuresDimensions.data(), nullptr);
    if (memorySpace == H5I_INVALID_HID)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to create features memory space");
    }
    auto dataSet = H5Dcreate2(knnGroup, FEATURES_NAME,
                              H5T_NATIVE_DOUBLE, memorySpace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataSet == H5I_INVALID_HID)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to create features data space");
    }
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, memorySpace,
                      H5P_DEFAULT, pImpl->mTrainingFeaturesMatrix.data());
    if (status < 0)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to write features");
    }
    H5Dclose(dataSet);
    H5Sclose(memorySpace);
    // Write targets
    memorySpace = H5Screate_simple(2, targetsDimensions.data(), nullptr);
    if (memorySpace == H5I_INVALID_HID)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to create targets memory space");
    }
    dataSet = H5Dcreate2(knnGroup, TARGETS_NAME,
                         H5T_NATIVE_DOUBLE, memorySpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataSet == H5I_INVALID_HID)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to create targets data space");
    }
    status = H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, memorySpace,
                      H5P_DEFAULT, pImpl->mTrainingTargetsMatrix.data());
    if (status < 0)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to write targets");
    }
    H5Dclose(dataSet);
    H5Sclose(memorySpace);
    // Write scalars
    if (H5Lexists(knnGroup, NEIGHBORS_NAME, H5P_DEFAULT))
    {
        status = H5Ldelete(knnGroup, NEIGHBORS_NAME, H5P_DEFAULT);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            std::runtime_error("Failed to delete "
                             + std::string {NEIGHBORS_NAME});
        }
    }
    if (H5Lexists(knnGroup, EVALUATION_METHOD_NAME, H5P_DEFAULT))
    {
        status = H5Ldelete(knnGroup, EVALUATION_METHOD_NAME, H5P_DEFAULT);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            std::runtime_error("Failed to delete "
                             + std::string {EVALUATION_METHOD_NAME});
        }
    }
    if (H5Lexists(knnGroup, MAXIMUM_DISTANCE_NAME, H5P_DEFAULT))
    {
        status = H5Ldelete(knnGroup, MAXIMUM_DISTANCE_NAME, H5P_DEFAULT);
        if (status < 0)
        {
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            std::runtime_error("Failed to delete "
                             + std::string {MAXIMUM_DISTANCE_NAME});
        }
    }
    const std::array<size_t, 1> onesDimension{1};
    auto scalarSpace = H5Screate_simple(1, onesDimension.data(), nullptr);
    dataSet = H5Dcreate2(knnGroup, NEIGHBORS_NAME,
                         H5T_NATIVE_INT, scalarSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int nNeighbors = getNumberOfNeighbors();
    H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL, scalarSpace,
             H5P_DEFAULT, &nNeighbors);
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(knnGroup, EVALUATION_METHOD_NAME,
                         H5T_NATIVE_INT, scalarSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int evaluationMethod = static_cast<int> (getEvaluationMethod());
    H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL, scalarSpace,
             H5P_DEFAULT, &evaluationMethod);
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(knnGroup, OPTIMIZATION_ALGORITHM_NAME,
                         H5T_NATIVE_INT, scalarSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int optimizationMethod = static_cast<int> (getOptimizationAlgorithm());
    H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL, scalarSpace,
             H5P_DEFAULT, &optimizationMethod);
    H5Dclose(dataSet);

    dataSet = H5Dcreate2(knnGroup, MAXIMUM_DISTANCE_NAME,
                         H5T_NATIVE_DOUBLE, scalarSpace,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    double maximumDistance{-1};
    if (haveMaximumDistance()){maximumDistance = getMaximumDistance();}
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, scalarSpace,
             H5P_DEFAULT, &maximumDistance);
    H5Dclose(dataSet);
    H5Sclose(scalarSpace);
    // Close everything else
    H5Gclose(knnGroup);
    H5Gclose(stationGroup);
    H5Fclose(file);
#endif
}

/// Load the model
void SourceSpecific::load(const std::string &fileName)
{
    if (!compiledWithHDF5()){throw std::runtime_error("Recompile with HDF5");}
    if (!haveStationNameAndPhase())
    {   
        throw std::runtime_error("Station name/phase not set");
    }
    if (!std::filesystem::exists(fileName))
    {
        throw std::runtime_error(fileName + " does not exist");
    }
#ifdef WITH_HDF5
    auto file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    auto groupName = getNetwork() + "." + getStation() + "." + getPhase();
    herr_t status;
    if (!H5Lexists(file, groupName.c_str(), H5P_DEFAULT))
    {
        H5Fclose(file);
        throw std::runtime_error(groupName + " station doesn't exist");
    }
    auto stationGroup = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
    if (!H5Lexists(stationGroup, KNN_GROUP_NAME, H5P_DEFAULT))
    {
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error(std::string{KNN_GROUP_NAME}
                             + " group doesn't exist");
    }
    auto knnGroup = H5Gopen(stationGroup, KNN_GROUP_NAME, H5P_DEFAULT);
    // Load the training features
    auto dataSet = H5Dopen2(knnGroup, FEATURES_NAME, H5P_DEFAULT);
    auto dataSpace = H5Dget_space(dataSet);
#ifndef NDEBUG
    assert(H5Sget_simple_extent_ndims(dataSpace) == 2);
#endif
    std::array<hsize_t, 2> dimensions;
    H5Sget_simple_extent_dims(dataSpace, dimensions.data(), nullptr);
    hsize_t length = 1;
    for (int i = 0; i < static_cast<int> (dimensions.size()); ++i)
    {
        length = length*dimensions.at(i);
    }
    pImpl->mTrainingFeaturesMatrix.resize(length, 0);
    auto memorySpace = H5Screate_simple(dimensions.size(),
                                        dimensions.data(), nullptr);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memorySpace, dataSpace,
                     H5P_DEFAULT, pImpl->mTrainingFeaturesMatrix.data());
    H5Sclose(memorySpace);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    if (status < 0)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to read "
                               + std::string{FEATURES_NAME});
    }
    // Load the training features
    dataSet = H5Dopen2(knnGroup, TARGETS_NAME, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
#ifndef NDEBUG
    assert(H5Sget_simple_extent_ndims(dataSpace) == 2);
#endif
    H5Sget_simple_extent_dims(dataSpace, dimensions.data(), nullptr);
    length = 1;
    for (int i = 0; i < static_cast<int> (dimensions.size()); ++i)
    {
        length = length*dimensions.at(i);
    }
    pImpl->mTrainingTargetsMatrix.resize(length);
    memorySpace = H5Screate_simple(dimensions.size(),
                                   dimensions.data(), nullptr);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, memorySpace, dataSpace,
                     H5P_DEFAULT, pImpl->mTrainingTargetsMatrix.data());
    H5Sclose(memorySpace);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    if (status < 0)
    {
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to read "
                               + std::string{TARGETS_NAME});
    }
    // Load the scalars
    const std::array<size_t, 1> onesDimension{1};
    auto scalarSpace = H5Screate_simple(1, onesDimension.data(), nullptr);

    dataSet = H5Dopen2(knnGroup, NEIGHBORS_NAME, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    int nNeighbors;
    status = H5Dread(dataSet, H5T_NATIVE_INT, scalarSpace, dataSpace,
                     H5P_DEFAULT, &nNeighbors);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    if (status < 0)
    {
        H5Sclose(scalarSpace);
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to read number of neighbors");
    }
    setNumberOfNeighbors(nNeighbors);

    dataSet = H5Dopen2(knnGroup, EVALUATION_METHOD_NAME, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    int iEvaluationMethod;
    status = H5Dread(dataSet, H5T_NATIVE_INT, scalarSpace, dataSpace,
                     H5P_DEFAULT, &iEvaluationMethod);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    if (status < 0)
    {
        H5Sclose(scalarSpace);
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to read evaluation method");
    }
    setEvaluationMethod(static_cast<EvaluationMethod> (iEvaluationMethod));

    if (H5Lexists(knnGroup, OPTIMIZATION_ALGORITHM_NAME, H5P_DEFAULT))
    {
        dataSet = H5Dopen2(knnGroup, OPTIMIZATION_ALGORITHM_NAME, H5P_DEFAULT);
        dataSpace = H5Dget_space(dataSet);
        int iOptimizationAlgorithm{0};
        status = H5Dread(dataSet, H5T_NATIVE_INT, scalarSpace, dataSpace,
                         H5P_DEFAULT, &iOptimizationAlgorithm);
        H5Sclose(dataSpace);
        H5Dclose(dataSet);
        if (status < 0) 
        {    
            H5Sclose(scalarSpace);
            H5Gclose(knnGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            throw std::runtime_error("Failed to read optimization algorithm");
        }    
        setOptimizationAlgorithm(static_cast<OptimizationAlgorithm>
                                 (iOptimizationAlgorithm));
    }

    double maximumDistance{-1};
    dataSet = H5Dopen2(knnGroup, MAXIMUM_DISTANCE_NAME, H5P_DEFAULT);
    dataSpace = H5Dget_space(dataSet);
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, scalarSpace, dataSpace,
                     H5P_DEFAULT, &maximumDistance);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    if (status < 0)
    {
        H5Sclose(scalarSpace);
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Failed to read maximum distance");
    }
    pImpl->mMaximumDistance =-1;
    if (maximumDistance >= 0){setMaximumDistance(maximumDistance);}

    H5Sclose(scalarSpace);
    // Close the rest up
    H5Gclose(knnGroup);
    H5Gclose(stationGroup);
    H5Fclose(file);
    // And train it
    pImpl->train();
#ifndef NDEBUG
    assert(haveModel());
#endif
#endif
}

/// Evaluate
double SourceSpecific::evaluate(
     const double xSource,
     const double ySource,
     const double zSource,
     const double predictedTime) const
{
    return evaluate(xSource, ySource, zSource,
                    nullptr, nullptr, nullptr, nullptr,
                    predictedTime);
}

double SourceSpecific::evaluate(
     const double xSource,
     const double ySource,
     const double zSource,
     double *dtdt0, double *dtdx, double *dtdy, double *dtdz,
     const double predictedTime) const
{
    if (dtdt0 != nullptr){*dtdt0 = 0;} // Always zero
    std::vector<double> xSources{xSource};
    std::vector<double> ySources{ySource};
    std::vector<double> zSources{zSource};
    constexpr bool isTraining{false};
    pImpl->fillFeatureVector(xSources, ySources, zSources, isTraining);
    double correction{0};
    pImpl->predict(&correction, dtdx, dtdy, dtdz);
    return predictedTime + correction;
}
 
/// Train
void SourceSpecific::train(
    const std::vector<double> &xSources,
    const std::vector<double> &ySources,
    const std::vector<double> &zSources,
    const std::vector<double> &observedArrivalTimes,
    const std::vector<double> &predictedArrivalTimes,
    const std::vector<double> &weights)
{
    auto nNeighbors = getNumberOfNeighbors();
    if (static_cast<int> (observedArrivalTimes.size()) < nNeighbors)
    {
        throw std::invalid_argument("At least " + std::to_string(nNeighbors)
                                  + " observations required");
    }
    if (observedArrivalTimes.size() != xSources.size())
    {   
        throw std::invalid_argument("Observed/xSources sizes inconsistent");
    }   
    if (observedArrivalTimes.size() != ySources.size())
    {
        throw std::invalid_argument("Observed/ySources sizes inconsistent");
    } 
    if (observedArrivalTimes.size() != zSources.size())
    {
        throw std::invalid_argument("Observed/zSources sizes inconsistent");
    }
    if (observedArrivalTimes.size() != predictedArrivalTimes.size())
    {
        throw std::invalid_argument("Observed/predicted sizes inconsistent");
    }
    if (observedArrivalTimes.size() != weights.size())
    {
        throw std::invalid_argument("Observed/weights sizes inconsistent");
    }
    if (nNeighbors > observedArrivalTimes.size())
    {
        throw std::invalid_argument("More neighbors than observations");
    }
    constexpr bool isTraining{true};
    pImpl->fillFeatureVector(xSources, ySources, zSources, isTraining);
    pImpl->fillTargetVector(observedArrivalTimes,
                            predictedArrivalTimes, weights);
    pImpl->mNeighbors = nNeighbors;
    pImpl->train();
}

bool SourceSpecific::haveModel() const noexcept
{
    return pImpl->haveModel();
}

/// Evaluation method
void SourceSpecific::setEvaluationMethod(
    const EvaluationMethod method) noexcept
{
    pImpl->mEvaluationMethod = method;
}

SourceSpecific::EvaluationMethod
SourceSpecific::getEvaluationMethod() const noexcept
{
    return pImpl->mEvaluationMethod;
}

/// Number of neighbors
void SourceSpecific::setNumberOfNeighbors(const int nNeighbors)
{
    if (nNeighbors < 1)
    {
        throw std::invalid_argument("Number of neighbors must be positive");
    }
    pImpl->mNeighbors = nNeighbors;
}

int SourceSpecific::getNumberOfNeighbors() const noexcept
{
    return static_cast<int> (pImpl->mNeighbors);
}

/// Optimization algorithm
void SourceSpecific::setOptimizationAlgorithm(
    const OptimizationAlgorithm algorithm) noexcept
{
    pImpl->mOptimizationAlgorithm = algorithm;
}

SourceSpecific::OptimizationAlgorithm 
    SourceSpecific::getOptimizationAlgorithm() const noexcept
{
    return pImpl->mOptimizationAlgorithm;
}

/// Maximum distance
void SourceSpecific::setMaximumDistance(const double maximumDistance) 
{
    if (maximumDistance <= 0)
    {
        throw std::invalid_argument("Maximum distance must be positive");
    }
    pImpl->mMaximumDistance = maximumDistance;
}

double SourceSpecific::getMaximumDistance() const
{
    if (!haveMaximumDistance())
    {
        throw std::runtime_error("Maximum distance not set");
    }
    return pImpl->mMaximumDistance;
}

bool SourceSpecific::haveMaximumDistance() const noexcept
{
    return (pImpl->mMaximumDistance > 0);
}

/// HDF5 support?
bool SourceSpecific::compiledWithHDF5() noexcept
{
#ifdef WITH_HDF5
    return true;
#else
    return false;
#endif
}

