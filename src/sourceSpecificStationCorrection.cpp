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
#include <daal.h>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uLocator/sourceSpecificStationCorrection.hpp"
#include "uLocator/position/wgs84.hpp"

#define KNN_GROUP_NAME "KNN"
#define FEATURES_NAME "KNNTrainingFeaturesXYZ"
#define TARGETS_NAME  "KNNTargetsResidualsWeights"
#define NEIGHBORS_NAME "NumberOfNeighbors"
#define MAXIMUM_DISTANCE_NAME "MaximumDistance"
#define UTM_ZONE_NAME "UTMZone"

using namespace ULocator;

class SourceSpecificStationCorrection::SourceSpecificStationCorrectionImpl
{
public:
    void fillFeatureVector(
        const std::vector<double> &latitudes,
        const std::vector<double> &longitudes,
        const std::vector<double> &depths,
        const bool isTraining)
    {
        bool doDepth = true;
        if (mFeatures != 3){doDepth = false;}
#ifndef NDEBUG
        assert(latitudes.size() == longitudes.size());
        if (doDepth){assert(latitudes.size() == depths.size());}
#endif
        auto nObservations = static_cast<int> (latitudes.size());
        ULocator::Position::WGS84 position;
        position.setAlternateUTMZone(mUTMZone);
        std::vector<double> featuresMatrix(mFeatures*nObservations);
        for (int i = 0; i < nObservations; ++i)
        {
            auto i0 = mFeatures*i;
            position.setLatitudeAndLongitude(latitudes[i], longitudes[i]); 
            auto x = position.getEasting();
            auto y = position.getNorthing();
            featuresMatrix[i0] = x;
            featuresMatrix[i0 + 1] = y;
            if (doDepth)
            {
                featuresMatrix[i0 + 2] = depths[i];
            }
        }
/*
for (int i = 0; i < nObservations; ++i)
{
 std::cout << featuresMatrix[3*i] << " " << featuresMatrix[3*i + 1] << " " << featuresMatrix[3*i + 2] << std::endl;
} 
*/
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
        mDescriptor
            = std::make_unique<oneapi::dal::knn::descriptor
              <
                 double,
                 oneapi::dal::knn::method::kd_tree,
                 oneapi::dal::knn::task::search,
                 oneapi::dal::minkowski_distance::descriptor<double>>
              > (nNeighbors);
    }
    /// Save the model
    void save(const std::string &fileName)
    {
#ifndef NDEBUG
        assert(mTrainingFeaturesMatrix.size()%3 == 0);
        assert(mTrainingTargetsMatrix.size()%2 == 0);
#endif
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
        auto groupName = mNetwork + "." + mStation + "." + mPhase; 
        herr_t status;
        hid_t stationGroup;
        if (!H5Lexists(file, groupName.c_str(), H5P_DEFAULT))
        {
            stationGroup = H5Gcreate2(file, groupName.c_str(),
                                      H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (stationGroup == H5I_INVALID_HID)
            {
                 H5Fclose(file);
                 throw std::runtime_error("Could not station group "
                                        + groupName);
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
        hsize_t nObservations = mTrainingFeaturesMatrix.size()/3;
        constexpr hsize_t nFeatures{3};
        std::array<hsize_t, 2> featuresDimensions{nObservations, nFeatures};
        std::array<hsize_t, 2> targetsDimensions{nObservations, 2};
#ifndef NDEBUG
        assert(nObservations*2 == mTrainingTargetsMatrix.size());
        assert(nObservations*nFeatures == mTrainingFeaturesMatrix.size());
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
                          H5P_DEFAULT, mTrainingFeaturesMatrix.data());
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
                          H5P_DEFAULT, mTrainingTargetsMatrix.data());
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
        if (H5Lexists(knnGroup, UTM_ZONE_NAME, H5P_DEFAULT))
        {
            status = H5Ldelete(knnGroup, UTM_ZONE_NAME, H5P_DEFAULT);
            if (status < 0)
            {
                H5Gclose(knnGroup);
                H5Gclose(stationGroup);
                H5Fclose(file);
                std::runtime_error("Failed to delete "
                                 + std::string {UTM_ZONE_NAME});
            }   
        }
        const std::array<size_t, 1> onesDimension{1};
        auto scalarSpace = H5Screate_simple(1, onesDimension.data(), nullptr);
        dataSet = H5Dcreate2(knnGroup, NEIGHBORS_NAME,
                             H5T_NATIVE_INT, scalarSpace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        auto nNeighbors = static_cast<int> (mNeighbors);
        H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL, scalarSpace,
                 H5P_DEFAULT, &nNeighbors);
        H5Dclose(dataSet);
        dataSet = H5Dcreate2(knnGroup, MAXIMUM_DISTANCE_NAME,
                             H5T_NATIVE_DOUBLE, scalarSpace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, scalarSpace,
                 H5P_DEFAULT, &mMaximumDistance);
        H5Dclose(dataSet); 
        dataSet = H5Dcreate2(knnGroup, UTM_ZONE_NAME,
                             H5T_NATIVE_INT, scalarSpace,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Dwrite(dataSet, H5T_NATIVE_INT, H5S_ALL, scalarSpace,
                 H5P_DEFAULT, &mUTMZone); 
        H5Dclose(dataSet);
        H5Sclose(scalarSpace);
        // Close everything else
        H5Gclose(knnGroup);
        H5Gclose(stationGroup);
        H5Fclose(file); 
    }
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
        mKNNModel = oneapi::dal::train(*mDescriptor, XTrain).get_model();
        mHaveModel = true;
    }
    /// Predict
    std::vector<double> predict(const EvaluationMethod method) const
    {
        if (!haveModel())
        {
            throw std::runtime_error("Model not yet set");
        } 
        // Initialize result
        auto nObservations
            = static_cast<int64_t> (mFeaturesMatrix.size()/mFeatures);
        std::vector<double> knnResiduals(nObservations, 0.0);
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
        const auto testResult = oneapi::dal::infer(*mDescriptor, X, mKNNModel);
        // Extract result
        const auto &trainingIndices = testResult.get_indices();
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
            double numerator = 0;
            double denominator = 0;
            for (int64_t j = 0; j < mNeighbors; ++j)
            {
                auto trainingRow = static_cast<int> (std::round(rowValues[j]));
                auto iSrcFeatures = mFeatures*trainingRow;
                auto iSrcTargets = 2*trainingRow;
                double dx = xi - mTrainingFeaturesMatrix.at(iSrcFeatures);
                double dy = yi - mTrainingFeaturesMatrix.at(iSrcFeatures + 1);
                double dz = 0; 
                if (doDepth)
                {
                    dz = zi - mTrainingFeaturesMatrix.at(iSrcFeatures + 2);
                }
                double distance = std::sqrt( dx*dx + dy*dy + dz*dz );
                double residual = mTrainingTargetsMatrix.at(iSrcTargets);
                double weight   = mTrainingTargetsMatrix.at(iSrcTargets + 1);
                if (method == EvaluationMethod::MaximumDistanceWeighted)
                {
                    if (distance > mMaximumDistance)
                    {
                        weight = 0;
                    }
                }
                else if (method == EvaluationMethod::InverseDistanceWeighted)
                {
                    double distanceWeight = 1./std::max(0.1, distance);
                    weight = distanceWeight*weight;
                    if (distance > mMaximumDistance){weight = 0;}
                }
                numerator = numerator + weight*residual;
                denominator = denominator + weight;
                if (weight > 0){allZeroWeights = false;}
            }
            knnResiduals[row] = 0;
            if (!allZeroWeights)
            {
                knnResiduals[row] = numerator/denominator;
            }
        }
        return knnResiduals;
    }
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
    > mDescriptor{nullptr};
    std::string mNetwork;
    std::string mStation;
    std::string mPhase;
    std::vector<double> mTrainingFeaturesMatrix;
    std::vector<double> mFeaturesMatrix;
    std::vector<double> mTrainingTargetsMatrix;
    oneapi::dal::knn::model<oneapi::dal::knn::task::search> mKNNModel;
    double mMaximumDistance{35000}; // Contribution 0 beyond this distance (m)
    int mUTMZone{12};  // By default works for Utah/Yellowstone
    int mFeatures{3};  // UTM X / UTM Y / Depth
    int64_t mNeighbors{5}; // Number of nearest neighbors
    bool mInitialized{false};
    bool mHaveModel{false};
};

/// Constructor
SourceSpecificStationCorrection::SourceSpecificStationCorrection() :
    pImpl(std::make_unique<SourceSpecificStationCorrectionImpl> ())
{
}

/// Move constructor
SourceSpecificStationCorrection::SourceSpecificStationCorrection(
    SourceSpecificStationCorrection &&correction) noexcept
{
    *this = std::move(correction);
}

/// Move assignment
SourceSpecificStationCorrection& SourceSpecificStationCorrection::operator=(
    SourceSpecificStationCorrection &&correction) noexcept
{
    if (&correction == this){return *this;}
    pImpl = std::move(correction.pImpl);
    return *this;
}

/// Destructor
SourceSpecificStationCorrection::~SourceSpecificStationCorrection() = default;

/// Initialize model
void SourceSpecificStationCorrection::initialize(const std::string &network,
                                                 const std::string &station,
                                                 const std::string &phase)
{
    if (network.empty()){throw std::invalid_argument("Network cannot be empty");}  
    if (station.empty()){throw std::invalid_argument("Station cannot be empty");}
    if (phase.empty()){throw std::invalid_argument("Phase cannot be empty");}
    pImpl->mNetwork = network;
    pImpl->mStation = station;
    pImpl->mPhase = phase;
    pImpl->mInitialized = true;
}

/// Initialized?
bool SourceSpecificStationCorrection::isInitialized() const noexcept
{
    return pImpl->mInitialized;
}

/// Save the model
void SourceSpecificStationCorrection::save(const std::string &fileName) const
{
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
    pImpl->save(fileName); 
}

/// Load the model
void SourceSpecificStationCorrection::load(const std::string &fileName)
{
    if (!std::filesystem::exists(fileName))
    {
        throw std::runtime_error(fileName + " does not exist");
    }
    pImpl->load(fileName);
}

/// Evaluate
double SourceSpecificStationCorrection::evaluate(
     const double latitude,
     const double longitude,
     const double depth,
     const EvaluationMethod method) const
{
    std::vector<double> latitudes{latitude};
    std::vector<double> longitudes{longitude};
    std::vector<double> depths{depth};
    return evaluate(latitudes, longitudes, depths, method)[0];
}

std::vector<double>
     SourceSpecificStationCorrection::evaluate(
     const std::vector<double> &latitudes,
     const std::vector<double> &longitudes,
     const std::vector<double> &depths,
     const EvaluationMethod method) const
{
    std::vector<double> result;
    if (!haveModel()){throw std::runtime_error("Model not set");}
    if (latitudes.size() != depths.size())
    {
        throw std::invalid_argument("latitudes/depths size unequal");
    }
    if (longitudes.size() != depths.size())
    {
        throw std::invalid_argument("longitudes/depths size unequal");
    }
    if (depths.empty()){return result;}
    constexpr bool isTraining{false};
    pImpl->fillFeatureVector(latitudes, longitudes, depths, isTraining);
    return pImpl->predict(method);
}

/// Train
void SourceSpecificStationCorrection::train(
    const std::vector<double> &latitudes,
    const std::vector<double> &longitudes,
    const std::vector<double> &depths,
    const std::vector<double> &observedArrivalTimes,
    const std::vector<double> &predictedArrivalTimes,
    const std::vector<double> &weights,
    const int nNeighbors,
    const double maximumDistance)
{
    if (nNeighbors < 1)
    {
        throw std::invalid_argument("Number of neighbors must be positive");
    }
    if (static_cast<int> (observedArrivalTimes.size()) < nNeighbors)
    {
        throw std::invalid_argument("At least " + std::to_string(nNeighbors)
                                  + " observations required");
    }
    if (observedArrivalTimes.size() != latitudes.size())
    {   
        throw std::invalid_argument("Observed/latitudes sizes inconsistent");
    }   
    if (observedArrivalTimes.size() != longitudes.size())
    {
        throw std::invalid_argument("Observed/longitudes sizes inconsistent");
    } 
    if (observedArrivalTimes.size() != depths.size())
    {
        throw std::invalid_argument("Observed/depths sizes inconsistent");
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
    if (maximumDistance > 0)
    {
        pImpl->mMaximumDistance = maximumDistance;
    }
    constexpr bool isTraining{true};
    pImpl->fillFeatureVector(latitudes, longitudes, depths, isTraining);
    pImpl->fillTargetVector(observedArrivalTimes,
                            predictedArrivalTimes, weights);
    pImpl->mNeighbors = nNeighbors;
    pImpl->train();
}

bool SourceSpecificStationCorrection::haveModel() const noexcept
{
    return pImpl->haveModel();
}

/// Number of neighbors
int SourceSpecificStationCorrection::getNumberOfNeighbors() const
{
    return pImpl->mNeighbors;
}

/// Maximum distance
double SourceSpecificStationCorrection::getMaximumDistance() const
{
    return pImpl->mMaximumDistance;
}
