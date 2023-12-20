#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#ifdef WITH_HDF5
#include <hdf5.h>
#endif
#include "uLocator/corrections/static.hpp"
#include "weightedMedian.hpp"
#include "weightedMean.hpp"

#define STATIC_CORRECTION_GROUP_NAME "StaticCorrection"
#define STATIC_CORRECTION_NAME "Correction"

using namespace ULocator::Corrections;

class Static::StaticImpl
{
public:
    std::string mNetwork;
    std::string mStation;
    std::string mPhase;
    double mCorrection{0}; 
    bool mHaveCorrection{false};
};

/// Constructor
Static::Static() :
    pImpl(std::make_unique<StaticImpl> ())
{
}

/// Copy constructor
Static::Static(const Static &correction)
{
    *this = correction;
}

/// Move constructor
Static::Static(Static &&correction) noexcept
{
    *this = std::move(correction);
}

/// Copy assignment
Static& Static::operator=(const Static &correction)
{
    if (&correction == this){return *this;}
    pImpl = std::make_unique<StaticImpl> (*correction.pImpl);
    return *this;
}

/// Move assignment
Static& Static::operator=(Static &&correction) noexcept
{
    if (&correction == this){return *this;}
    pImpl = std::move(correction.pImpl);
    return *this;
}

/// Destructor
Static::~Static() = default;

/// Initialization
void Static::setStationNameAndPhase(const std::string &network,
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

std::string Static::getNetwork() const
{
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    return pImpl->mNetwork;
}

std::string Static::getStation() const
{
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    return pImpl->mStation;
}

std::string Static::getPhase() const
{
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    return pImpl->mPhase;
}

bool Static::haveStationNameAndPhase() const noexcept
{
    return !pImpl->mNetwork.empty();
}

void Static::setCorrection(const double correction) noexcept
{
    pImpl->mCorrection = correction;
    pImpl->mHaveCorrection = true;
}

double Static::getCorrection() const
{
    constexpr double zero{0};
    return evaluate(zero);
}

bool Static::haveCorrection() const noexcept
{
    return pImpl->mHaveCorrection;
}

void Static::train(const std::vector<double> &observedArrivalTimes,
                   const std::vector<double> &predictedArrivalTimes,
                   const Method method)
{
    constexpr double one{1};
    std::vector<double> weights(observedArrivalTimes.size(), one);
    train(observedArrivalTimes, predictedArrivalTimes, weights, method);
}

void Static::train(const std::vector<double> &observedArrivalTimes,
                   const std::vector<double> &predictedArrivalTimes,
                   const std::vector<double> &weights,
                   const Method method)
{
    if (observedArrivalTimes.empty())
    {
        throw std::invalid_argument("No observations");
    }
    if (observedArrivalTimes.size() != predictedArrivalTimes.size())
    {
        throw std::invalid_argument("Inconsistent observed/predicted sizes");
    }
    if (observedArrivalTimes.size() != weights.size())
    {
        throw std::invalid_argument("Inconsistent observed/weights sizes");
    }
    auto nObservations = static_cast<int> (observedArrivalTimes.size());
    // Check the weights
    bool allZero{true};
    for (const auto &weight : weights)
    {
        if (weight < 0)
        {
            throw std::invalid_argument("All weights must be positive");
        }
        if (weight > 0){allZero = false;}
    }
    if (allZero){throw std::invalid_argument("Weights cannot all be zero");}
    // Tabulate residuals
    std::vector<double> residuals(nObservations);
        std::transform(observedArrivalTimes.begin(), observedArrivalTimes.end(),
                       predictedArrivalTimes.begin(), residuals.begin(),
                       std::minus<double> ());
    // Basically, I am going to compute
    //      obs - est = bias
    // so as to learn the systematic bias.
    // Say this bias is 1, i.e.,
    //    obs - est = 2 - 1 = 1
    // Now, to remove the bias I do the following:
    //    obs - est - bias = obs - (est + bias) = 0
    // Hence, the correction is simply the bias.  All that is left to do is to
    // now define the bias.
    if (method == Static::Method::Median)
    {
        std::vector<std::pair<double, int>> workSpace(nObservations);
        auto medianCorrection
            = ::weightedMedian(residuals, weights, workSpace);
        setCorrection(medianCorrection);
    }
    else
    {
        auto meanCorrection = ::weightedMean(residuals, weights);
        setCorrection(meanCorrection);
    }
}

/// Applies the model
double Static::evaluate(const double predictedTime) const
{
    if (!haveCorrection()){throw std::runtime_error("Correction not set");}
    return evaluate(predictedTime, nullptr, nullptr, nullptr, nullptr);
}

double Static::evaluate(
    const double predictedTime,
    double *dtdt0, double *dtdx, double *dtdy, double *dtdz) const
{
    if (!haveCorrection()){throw std::runtime_error("Correction not set");}
    if (dtdt0 != nullptr){*dtdt0 = 0;}
    if (dtdx  != nullptr){*dtdx  = 0;}
    if (dtdy  != nullptr){*dtdy  = 0;}
    if (dtdz  != nullptr){*dtdz  = 0;}
    return predictedTime + pImpl->mCorrection;
}

double Static::operator()(const double predictedTime) const
{
    return evaluate(predictedTime);
}

double Static::operator()(
    const double predictedTime,
    double *dtdt0, double *dtdx, double *dtdy, double *dtdz) const
{
    return evaluate(predictedTime, dtdt0, dtdx, dtdy, dtdz);
}

/// HDF5 support?
bool Static::compiledWithHDF5() noexcept
{
#ifdef WITH_HDF5
    return true;
#else
    return false;
#endif
}

/// Saves the model
void Static::save(const std::string &fileName) const
{
    if (!compiledWithHDF5())
    {
        throw std::runtime_error("Library not compiled with HDF5");
    }
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    if (!haveCorrection())
    {
        throw std::runtime_error("Correction not set or computed");
    }
#ifdef WITH_HDF5
    double correction = getCorrection();
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
             throw std::runtime_error("Could not create station group "
                                    + groupName);
        }
    }
    else
    {
        stationGroup = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
    }
    if (stationGroup == H5I_INVALID_HID)
    {
        H5Fclose(file);
        throw std::runtime_error("Could not open station group");
    }
    // Create a static correction group in the station group
    hid_t staticCorrectionGroup;
    if (!H5Lexists(stationGroup, STATIC_CORRECTION_GROUP_NAME, H5P_DEFAULT))
    {
        staticCorrectionGroup = H5Gcreate2(stationGroup,
                                           STATIC_CORRECTION_GROUP_NAME,
                                           H5P_DEFAULT, H5P_DEFAULT,
                                           H5P_DEFAULT);
        if (staticCorrectionGroup == H5I_INVALID_HID)
        {
            H5Gclose(stationGroup);
            H5Fclose(file);
            throw std::runtime_error(
                                   "Could not make static correction group "
                                   + std::string(STATIC_CORRECTION_GROUP_NAME));
        }
    }
    else
    {
        staticCorrectionGroup
            = H5Gopen(stationGroup, STATIC_CORRECTION_GROUP_NAME, H5P_DEFAULT);
    }
    if (staticCorrectionGroup == H5I_INVALID_HID)
    {
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error("Could not open static correction group");
    }
    // Remove
    if (H5Lexists(staticCorrectionGroup, STATIC_CORRECTION_NAME, H5P_DEFAULT))
    {
        status = H5Ldelete(staticCorrectionGroup, STATIC_CORRECTION_NAME,
                           H5P_DEFAULT);
        if (status < 0)
        {
            H5Gclose(staticCorrectionGroup);
            H5Gclose(stationGroup);
            H5Fclose(file);
            std::runtime_error("Failed to delete "
                              + std::string {STATIC_CORRECTION_NAME});
        }
    }
    // Write scalars
    const std::array<size_t, 1> onesDimension{1};
    auto scalarSpace = H5Screate_simple(1, onesDimension.data(), nullptr);
    auto dataSet = H5Dcreate2(staticCorrectionGroup, STATIC_CORRECTION_NAME,
                              H5T_NATIVE_DOUBLE, scalarSpace,
                              H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Dwrite(dataSet, H5T_NATIVE_DOUBLE, H5S_ALL, scalarSpace,
             H5P_DEFAULT, &correction);
    H5Dclose(dataSet);
    H5Sclose(scalarSpace);
    H5Gclose(staticCorrectionGroup);
    H5Gclose(stationGroup);
    H5Fclose(file);
#endif
}

/// Loads the model
void Static::load(const std::string &fileName)
{
    if (!compiledWithHDF5())
    {
        throw std::runtime_error("Library not compiled with HDF5");
    }
    if (!haveStationNameAndPhase())
    {
        throw std::runtime_error("Station name and phase not set");
    }
    if (!std::filesystem::exists(fileName))
    {
        throw std::invalid_argument(fileName + " does not exist");
    }
#ifdef WITH_HDF5
    // Open
    auto file = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    auto groupName = getNetwork() + "."
                   + getStation() + "."
                   + getPhase();
    herr_t status;
    if (!H5Lexists(file, groupName.c_str(), H5P_DEFAULT))
    {
        H5Fclose(file);
        throw std::runtime_error(groupName + " group doesn't exist");
    }
    auto stationGroup = H5Gopen(file, groupName.c_str(), H5P_DEFAULT);
    if (!H5Lexists(stationGroup, STATIC_CORRECTION_GROUP_NAME, H5P_DEFAULT))
    {
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error(std::string{STATIC_CORRECTION_GROUP_NAME}
                               + " correction group doesn't exist");
    }
    auto staticCorrectionGroup
        = H5Gopen(stationGroup, STATIC_CORRECTION_GROUP_NAME, H5P_DEFAULT);
    // Load the scalar
    if (!H5Lexists(staticCorrectionGroup, STATIC_CORRECTION_NAME, H5P_DEFAULT))
    {
        H5Gclose(staticCorrectionGroup);
        H5Gclose(stationGroup);
        H5Fclose(file);
        throw std::runtime_error(std::string{STATIC_CORRECTION_NAME}
                              + " dataset doesn't exist");
    }
    const std::array<size_t, 1> onesDimension{1};
    auto scalarSpace = H5Screate_simple(1, onesDimension.data(), nullptr);

    auto dataSet
        = H5Dopen2(staticCorrectionGroup, STATIC_CORRECTION_NAME, H5P_DEFAULT);
    auto dataSpace = H5Dget_space(dataSet);
    double staticCorrection;
    status = H5Dread(dataSet, H5T_NATIVE_DOUBLE, scalarSpace, dataSpace,
                     H5P_DEFAULT, &staticCorrection);
    H5Sclose(dataSpace);
    H5Dclose(dataSet);
    H5Gclose(staticCorrectionGroup);
    H5Gclose(stationGroup);
    H5Fclose(file);
    if (status != 0)
    {
        throw std::runtime_error("Failed t load static correction");
    }
    setCorrection(staticCorrection); 
#endif
}

void Static::clear() noexcept
{
    pImpl = std::make_unique<StaticImpl> ();
}
