#include <iostream>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#ifndef NDEBUG
#include <cassert>
#endif
//#include <ttimes/ak135.hpp>
#include <ttimes/phase.hpp>
#include <eikonalxx/ray/layerSolver.hpp>
#include <eikonalxx/ray/path2d.hpp>
#include <eikonalxx/ray/segment2d.hpp>
#include <eikonalxx/ray/point2d.hpp>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
#include "uLocator/uussTravelTimeCalculator.hpp"
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"
#include "uLocator/staticCorrection.hpp"
#include "uLocator/sourceSpecificStationCorrection.hpp"
#include "uLocator/firstArrivalRayTracer.hpp"
#include "h5io.hpp"

using namespace ULocator;

namespace
{
std::string convertString(const std::string &s) 
{
    auto temp = s;
    temp.erase(std::remove(temp.begin(), temp.end(), ' '), temp.end());
    std::transform(temp.begin(), temp.end(), temp.begin(), ::toupper);
    return temp;
}
double interpolate(const double xi, const double yi,
                   const double x0, const double x1,
                   const double y0, const double y1,
                   const std::array<double, 4> &f,
                   const double dxi, const double dyi)
{
    auto fy0 = dxi*( (x1 - xi)*f[0] + (xi - x0)*f[1] );
    auto fy1 = dxi*( (x1 - xi)*f[2] + (xi - x0)*f[3] );
    return dyi*( (y1 - yi)*fy0 + (yi - y0)*fy1 );
}
}

class UUSSTravelTimeCalculator::UUSSTravelTimeCalculatorImpl
{
public:
    UUSSTravelTimeCalculatorImpl(
        std::shared_ptr<UMPS::Logging::ILog> logger = nullptr) :
        mLogger(logger)
    {
        if (mLogger == nullptr)
        {
            mLogger = std::make_shared<UMPS::Logging::StandardOut> ();
        }
    }
/*
    void evaluateAK135(const double epicentralDistance, 
                       const double depthInKilometers,
                       double *travelTime, double *dtdz)
    {
        mAK135.compute(depthInKilometers, epicentralDistance);
        TTimes::Phase phase;
        if (mPhase == "P")
        {
            phase = mAK135.getFirstPPhase();
        }
        else if (mPhase == "S")
        {
            phase = mAK135.getFirstSPhase();
        }
        if (travelTime != nullptr){*travelTime = phase.getTime();}
        if (dtdz != nullptr)
        {
            // seconds/km -> seconds/m
            *dtdz = phase.getChangeInTravelTimeWithDepth()*1.e-3;
        }
    }
*/
    // Evaluate look-up table
    void evaluateLookUpTable(const double offset, const double depth,
                             double *travelTime,
                             double *dtdx, double *dtdz) const
    {
        /// Find the right cell
        double xOffset = std::min(std::max(0., offset - mOriginInX), mWidthInX);
        double zOffset = std::min(std::max(0., depth  - mOriginInZ), mWidthInZ);
        auto iCellX
            = std::min(mNumberOfCellsInX - 1,
                       static_cast<int> (xOffset/mGridSpacingInX));
        auto iCellZ
            = std::min(mNumberOfCellsInZ - 1,
                       static_cast<int> (zOffset/mGridSpacingInZ));
#ifndef NDEBUG
        assert(iCellX >= 0 && iCellX < mNumberOfCellsInX); 
        assert(iCellZ >= 0 && iCellZ < mNumberOfCellsInZ);
#endif
        if (iCellX != mCellX0 || iCellZ != mCellZ0)
        {
            auto i0 = 3*(iCellZ*mNumberOfGridPointsInX + iCellX);
            auto i1 = i0 + 3;
            auto i2 = i0 + 3*mNumberOfGridPointsInX;
            auto i3 = i2 + 3;
            mTravelTimes0[0]   = mTravelTimeAndGradientField[i0]; 
            mTravelTimesDx0[0] = mTravelTimeAndGradientField[i0 + 1];
            mTravelTimesDz0[0] = mTravelTimeAndGradientField[i0 + 2];

            mTravelTimes0[1]   = mTravelTimeAndGradientField[i1];
            mTravelTimesDx0[1] = mTravelTimeAndGradientField[i1 + 1];
            mTravelTimesDz0[1] = mTravelTimeAndGradientField[i1 + 2];

            mTravelTimes0[2]   = mTravelTimeAndGradientField[i2];
            mTravelTimesDx0[2] = mTravelTimeAndGradientField[i2 + 1];
            mTravelTimesDz0[2] = mTravelTimeAndGradientField[i2 + 2];

            mTravelTimes0[3]   = mTravelTimeAndGradientField[i3];
            mTravelTimesDx0[3] = mTravelTimeAndGradientField[i3 + 1];
            mTravelTimesDz0[3] = mTravelTimeAndGradientField[i3 + 2];

            mCellX0 = iCellX; 
            mCellZ0 = iCellZ;
        }
        // Interpolate
        double xi = xOffset - iCellX*mGridSpacingInX;
        double zi = zOffset - iCellZ*mGridSpacingInZ;
        *travelTime = ::interpolate(xi, zi,
                                    0, mGridSpacingInX,
                                    0, mGridSpacingInZ,
                                    mTravelTimes0,
                                    mInverseGridSpacingInX,
                                    mInverseGridSpacingInZ);
        if (dtdx != nullptr)
        {
            *dtdx = ::interpolate(xi, zi, 
                                  0, mGridSpacingInX,
                                  0, mGridSpacingInZ,
                                  mTravelTimesDx0,
                                  mInverseGridSpacingInX,
                                  mInverseGridSpacingInZ);
        }
        if (dtdz != nullptr)
        {
            *dtdz = ::interpolate(xi, zi,
                                  0, mGridSpacingInX,
                                  0, mGridSpacingInZ,
                                  mTravelTimesDz0,
                                  mInverseGridSpacingInX,
                                  mInverseGridSpacingInZ);
            *dtdz =-*dtdz;
        }
    }
    [[nodiscard]] double applyStaticCorrection(const double time) const noexcept
    {
        if (!mStaticCorrection.haveCorrection()){return time;}
        return mStaticCorrection.evaluate(time);
    }
    [[nodiscard]] double applySourceSpecificStationCorrection(
        const double latitude, const double longitude, const double depth,
        const double time) const noexcept
    {
        if (!mSourceSpecificStationCorrection.isInitialized())
        {
            return time;
        }
        if (!mSourceSpecificStationCorrection.haveModel())
        {
            return time;
        }
        double correction = 0;
        try
        {
            correction
                = mSourceSpecificStationCorrection.evaluate(latitude,
                                                            longitude,
                                                            depth,
                                                            ULocator::SourceSpecificStationCorrection::EvaluationMethod::InverseDistanceWeighted);
        }
        catch (const std::exception &e) 
        {
            mLogger->warn("Failed to evaluate SSSC: " + std::string {e.what()}
                        + " - setting correction to 0");
        }
        return correction + time; 
    }
    void loadStaticCorrection(const std::string &fileName)
    {
#ifndef NDEBUG
        assert(mStaticCorrection.isInitialized());
#endif
        try
        {
            mStaticCorrection.load(fileName);
            mLogger->debug("Loaded static correction: "
                         + std::to_string(mStaticCorrection.getCorrection())
                         + " from " + fileName);
        }
        catch (const std::exception &e)
        {
            mLogger->warn("Failed to load static correction:  "
                         + std::string {e.what()}
                         + " -  Setting correction to 0");
            mStaticCorrection.setCorrection(0);
        }
    }
    void loadSourceSpecificStationCorrection(const std::string &fileName)
    {
#ifndef NDEBUG
        assert(mSourceSpecificStationCorrection.isInitialized());
#endif
        try
        {
            mSourceSpecificStationCorrection.load(fileName);
            mLogger->debug(
                "Loaded source-specific station correction from "
              + fileName);
        }
        catch (const std::exception &e)
        {
            mLogger->warn(
                "Failed to load source-specific station correction from "
              + fileName + " failed with: " + std::string(e.what()));
        }
    }
    void load(const std::string &fileName,
              const std::string &network, const std::string &station,
              const std::string &phase,
              const bool isYellowstone = false)
    {
        if (!std::filesystem::exists(fileName))
        {
            throw std::runtime_error("H5 file: "
                                   + fileName + " does not exist");
        }
        auto fileID = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        auto groupName = ::convertString(network + "." + station);
        if (!H5Lexists(fileID, groupName.c_str(), H5P_DEFAULT))
        {
            H5Fclose(fileID);
            throw std::runtime_error("Group " + groupName + " does not exist");
        }
        auto groupID = H5Gopen2(fileID, groupName.c_str(), H5P_DEFAULT);
        std::vector<float> travelTimes;
        std::vector<float> dTravelTimesDx;
        std::vector<float> dTravelTimesDz;
        std::vector<double> scalar;
        std::vector<size_t> dimensions;
        // Things to always read
        try
        {
            ::readDataset(groupID, "x0", &scalar);
            mOriginInX = scalar.at(0);
            ::readDataset(groupID, "z0", &scalar);
            mOriginInZ = scalar.at(0); 
            ::readDataset(groupID, "dx", &scalar);
            mGridSpacingInX = scalar.at(0); 
            ::readDataset(groupID, "dz", &scalar);
            mGridSpacingInZ = scalar.at(0);
            ::readDataset(groupID, "station_elevation", &scalar);
            mStationDepth =-scalar.at(0);
            mLogger->debug("Station depth: "
                         + std::to_string(mStationDepth));
        }
        catch (const std::exception &e) 
        {
            H5Gclose(groupID);
            H5Fclose(fileID);
            throw std::invalid_argument("Failed to read scalar values: "
                                      + std::string(e.what()));
        }
        // P information
        if (phase == "P")
        {
            if (!H5Lexists(groupID, "p_travel_times", H5P_DEFAULT))
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::runtime_error("p_travel_times does not exist");
            }
            if (!H5Lexists(groupID, "p_dtdx", H5P_DEFAULT))
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::runtime_error("p_dtdx does not exist");
            }
            if (!H5Lexists(groupID, "p_dtdz", H5P_DEFAULT))
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::runtime_error("p_dtdz does not exist");
            }
            /*
            if (H5Lexists(groupID, "p_static_correction", H5P_DEFAULT))
            {
                std::vector<double> correction;
                ::readDataset(groupID, "p_static_correction", &correction);
                //mStaticCorrection = correction[0];
            }
            */
            try
            {
                ::readDataset(groupID, "p_travel_times", &travelTimes, &dimensions);
                ::readDataset(groupID, "p_dtdx", &dTravelTimesDx);
                ::readDataset(groupID, "p_dtdz", &dTravelTimesDz);
            }
            catch (const std::exception &e)
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::invalid_argument("Failed to read P field: "
                                          + std::string(e.what()));
            }
            if (dimensions.size() != 2)
            {
                H5Gclose(groupID);
                H5Fclose(fileID);                 
                throw std::invalid_argument("Field should have two dimensions");
            }
            mNumberOfGridPointsInX = static_cast<int> (dimensions[1]);
            mNumberOfGridPointsInZ = static_cast<int> (dimensions[0]);
        }
        else if (phase == "S")
        {
            if (!H5Lexists(groupID, "s_travel_times", H5P_DEFAULT))
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::runtime_error("p_travel_times does not exist");
            }
            if (!H5Lexists(groupID, "s_dtdx", H5P_DEFAULT))
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::runtime_error("s_dtdx does not exist");
            }
            if (!H5Lexists(groupID, "s_dtdz", H5P_DEFAULT))
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::runtime_error("s_dtdz does not exist");
            }
            /*
            if (H5Lexists(groupID, "s_static_correction", H5P_DEFAULT))
            {
                std::vector<double> correction;
                ::readDataset(groupID, "s_static_correction", &correction);
                //mStaticCorrection = correction[0];
            }
            */
            try
            {
                ::readDataset(groupID, "s_travel_times", &travelTimes, &dimensions);
                ::readDataset(groupID, "s_dtdx", &dTravelTimesDx);
                ::readDataset(groupID, "s_dtdz", &dTravelTimesDz);
            }
            catch (const std::exception &e) 
            {
                H5Gclose(groupID);
                H5Fclose(fileID);
                throw std::invalid_argument("Failed to read S field: "
                                          + std::string(e.what()));
            }
            if (dimensions.size() != 2)
            {
                H5Gclose(groupID);
                H5Fclose(fileID);    
                throw std::invalid_argument("Field should have two dimensions");
            }
            mNumberOfGridPointsInX = static_cast<int> (dimensions[1]);
            mNumberOfGridPointsInZ = static_cast<int> (dimensions[0]);
        }
#ifndef NDEBUG
        else
        {
           assert(false);
        }
#endif
        H5Gclose(groupID);
        H5Fclose(fileID);
        // Unpack
        mNumberOfCellsInX = mNumberOfGridPointsInX - 1;
        mNumberOfCellsInZ = mNumberOfGridPointsInZ - 1;
        auto nGrid = mNumberOfGridPointsInX*mNumberOfGridPointsInZ;
        mTravelTimeAndGradientField.resize(3*nGrid);
        for (int k = 0; k < mNumberOfGridPointsInZ; ++k)
        {
            for (int i = 0; i < mNumberOfGridPointsInX; ++i)
            {
                auto indx = k*mNumberOfGridPointsInX + i;
                mTravelTimeAndGradientField[3*indx]     = travelTimes.at(indx);
                mTravelTimeAndGradientField[3*indx + 1] = dTravelTimesDx[indx];
                mTravelTimeAndGradientField[3*indx + 2] = dTravelTimesDz[indx];
            }
        }
        mWidthInX = mNumberOfCellsInX*mGridSpacingInX;
std::cerr << "fix this" << std::endl;
mWidthInX = 0;
        mWidthInZ = mNumberOfCellsInZ*mGridSpacingInZ;
        mInverseGridSpacingInX = 1./mGridSpacingInX;
        mInverseGridSpacingInZ = 1./mGridSpacingInZ;
        mPhase = phase;
        mName = network + "." + station + "." + mPhase;
        initializeLayerSolver(isYellowstone);
    }
    void initializeLayerSolver(const bool isYellowstone)
    {
        const std::vector<double> utahInterfaces{-4500,  40,  15600, 26500, 40500};
        //const std::vector<double> utahPVelocities{3400, 5900,  6400,  7500,  7900};
        //const std::vector<double> utahPVelocities{3400, 5950, 6450,  7550,  7900};
        const std::vector<double> utahPVelocities{3685, 5972, 6553, 7106, 7672};
        //const std::vector<double> utahSVelocities{1950, 3390, 3680,   4310,  4540};
        //const std::vector<double> utahSVelocities{2000, 3425, 3700,   4400,  4550};
        const std::vector<double> utahSVelocities{2169, 3443, 3653, 4486, 5022};
        const std::vector<double> ynpInterfaces{-4500, -1000,  2000,  5000,  8000,
                                                12000, 16000, 21000, 50000};
        //const std::vector<double> ynpPVelocities{2720, 2790, 5210, 5560, 5770,
        //                                         6070, 6330, 6630, 8000};
        const std::vector<double> ynpPVelocities{2512, 3398, 4689, 5456, 5674,
                                                 6250, 6398, 6574, 8200};
        //const std::vector<double> ynpSVelocities{1950, 2000, 3400, 3420, 3490,
        //                                         3680, 3780, 4000, 4850};
        const std::vector<double> ynpSVelocities{1725, 2343, 3064, 3425, 3569,
                                                 3690, 3705, 3975, 4950};
        /*
        const std::vector<double> ynpPVelocities{2750, 2810, 5210, 5560, 5770,
                                                 6070, 6330, 6630, 8000};
        const std::vector<double> ynpSVelocities{1990, 2190, 3150, 3300, 3400,
                                                 3550, 3625, 3750, 4700};
        */
        std::vector<double> velocities;
        std::vector<double> interfaces;
        if (isYellowstone)
        {
            interfaces = ynpInterfaces;
            if (mPhase == "S")
            {
                velocities = ynpSVelocities;
                mLogger->info("Initializing YNP S velocities for " + mName
                            + " at depth " + std::to_string(mStationDepth));
            }
            else
            {
                velocities = ynpPVelocities;
                mLogger->info("Initializing YNP P velocities for " + mName
                            + " at depth " + std::to_string(mStationDepth));
            }
        }
        else
        {
            interfaces = utahInterfaces;
            if (mPhase == "P")
            {
                velocities = utahPVelocities;
                mLogger->info("Initializing Utah P velocities for " + mName
                            + " at depth " + std::to_string(mStationDepth));
            }
            else
            {
                velocities = utahSVelocities;
                mLogger->info("Initializing Utah S velocities for " + mName
                            + " at depth " + std::to_string(mStationDepth));
            }
        }
#ifndef NDEBUG
        assert(!velocities.empty());
        assert(!interfaces.empty());
#endif
        mMinimumDepth = interfaces.at(0);
        mLayerSolver.setVelocityModel(interfaces, velocities);
    }
    double evaluateLayerSolver(const double sourceDepth, const double offset,
                               double *dtdx, double *dtdz) const
    {
        double travelTime{0};
        if (dtdx != nullptr){*dtdx = 0;}
        if (dtdz != nullptr){*dtdz = 0;}
        try
        {
            auto sourceDepthToUse = sourceDepth;
            if (sourceDepth < mMinimumDepth)
            {
                sourceDepthToUse = mMinimumDepth;
                mLogger->warn("Overriding source depth "
                            + std::to_string(sourceDepth)
                            + " to " + std::to_string(mMinimumDepth));
            }
            //std::cout << offset << "," << mStationDepth << "," << sourceDepthToUse << std::endl;
            mLayerSolver.setStationOffsetAndDepth(offset, mStationDepth);
            mLayerSolver.setSourceDepth(sourceDepthToUse);
            mLayerSolver.solve();
        }
        catch (const std::exception &e)
        {
            mLogger->error("Failed to compute ray paths for "
                         + mStation.getNetwork() + "." 
                         + mStation.getName() + "."
                         + mPhase
                         + ".  Failed with: "
                         + std::string {e.what()}
                         +  ".  Source depth is: " 
                         + std::to_string(sourceDepth));
            return travelTime;
        }
        const auto &rayPaths = mLayerSolver.getRayPaths(); 
        if (!rayPaths.empty())
        {
            travelTime = rayPaths[0].getTravelTime();
            // Decompose segments in first leg into slowness in x and z
            // as these are the sensitivities at the source.
            if (dtdx != nullptr || dtdz != nullptr)
            {
                auto takeOffAngle = rayPaths[0].getTakeOffAngle()*(M_PI/180);
                auto slowness = rayPaths[0].at(0).getSlowness();
                if (dtdx != nullptr){*dtdx = slowness*std::sin(takeOffAngle);}
                if (dtdz != nullptr){*dtdz = slowness*std::cos(takeOffAngle);}
            }
        } 
        return travelTime;
    }
    std::shared_ptr<UMPS::Logging::ILog> mLogger{nullptr};
    mutable EikonalXX::Ray::LayerSolver mLayerSolver;
    StaticCorrection mStaticCorrection;
    SourceSpecificStationCorrection mSourceSpecificStationCorrection;
    //TTimes::AK135 mAK135;
    Station mStation;
    Position::WGS84 mStationLocation;
    std::string mPhase;
    std::string mName;
    std::vector<float> mTravelTimeAndGradientField;
    //double mStaticCorrection{0};
    double mOriginInX{0};
    double mOriginInZ{-5000};
    double mWidthInX{0};
    double mWidthInZ{0};
    double mStationDepth{0};
    double mGridSpacingInX{125};
    double mGridSpacingInZ{125};
    double mInverseGridSpacingInX{1.0/mGridSpacingInX};
    double mInverseGridSpacingInZ{1.0/mGridSpacingInZ};
    double mMinimumDepth{-4500};
    mutable std::array<double, 4> mTravelTimes0;
    mutable std::array<double, 4> mTravelTimesDx0;
    mutable std::array<double, 4> mTravelTimesDz0;
    int mNumberOfCellsInX{0};
    int mNumberOfCellsInZ{0};
    int mNumberOfGridPointsInX{0};
    int mNumberOfGridPointsInZ{0};
    mutable int mCellX0{-1};
    mutable int mCellZ0{-1};
    bool mInitialized{false};
};

UUSSTravelTimeCalculator::UUSSTravelTimeCalculator() :
    pImpl(std::make_unique<UUSSTravelTimeCalculatorImpl> ())
{
}

UUSSTravelTimeCalculator::UUSSTravelTimeCalculator(
    std::shared_ptr<UMPS::Logging::ILog> &logger) :
    pImpl(std::make_unique<UUSSTravelTimeCalculatorImpl> (logger))
{
}

UUSSTravelTimeCalculator::~UUSSTravelTimeCalculator() = default;

void UUSSTravelTimeCalculator::clear() noexcept
{
    pImpl = std::make_unique<UUSSTravelTimeCalculatorImpl> ();
}

double UUSSTravelTimeCalculator::evaluate(
    const Position::WGS84 &sourcePosition, const double depth,
    const bool applyCorrection) const
{
    // Compute the distance and source-to-receiver azimuth
    double travelTime;
    double distance, epicentralDistance;
    ULocator::Position::computeDistanceAzimuth(sourcePosition,
                                               pImpl->mStationLocation,
                                               &epicentralDistance,
                                               &distance,
                                               nullptr,
                                               nullptr);
    if (distance > pImpl->mWidthInX)
    {
        travelTime = pImpl->evaluateLayerSolver(depth, distance, nullptr, nullptr);
//std::cout << " layer " << travelTime;
        if (applyCorrection)
        {
            auto sourceLatitude = sourcePosition.getLatitude();
            auto sourceLongitude = sourcePosition.getLongitude();
            travelTime = pImpl->applyStaticCorrection(travelTime);// + pImpl->evaluateStaticCorrection();
            travelTime
                = pImpl->applySourceSpecificStationCorrection(
                      sourceLatitude, sourceLongitude, depth, travelTime);
        }
        return travelTime;
/*
        auto depthInKilometers = std::max(0.0, depth*1.e-3);
        if (depthInKilometers > TTimes::AK135::getMaximumDepth())
        {
            throw std::invalid_argument("Input depth: "
                        + std::to_string(depthInKilometers)
                        + " exceeds maximum depth of " 
                        + std::to_string(TTimes::AK135::getMaximumDepth()));
        }
double ak135TravelTime;
        pImpl->evaluateAK135(epicentralDistance, depthInKilometers,
                             &ak135TravelTime, nullptr);
        std::cout << " ak135 " << ak135TravelTime << std::endl << std::endl;
        if (applyCorrection)
        {
            travelTime = travelTime + pImpl->evaluateStaticCorrection();
        }
        return travelTime;
*/
    }
    auto zMin = pImpl->mOriginInZ;
    auto zMax = zMin + pImpl->mWidthInZ;
    if (depth < zMin || depth > zMax)
    {
        throw std::invalid_argument("Depth must be in range ["
                                  + std::to_string(zMin) + ","
                                  + std::to_string(zMax) + "]");
    }
    // Look up the travel time
    pImpl->evaluateLookUpTable(distance, depth,
                               &travelTime,
                               nullptr, nullptr);
    //std::cout << " eikonal: " << travelTime << std::endl; 
    if (applyCorrection)
    {
        auto sourceLatitude = sourcePosition.getLatitude();
        auto sourceLongitude = sourcePosition.getLongitude();
        travelTime = pImpl->applyStaticCorrection(travelTime); //travelTime + pImpl->evaluateStaticCorrection();
        travelTime
            = pImpl->applySourceSpecificStationCorrection(
                 sourceLatitude, sourceLongitude, depth, travelTime);
    }
    return travelTime;
}

double UUSSTravelTimeCalculator::evaluate(
    const Position::WGS84 &sourcePosition, const double depth,
    double *dtdx, double *dtdy, double *dtdz,
    const bool applyCorrection) const
{
    //std::cout << sourcePosition.getLongitude() << "," << sourcePosition.getLatitude() << std::endl;
    //std::cout << pImpl->mStationLocation.getLongitude() << "," << pImpl->mStationLocation.getLatitude() << std::endl;
    constexpr double h{2.5}; // Finite difference in meters
    constexpr bool finiteDifference{false};
    constexpr double degreesToRadians{M_PI/180.};
    // Compute the distance and source-to-receiver azimuth
    double travelTime;
    double epicentralDistance, distance, azimuth;
    ULocator::Position::computeDistanceAzimuth(sourcePosition,
                                               pImpl->mStationLocation,
                                               &epicentralDistance,
                                               &distance,
                                               &azimuth,
                                               nullptr);
    double dtdr;
/*
    travelTime = pImpl->evaluateLayerSolver(depth, distance, &dtdr, dtdz);
//std::cout << depth << "," << distance << std::endl;
//std::cout << "layer dtdr, dtdz, ttime: " << dtdr << "," << *dtdz << "," << travelTime << " " << pImpl->mPhase << std::endl;
//std::cout << (pImpl->evaluateLayerSolver(depth + 1, distance, nullptr, nullptr) - travelTime)/1 << std::endl;
    if (applyCorrection)
    {
        auto sourceLatitude = sourcePosition.getLatitude();
        auto sourceLongitude = sourcePosition.getLongitude();
        travelTime = pImpl->applyStaticCorrection(travelTime);// + pImpl->evaluateStaticCorrection();
        travelTime
            = pImpl->applySourceSpecificStationCorrection(
                 sourceLatitude, sourceLongitude, depth, travelTime);
    }
*/
//travelTime = pImpl->evaluateLayerSolver(depth, distance, &dtdr, dtdz);
//std::cout << *dtdz << std::endl;
    if (distance > pImpl->mWidthInX)
    {
        double dtdr;
        travelTime = pImpl->evaluateLayerSolver(depth, distance, &dtdr, dtdz);
        if (applyCorrection)
        {
            auto sourceLatitude = sourcePosition.getLatitude();
            auto sourceLongitude = sourcePosition.getLongitude();
            travelTime = pImpl->applyStaticCorrection(travelTime);
            travelTime
                = pImpl->applySourceSpecificStationCorrection(
                     sourceLatitude, sourceLongitude, depth, travelTime);
        }
        *dtdx = dtdr*std::sin(azimuth*degreesToRadians);
        *dtdy = dtdr*std::cos(azimuth*degreesToRadians);
        /*
        auto depthInKilometers = std::max(0.0, depth*1.e-3);
        if (depthInKilometers > TTimes::AK135::getMaximumDepth())
        {
            throw std::invalid_argument("Input depth: "
                        + std::to_string(depthInKilometers)
                        + " (km) exceeds maximum depth of " 
                        + std::to_string(TTimes::AK135::getMaximumDepth()));
        }
        pImpl->evaluateAK135(epicentralDistance, depthInKilometers,
                             &ak135TravelTime, nullptr);
        */
        //*dtdx = 0;
        //*dtdy = 0;
        return travelTime;
    }   
    auto zMin = pImpl->mOriginInZ;
    auto zMax = zMin + pImpl->mWidthInZ;
    if (depth < zMin || depth > zMax)
    {   
        throw std::invalid_argument("Depth must be in range ["
                                  + std::to_string(zMin) + "," 
                                  + std::to_string(zMax) + "]");
    }
    // Look up table
    pImpl->evaluateLookUpTable(distance, depth,
                               &travelTime,
                               &dtdr, dtdz);
//std::cout << *dtdz << std::endl;
//getchar();
//std::cout << "eikonal dtdr, dtdz, ttime " << dtdr << "," << *dtdz << "," << travelTime << std::endl;
    if (applyCorrection)
    {
        auto sourceLatitude = sourcePosition.getLatitude();
        auto sourceLongitude = sourcePosition.getLongitude();
        travelTime = pImpl->applyStaticCorrection(travelTime);// + pImpl->evaluateStaticCorrection();
        travelTime
            = pImpl->applySourceSpecificStationCorrection(
                 sourceLatitude, sourceLongitude, depth, travelTime);
    }
    // Use a finite difference, otherwise, use math
    if (finiteDifference)
    {
        auto utmZone = sourcePosition.getUTMZone();
        auto isNorth = sourcePosition.isNorth();
        auto utmX = sourcePosition.getEasting();
        auto utmY = sourcePosition.getNorthing();
        Position::WGS84 dSourceX(utmZone, isNorth, utmX + h, utmY);
        Position::WGS84 dSourceY(utmZone, isNorth, utmX, utmY + h);
        auto dTimeX = evaluate(dSourceX, depth, applyCorrection);
        auto dTimeY = evaluate(dSourceY, depth, applyCorrection);
        *dtdx = (dTimeX - travelTime)/h;
        *dtdy = (dTimeY - travelTime)/h;
        return travelTime;
    }
    // Nominally, we have a look up table of the form T(r,z)
    // {dT/dr}   = [   cos(theta)   sin(theta) 0 ] {dT/dx}
    // {dT/dphi} = [-r sin(theta) r cos(theta) 0 ] {dT/dy}
    // {dT/dz}   = [     0              0      1 ] {dT/dz}
    // The inverse relations are 
    // {dT/dx}   = [   cos(theta)  -sin(theta)/r 0 ] {dT/dr}
    // {dT/dz}   = [   sin(theta)   cos(theta)/r 0 ] {dT/dtheta}
    // {dT/dz}   = ]      0             0        1 ] {dT/dz}
    // Note, the models are 1D so dT/dtheta = 0.  Thus:
    //   dT/dx = cos(theta) dT/dr
    //   dT/dy = sin(theta) dT/dr
    //   dT/dz = dT/dz
    // However, theta increases + from east whereas our azimuth
    // is measured positive east of north.  So we have to flip
    // theta to 90 - theta.  Then
    //   cos(90 - azimuth) = cos(90) cos(theta) + sin(90) sin(theta)
    //                     = sin(theta)
    //   sin(90 - azimuth) = sin(90) cos(theta) - cos(90) sin(theta)
    //                     = cos(theta)
    //std::cout << azimuth << std::endl;
    if (distance > 0)
    {
        // Note, these are reciprocity grids so we are pointing
        // the wrong way.
        *dtdx =-dtdr*std::sin(azimuth*degreesToRadians);
        *dtdy =-dtdr*std::cos(azimuth*degreesToRadians);
    }
    else
    {
        //auto theta = (90 - azimuth)*degreesToRadians;
        //*dtdx = std::cos(theta)*dtdr
        //      - 1./distance*std::sin(theta)*  (90 - azimuth)
        *dtdx = 0;
        *dtdy = 0;
    }
    return travelTime;
}

void UUSSTravelTimeCalculator::load(
    const std::string &fileName,
    const Station &station,
    const std::string &phase,
    const std::string &region)
{
    if (!std::filesystem::exists(fileName))
    {
        throw std::invalid_argument(fileName + " does not exist");
    }
    if (!station.haveNetwork())
    {
        throw std::invalid_argument("Network code not set");
    }
    if (!station.haveName())
    {
        throw std::invalid_argument("Station name not set");
    }
    if (!station.haveGeographicPosition())
    {
        throw std::invalid_argument("Station's geographic position not set");
    }
    if (phase != "P" && phase != "S")
    {
        throw std::invalid_argument("Phase must be P or S");
    }
    bool isYellowstone = false;
    if (region == "ynp"){isYellowstone = true;}
    clear();
    pImpl->load(fileName, station.getNetwork(), station.getName(), phase, isYellowstone);
    pImpl->mStation = station;
    pImpl->mStationLocation = pImpl->mStation.getGeographicPosition();
    pImpl->mStaticCorrection.initialize(station.getNetwork(),
                                        station.getName(),
                                        phase);
    pImpl->mSourceSpecificStationCorrection.initialize(station.getNetwork(),
                                                       station.getName(),
                                                       phase);
    pImpl->mInitialized = true;
}

void UUSSTravelTimeCalculator::loadStaticCorrection(const std::string &fileName)
{
    if (!pImpl->mInitialized)//isInitialized())
    {
        throw std::runtime_error("Calculator not initialized");
    }
    if (!std::filesystem::exists(fileName))
    {
        throw std::runtime_error(fileName + " does not exist");
    }
    pImpl->loadStaticCorrection(fileName);
}

void UUSSTravelTimeCalculator::loadSourceSpecificStationCorrections(
    const std::string &fileName)
{
    if (!pImpl->mInitialized)//isInitialized())
    {
        throw std::runtime_error("Calculator not initialized");
    }
    if (!std::filesystem::exists(fileName))
    {
        throw std::runtime_error(fileName + " does not exist");
    }
    pImpl->loadSourceSpecificStationCorrection(fileName);
}
