#include <string>
#include <iostream>
#include <array>
#include <vector>
#include <cmath>
#include "uLocator/nloptOptions.hpp"
#include "uLocator/position/wgs84.hpp"
#include "defaultUtahQuarries.hpp"

using namespace ULocator;

class NLOptOptions::NLOptOptionsImpl
{
public:
    NLOptOptions::ObjectiveFunction mObjectiveFunction
    {
        NLOptOptions::ObjectiveFunction::LeastSquares
    };
    std::vector<Quarry> mQuarries;
    std::vector<double> mSearchDepths{-1000, 5000, 30000, 50000};
    std::array<double, 2> mLatitudeBoundaries;
    std::array<double, 2> mLongitudeBoundaries;
    std::array<double, 2> mDepthBoundaries;
    double mInitialEarthquakeDepth{6500}; // Utah avg is 6400, YNP avg is 6600
    double mDefaultElevation{1500};
    double mInitialQuarryBlastSearchDepth{-1500};
    double mInitialAbsoluteModelTolerance{1.e-1};
    double mAbsoluteModelTolerance{1.e-1};
    double mLatitudeRefinement{50000};
    double mLongitudeRefinement{50000};
    int mUTMZone{12};
    int mInitialMaximumFunctionEvaluations{2500};
    int mMaximumFunctionEvaluations{1000};
    bool mNorth{true};
    bool mHaveLatitudeBoundaries{false};
    bool mHaveLongitudeBoundaries{false};
    bool mHaveDepthBoundaries{false};
};

/// Constructor 
NLOptOptions::NLOptOptions() :
    pImpl(std::make_unique<NLOptOptionsImpl> ())
{
}

/// Constructor in a given region
NLOptOptions::NLOptOptions(const Region region) :
    pImpl(std::make_unique<NLOptOptionsImpl> ())
{
    if (region == Region::Utah)
    {
        setInitialEarthquakeSearchDepth(6400);
        setInitialQuarryBlastSearchDepth(2000); 
        setLatitudeBoundaries(std::array<double, 2> {  36.25,   43.0});
        setLongitudeBoundaries(std::array<double, 2> {-114.75, -108.25});
        setDepthBoundaries(std::array<double, 2> {-4500, 65000});
        pImpl->mQuarries = ::defaultUtahQuarries();
    }   
    else if (region == Region::Yellowstone)
    {
        setInitialEarthquakeSearchDepth(6600);
        setLatitudeBoundaries(std::array<double, 2> {  43.5,   45.5});
        setLongitudeBoundaries(std::array<double, 2> {-111.5, -110.0});
        setDepthBoundaries(std::array<double, 2> {-4500, 35000});
    }
}

/// Copy constructor
NLOptOptions::NLOptOptions(const NLOptOptions &options)
{
    *this = options;
}

/// Move assignment
NLOptOptions::NLOptOptions(NLOptOptions &&options) noexcept
{
    *this = std::move(options);
}

/// Copy assignent
NLOptOptions& NLOptOptions::operator=(const NLOptOptions &options)
{
    if (&options == this){return *this;}
    pImpl = std::make_unique<NLOptOptionsImpl> (*options.pImpl);
    return *this;
}

/// Move assignment
NLOptOptions& NLOptOptions::operator=(NLOptOptions &&options) noexcept
{
    if (&options == this){return *this;}
    pImpl = std::move(options.pImpl);
    return *this;
}

/// Destructor
NLOptOptions::~NLOptOptions() = default;

/// Objective function
void NLOptOptions::setObjectiveFunction(
    const ObjectiveFunction objectiveFunction) noexcept
{
    pImpl->mObjectiveFunction = objectiveFunction;
}

NLOptOptions::ObjectiveFunction 
NLOptOptions::getObjectiveFunction() const noexcept
{
    return pImpl->mObjectiveFunction;
}

/// UTM zone
void NLOptOptions::setUTMZone(const int zone, const bool north)
{   
    if (zone < 0 || zone > 60)
    {   
        if (zone !=-1)
        {
            throw std::invalid_argument(
                "UTM zone must be -1 or in range [1,60]");
        }
    }
    pImpl->mUTMZone = zone;
    pImpl->mNorth = north;
}

int NLOptOptions::getUTMZone() const noexcept
{
    return pImpl->mUTMZone;
}

bool NLOptOptions::isNorthernHemisphere() const noexcept
{
    return pImpl->mNorth;
}

/// Elevation
void NLOptOptions::setDefaultElevation(const double elevation)
{
    if (elevation < -10000 || elevation > 10000)
    {
        throw std::invalid_argument(
            "Elevation must be in range [-10000,10000]");
    }
    pImpl->mDefaultElevation = elevation;
}

double NLOptOptions::getDefaultElevation() const noexcept
{
    return pImpl->mDefaultElevation;
}

/// Initial quarry blast depth
void NLOptOptions::setInitialQuarryBlastSearchDepth(const double depth)
{
    if (depth < -15000 || depth > 800000)
    {   
        throw std::invalid_argument("Invalid initial depth");
    }   
    pImpl->mInitialQuarryBlastSearchDepth = depth;
}

double NLOptOptions::getInitialQuarryBlastSearchDepth() const noexcept
{
    return pImpl->mInitialQuarryBlastSearchDepth;
}

/// Boundaries
void NLOptOptions::setLongitudeBoundaries(
    const std::array<double, 2> &boundaries)
{
    if (boundaries[0] >= boundaries[1])
    {
        throw std::invalid_argument(
            "First longitude boundary must be less than second boundary");
    }
    pImpl->mLongitudeBoundaries = boundaries;
    pImpl->mHaveLongitudeBoundaries = true;
}

std::array<double, 2> NLOptOptions::getLongitudeBoundaries() const
{
    if (!haveLongitudeBoundaries())
    {
        throw std::runtime_error("Longitude boundaries not set");
    }
    return pImpl->mLongitudeBoundaries;
}

bool NLOptOptions::haveLongitudeBoundaries() const noexcept
{
    return pImpl->mHaveLongitudeBoundaries;
}

void NLOptOptions::setLatitudeBoundaries(
    const std::array<double, 2> &boundaries)
{
    if (boundaries[0] >= boundaries[1])
    {   
        throw std::invalid_argument(
            "First latitude boundary must be less than second boundary");
    }
    if (boundaries[0] <-90)
    {
        throw std::invalid_argument("First latitude must be >= -90");
    }
    if (boundaries[1] > 90)
    {
        throw std::invalid_argument("Second latitude be <= 90");
    }
    pImpl->mLatitudeBoundaries = boundaries;
    pImpl->mHaveLatitudeBoundaries = true;
}

std::array<double, 2> NLOptOptions::getLatitudeBoundaries() const
{
    if (!haveLatitudeBoundaries())
    {
        throw std::runtime_error("Latitude boundaries not set");
    }
    return pImpl->mLatitudeBoundaries;
}

bool NLOptOptions::haveLatitudeBoundaries() const noexcept
{
    return pImpl->mHaveLatitudeBoundaries;
}

void NLOptOptions::setDepthBoundaries(
    const std::array<double, 2> &boundaries)
{
    if (boundaries[0] >= boundaries[1])
    {
        throw std::invalid_argument(
            "First depth boundary must be less than second boundary");
    }
    if (boundaries[0] <-10000)
    {
        throw std::invalid_argument(
            "First depth must be greater than -10000 m");
    }
    if (boundaries[1] > 800000)
    {
        throw std::invalid_argument(
            "Second depth must be less than 800000 m");
    }
    pImpl->mDepthBoundaries = boundaries;
    pImpl->mHaveDepthBoundaries = true;
}

std::array<double, 2> NLOptOptions::getDepthBoundaries() const
{
    if (!haveDepthBoundaries())
    {
        throw std::runtime_error("Depth boundaries not set");
    }
    return pImpl->mDepthBoundaries;
}

bool NLOptOptions::haveDepthBoundaries() const noexcept
{
    return pImpl->mHaveDepthBoundaries;
}

/// The default event depth
void NLOptOptions::setInitialEarthquakeSearchDepth(const double depth)
{
    if (depth < -15000 || depth > 800000)
    {   
        throw std::invalid_argument("Invalid initial depth");
    }
    pImpl->mInitialEarthquakeDepth = depth;
}

double NLOptOptions::getInitialEarthquakeSearchDepth() const noexcept
{
    return pImpl->mInitialEarthquakeDepth;
}

/// Latitude refinement
double NLOptOptions::getLatitudeRefinement() const noexcept
{
    return pImpl->mLatitudeRefinement;
}

/// Longitude refinement
double NLOptOptions::getLongitudeRefinement() const noexcept
{
    return pImpl->mLongitudeRefinement;
}

/// X criteria
void NLOptOptions::setAbsoluteModelTolerance(const double tolerance)
{
    if (tolerance < 0)
    {   
        throw std::invalid_argument("Tolerance cannot be negative");
    }   
    if (tolerance < std::sqrt(std::numeric_limits<double>::epsilon()))
    {   
        std::cerr << "Warning: tolerance < sqrt(mach_eps) isn't helpful"
                  << std::endl;
    }   
    pImpl->mAbsoluteModelTolerance = tolerance;
}

double NLOptOptions::getAbsoluteModelTolerance() const noexcept
{
    return pImpl->mAbsoluteModelTolerance;
}

void NLOptOptions::setInitialAbsoluteModelTolerance(const double tolerance)
{
    if (tolerance < 0)
    {   
        throw std::invalid_argument("Tolerance cannot be negative");
    }
    if (tolerance < std::sqrt(std::numeric_limits<double>::epsilon()))
    {
        std::cerr << "Warning: tolerance < sqrt(mach_eps) isn't helpful"
                  << std::endl;
    }
    pImpl->mInitialAbsoluteModelTolerance = tolerance;
}

double NLOptOptions::getInitialAbsoluteModelTolerance() const noexcept
{
    return pImpl->mInitialAbsoluteModelTolerance;
}

std::vector<double> NLOptOptions::getRefinementSearchDepths() const noexcept
{
    return pImpl->mSearchDepths;
}

std::vector<Quarry> NLOptOptions::getQuarries() const
{
    return pImpl->mQuarries;
}

