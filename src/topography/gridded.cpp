#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <filesystem>
#include "uLocator/topography/gridded.hpp"
#include "position/lonTo180.hpp"
#include "h5io.hpp"

using namespace ULocator::Topography;

namespace 
{
template<typename T>
void checkLatitudes(const int n, const T *x)
{
    if (x == nullptr){throw std::runtime_error("latitudes is NULL");}
    if (n < 2){throw std::runtime_error("At least two latitudes required");}
    for (int i = 0; i < n; ++i)
    {
        if (x[i] < -90 || x[i] > 90)
        {
            throw std::runtime_error("Latitudes must be in range [-90,90]");
        }
    }
    for (int i = 1; i < n; ++i)
    {
        if (x[i - 1] >= x[i])
        {
            throw std::runtime_error("Latitudes must be increasing");
        }
    }
}
template<typename T>
void checkLongitudes(const int n, const T *x) 
{
    if (x == nullptr){throw std::runtime_error("longitudes is NULL");}
    if (n < 2){throw std::runtime_error("At least two longitudes required");}
    for (int i = 1; i < n; ++i)
    {
        if (x[i - 1] >= x[i])
        {
            throw std::runtime_error("Longitudes must be increasing");
        }
    }
}
template<typename T>
void fixLongitudes(std::vector<T> &longitudes)
{
    for (int i = 0; i < static_cast<int> (longitudes.size()); ++i)
    {
        longitudes[i] = ::lonTo180(longitudes[i]);
    } 
}
}

class Gridded::GriddedImpl
{
public:
    std::vector<double> mElevations;
    std::vector<double> mLatitudes;
    std::vector<float> mLongitudes;
    double mZ00{0};
    double mZ01{0};
    double mZ10{0};
    double mZ11{0};
    int mLon0{-1};
    int mLat0{-1};
    bool mHaveTopography{false}; 
};

/// Constructor 
Gridded::Gridded() :
    pImpl(std::make_unique<GriddedImpl> ())
{
}

/// Copy constructor
Gridded::Gridded(const Gridded &topography)
{
    *this = topography;
}

/// Move constructor
Gridded::Gridded(Gridded &&topography) noexcept
{
    *this = std::move(topography);
}

/// Destructor
Gridded::~Gridded() = default;

/// Copy assignment
Gridded& Gridded::operator=(const Gridded &topography)
{
    if (&topography == this){return *this;}
    pImpl = std::make_unique<GriddedImpl> (*topography.pImpl);
    return *this;
}

/// Move assignment
Gridded& Gridded::operator=(Gridded &&topography) noexcept
{
    if (&topography == this){return *this;}
    pImpl = std::move(topography.pImpl);
    return *this;
}

/// Reset class and release memory
void Gridded::clear() noexcept
{
    pImpl = std::make_unique<GriddedImpl> ();
}

/// Load a model from HDF5
void Gridded::load(const std::string &fileName)
{
#ifdef WITH_HDF5
    if (!std::filesystem::exists(fileName))
    {
        throw std::runtime_error(fileName + " does not exist");
    }
    // Open
    auto mFile = H5Fopen(fileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    // Verify existence
    if (!H5Lexists(mFile, "Latitude", H5P_DEFAULT))
    {
        H5Fclose(mFile);
        throw std::runtime_error("Latitude dataset does not exist");
    }
    if (!H5Lexists(mFile, "Longitude", H5P_DEFAULT))
    {
        H5Fclose(mFile);
        throw std::runtime_error("Longitude dataset does not exist");
        H5Fclose(mFile);
    }
    if (!H5Lexists(mFile, "Elevation", H5P_DEFAULT))
    {
        H5Fclose(mFile);
        throw std::runtime_error("Elevation dataset does not exist");
    }
    try
    {
        ::readDataset(mFile, "Latitude",  &pImpl->mLatitudes);
        ::checkLatitudes(pImpl->mLatitudes.size(), pImpl->mLatitudes.data());
        ::readDataset(mFile, "Longitude", &pImpl->mLongitudes);
        ::checkLongitudes(pImpl->mLongitudes.size(), pImpl->mLongitudes.data());
        ::readDataset(mFile, "Elevation", &pImpl->mElevations);
        ::fixLongitudes(pImpl->mLongitudes);
    }
    catch (const std::exception &e)
    {
        pImpl->mLatitudes.clear();
        pImpl->mLongitudes.clear();
        pImpl->mElevations.clear();
        H5Fclose(mFile);
        throw e;
    }
    // Close it up
    H5Fclose(mFile); 
#else
    throw std::runtime_error("Library not compiled with HDF5");
#endif
    pImpl->mHaveTopography = true;
}

/// Have topography?
bool Gridded::haveTopography() const noexcept
{
    return pImpl->mHaveTopography;
}

/// Set the topography
template<typename U>
void Gridded::set(const int nLatitudes, const U *latitudes,
                  const int nLongitudes, const U *longitudes,
                  const int nGrid, const U *elevation)
{
    if (nLatitudes*nLongitudes != nGrid)
    {
        throw std::invalid_argument("nGrid != nLatitudes x nLongitudes");
    }
    if (elevation == nullptr)
    {
        throw std::invalid_argument("Elevation is NULL");
    }
    ::checkLatitudes(nLatitudes, latitudes);
    ::checkLongitudes(nLongitudes, longitudes);
    // Copy
    pImpl->mLatitudes.resize(nLatitudes);
    std::copy(latitudes, latitudes + nLatitudes, pImpl->mLatitudes.data());
    pImpl->mLongitudes.resize(nLongitudes);
    std::copy(longitudes, longitudes + nLongitudes, pImpl->mLongitudes.data());
    pImpl->mElevations.resize(nGrid);
    std::copy(elevation, elevation + nGrid, pImpl->mElevations.data());
    ::fixLongitudes(pImpl->mLongitudes);
    pImpl->mHaveTopography = true;
}

/// Interpolate
double Gridded::evaluate(const double latitude, const double longitudeIn) const
{
    if (latitude < -90 || latitude > 90)
    {
        throw std::invalid_argument("Latitude must be in range [-90,90]");
    }
    auto longitude = ::lonTo180(longitudeIn);
    double longitude0 = pImpl->mLongitudes.front();
    double longitude1 = pImpl->mLongitudes.back();
    double latitude0  = pImpl->mLatitudes.front();
    double latitude1  = pImpl->mLatitudes.back();
    auto xi = std::max(longitude0, std::min(longitude1, longitude));
    auto yi = std::max(latitude0,  std::min(latitude1,  latitude));
    auto lonPtr = std::lower_bound(pImpl->mLongitudes.begin(),
                                   pImpl->mLongitudes.end(),
                                   xi);
    auto iLon
        = static_cast<int> (std::distance(pImpl->mLongitudes.begin(), lonPtr)) - 1;
    iLon = std::max(0, std::min(iLon, static_cast<int> (pImpl->mLongitudes.size()) - 2));

    auto latPtr = std::lower_bound(pImpl->mLatitudes.begin(),
                                   pImpl->mLatitudes.end(),
                                   yi);
    auto iLat
        = static_cast<int> (std::distance(pImpl->mLatitudes.begin(), latPtr)) - 1;
    iLat = std::max(0, std::min(iLat, static_cast<int> (pImpl->mLatitudes.size()) - 2));
    // Interpolate
    double x0 = pImpl->mLongitudes[iLon];
    double x1 = pImpl->mLongitudes[iLon + 1];
    auto dxi = 1./(x1 - x0);
    double y0 = pImpl->mLatitudes[iLat];
    double y1 = pImpl->mLatitudes[iLat + 1];
    auto dyi = 1./(y1 - y0);
#ifndef NDEBUG
    assert(xi >= x0 && xi <= x1);
    assert(yi >= y0 && yi <= y1);
#endif
    auto nLongitudes = static_cast<int> (pImpl->mLongitudes.size());
    auto z00 = pImpl->mZ00;
    auto z01 = pImpl->mZ01;
    auto z10 = pImpl->mZ10;
    auto z11 = pImpl->mZ11;
    if (iLon != pImpl->mLon0 || iLat != pImpl->mLat0)
    {
        auto i00 = iLat*nLongitudes + iLon;
        auto i10 = i00 + 1;
        auto i01 = i00 + nLongitudes;
        auto i11 = i01 + 1;
        z00 = static_cast<double> (pImpl->mElevations.at(i00));
        z10 = static_cast<double> (pImpl->mElevations.at(i10));
        z01 = static_cast<double> (pImpl->mElevations.at(i01));
        z11 = static_cast<double> (pImpl->mElevations.at(i11));
        pImpl->mZ00 = z00;
        pImpl->mZ01 = z01;
        pImpl->mZ10 = z10;
        pImpl->mZ11 = z11;
        pImpl->mLon0 = iLon;
        pImpl->mLat0 = iLat;
    }
    auto fy0 = dxi*( (x1 - xi)*z00 + (xi - x0)*z10 );
    auto fy1 = dxi*( (x1 - xi)*z01 + (xi - x0)*z11 );
    return dyi*( (y1 - yi)*fy0 + (yi - y0)*fy1 );
}

///--------------------------------------------------------------------------///
///                            Template Instantiation                        ///
///--------------------------------------------------------------------------///
template void ULocator::Topography::Gridded::set(int, const double *,
                                                 int, const double *,
                                                 int, const double *);
template void ULocator::Topography::Gridded::set(int, const float *,
                                                 int, const float *,
                                                 int, const float *); 
