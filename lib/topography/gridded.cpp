#include <fstream>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <filesystem>
#include "uLocator/topography/gridded.hpp"
#include "uLocator/position/geographicRegion.hpp"
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
    std::vector<double> mX;
    std::vector<double> mY;
    double mZ00{0};
    double mZ01{0};
    double mZ10{0};
    double mZ11{0};
    double mMaximumElevation{0};
    double mMinimumElevation{0};
    int mIX0{-1};
    int mIY0{-1};
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
void Gridded::load(const std::string &fileName,
                   const ULocator::Position::IGeographicRegion &region)
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
        std::vector<double> latitudes, longitudes, elevations;
        ::readDataset(mFile, "Latitude",  &latitudes);
        ::checkLatitudes(latitudes.size(), latitudes.data());
        ::readDataset(mFile, "Longitude", &longitudes);
        ::checkLongitudes(longitudes.size(), longitudes.data());
        ::readDataset(mFile, "Elevation", &elevations);
        ::fixLongitudes(longitudes);
        // Make my custom (lon, lat) interpolator
        Gridded lonLatInterpolator;
        lonLatInterpolator.set(longitudes.size(), longitudes.data(),
                               latitudes.size(),  latitudes.data(),
                               elevations.size(), elevations.data());
        // Convert to a regular Cartesian grid
        auto nx = static_cast<int> (longitudes.size());
        auto ny = static_cast<int> (latitudes.size());
        std::vector<double> x(nx);
        std::vector<double> y(ny);
        auto [x0, x1] = region.getExtentInX();
        auto [y0, y1] = region.getExtentInY();
        double dx = (x1 - x0)/(nx - 1);
        double dy = (y1 - y0)/(ny - 1);
        std::vector<double> interpolatedElevations(elevations.size(), 0);
//std::ofstream ofl;
//ofl.open("interpElevation.txt");
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix) 
            {
                double xi = x0 + ix*dx;
                double yi = y0 + iy*dy;
                x[ix] = xi;
                y[iy] = yi;
                auto indx = iy*nx + ix;
                auto [latitude, longitude]
                    = region.localToGeographicCoordinates(xi, yi);
                interpolatedElevations[indx] = lonLatInterpolator(longitude, latitude);
            }
//ofl << std::endl;
        }
//ofl.close();
        set(x.size(), x.data(),
            y.size(), y.data(),
            interpolatedElevations.size(), interpolatedElevations.data());
/*
std::ofstream ofl;
ofl.open("recoveredElevation.txt");
        for (int iy = 0; iy < ny; ++iy)
        {
            for (int ix = 0; ix < nx; ++ix)
            {
                auto [xi, yi] = region.geographicToLocalCoordinates(latitudes[iy], longitudes[ix]);
//                auto elevation = evaluate(xi, yi);
                auto elevation = evaluate(x[ix], y[iy]);
ofl << x[ix] << " " << y[iy] << " " << elevation << std::endl;
            }
ofl << std::endl;
        }
ofl.close();
*/
    }
    catch (const std::exception &e)
    {
        pImpl->mX.clear();
        pImpl->mY.clear();
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
void Gridded::set(const int nx, const U *xs,
                  const int ny, const U *ys,
                  const int nGrid, const U *elevation)
{
    if (nx*ny != nGrid)
    {
        throw std::invalid_argument("nGrid != nx x ny");
    }
    if (elevation == nullptr)
    {
        throw std::invalid_argument("Elevation is NULL");
    }
    if (!std::is_sorted(xs, xs + nx, std::less<>{}))
    {
        throw std::invalid_argument("xs must be sorted");
    }
    if (!std::is_sorted(ys, ys + ny, std::less<>{}))
    {
        throw std::invalid_argument("ys must be sorted");
    }
    //::checkLatitudes(nLatitudes, latitudes);
    //::checkLongitudes(nLongitudes, longitudes);
    // Copy
    pImpl->mX.resize(nx);
    std::copy(xs, xs + nx, pImpl->mX.data());
    pImpl->mY.resize(ny);
    std::copy(ys, ys + ny, pImpl->mY.data());
    pImpl->mElevations.resize(nGrid);
    std::copy(elevation, elevation + nGrid, pImpl->mElevations.data());
    const auto [elevationMin, elevationMax]
       = std::minmax_element(pImpl->mElevations.begin(), pImpl->mElevations.end());
    pImpl->mMinimumElevation = *elevationMin;
    pImpl->mMaximumElevation = *elevationMax;
    pImpl->mHaveTopography = true;
}

/// Interpolate
double Gridded::evaluate(const double x, const double y) const
{
    return evaluate(x, y, nullptr, nullptr);
}

/// Interpolate
double Gridded::evaluate(const double x, const double y,
                         double *dElevationDx, double *dElevationDy) const
{
    double xOrigin = pImpl->mX.front();
    double xExtent = pImpl->mX.back();
    double yOrigin = pImpl->mY.front();
    double yExtent = pImpl->mY.back();
    auto xi = std::max(xOrigin, std::min(xExtent, x));
    auto yi = std::max(yOrigin, std::min(yExtent,  y));
    auto xPtr = std::lower_bound(pImpl->mX.begin(), pImpl->mX.end(), xi);
    auto ix = static_cast<int> (std::distance(pImpl->mX.begin(), xPtr)) - 1;
    ix = std::max(0, std::min(ix, static_cast<int> (pImpl->mX.size()) - 2));

    auto yPtr = std::lower_bound(pImpl->mY.begin(), pImpl->mY.end(), yi);
    auto iy = static_cast<int> (std::distance(pImpl->mY.begin(), yPtr)) - 1;
    iy = std::max(0, std::min(iy, static_cast<int> (pImpl->mY.size()) - 2));
    // Interpolate
    double x0 = pImpl->mX[ix];
    double x1 = pImpl->mX[ix + 1];
    auto dxi = 1./(x1 - x0);
    double y0 = pImpl->mY[iy];
    double y1 = pImpl->mY[iy + 1];
    auto dyi = 1./(y1 - y0);
#ifndef NDEBUG
    assert(xi >= x0 && xi <= x1);
    assert(yi >= y0 && yi <= y1);
#endif
    auto nx = static_cast<int> (pImpl->mX.size());
    auto z00 = pImpl->mZ00;
    auto z01 = pImpl->mZ01;
    auto z10 = pImpl->mZ10;
    auto z11 = pImpl->mZ11;
    if (ix != pImpl->mIX0 || iy != pImpl->mIY0)
    {
        auto i00 = iy*nx + ix;
        auto i10 = i00 + 1;
        auto i01 = i00 + nx;
        auto i11 = i01 + 1;
        z00 = static_cast<double> (pImpl->mElevations.at(i00));
        z10 = static_cast<double> (pImpl->mElevations.at(i10));
        z01 = static_cast<double> (pImpl->mElevations.at(i01));
        z11 = static_cast<double> (pImpl->mElevations.at(i11));
        pImpl->mZ00 = z00;
        pImpl->mZ01 = z01;
        pImpl->mZ10 = z10;
        pImpl->mZ11 = z11;
        pImpl->mIX0 = ix;
        pImpl->mIY0 = iy;
    }
    auto fy0 = dxi*( (x1 - xi)*z00 + (xi - x0)*z10 );
    auto fy1 = dxi*( (x1 - xi)*z01 + (xi - x0)*z11 );
    if (dElevationDx != nullptr)
    {
        auto dfy0dx = dxi*( z10 - z00 );
        auto dfy1dx = dxi*( z11 - z01 );
        *dElevationDx = dyi*( (y1 - yi)*dfy0dx + (yi - y0)*dfy1dx ); 
    }
    if (dElevationDy != nullptr)
    {
        *dElevationDy = dyi*(fy1 - fy0);
    }
    return dyi*( (y1 - yi)*fy0 + (yi - y0)*fy1 );
}

/// Min/max elevation
std::pair<double, double> Gridded::getMinimumAndMaximumElevation() const
{
    if (!haveTopography()){throw std::runtime_error("Topography not set");}
    return std::pair {pImpl->mMinimumElevation, pImpl->mMaximumElevation};
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
