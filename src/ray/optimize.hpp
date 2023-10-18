#ifndef PRIVATE_OPTIMIZE_HPP
#define PRIVATE_OPTIMIZE_HPP
/// @brief Utilities for optimizing ray paths in ray shooting.
#include <string>
#include <cmath>
#include <vector>
#include <functional>
#ifndef NDEBUG
#include <cassert>
#endif
#include <boost/math/tools/roots.hpp>
#include <boost/math/tools/minima.hpp>
#include "layerTracer.hpp"
#include "eikonalxx/ray/path2d.hpp"
namespace
{

struct IObjectiveFunction 
{
    virtual double operator()(const double &) const = 0;
};

struct DirectObjectiveFunction : public IObjectiveFunction
{
    DirectObjectiveFunction(const std::vector<double> &augmentedInterfacesIn,
                            const std::vector<double> &augmentedSlownessesIn,
                            const double sourceDepthIn,
                            const double stationDepthIn,
                            const double stationOffsetIn,
                            const int sourceLayerIn,
                            const int stationLayerIn) :
        augmentedInterfaces(augmentedInterfacesIn),
        augmentedSlownesses(augmentedSlownessesIn),
        sourceDepth(sourceDepthIn),
        stationDepth(stationDepthIn),
        stationOffset(stationOffsetIn),
        sourceLayer(sourceLayerIn),
        stationLayer(stationLayerIn)
    {
    }

    double operator()(const double &takeOffAngle) const override
    {
        nEvaluations = nEvaluations + 1;
        constexpr bool allowCriticalRefractions{false};
        constexpr bool keepOnlyHits{false};
        constexpr double rayHitTolerance{1};
        int lastLayer =-1; // Doesn't matter
        auto paths = ::shoot(takeOffAngle,
                             augmentedInterfaces,
                             augmentedSlownesses,
                             sourceLayer,
                             sourceDepth,
                             stationLayer,
                             stationOffset,
                             stationDepth,
                             allowCriticalRefractions,
                             rayHitTolerance,
                             keepOnlyHits,
                             lastLayer);
        if (!paths.empty())
        {
            // Otherwise the fastest ray
            if (paths.size() > 1 && std::abs(takeOffAngle - 180) > 1.e-10)
            {
                std::cerr << "Should be only 1 direct ray; size = "
                          << paths.size() << std::endl;
            }
            auto nSegments = paths[0].size();
            const auto lastSegment = paths[0].at(nSegments - 1); 
            auto xOffset = stationOffset
                         - lastSegment.getEndPoint().getPositionInX();
            if (takeAbsoluteValue){xOffset = std::abs(xOffset);}
            if (keepRayPath){rayPath = std::move(paths[0]);}
            return xOffset;
        }
        throw std::runtime_error("No path for direct wave with angle "
                               + std::to_string(takeOffAngle));
    }
    mutable EikonalXX::Ray::Path2D rayPath;
    std::vector<double> augmentedInterfaces;
    std::vector<double> augmentedSlownesses;
    double sourceDepth;
    double stationDepth;
    double stationOffset;
    int sourceLayer;
    int stationLayer;
    bool takeAbsoluteValue{false};
    bool keepRayPath{false};
    mutable int nEvaluations{0};
};

struct ReflectionObjectiveFunction : public IObjectiveFunction
{
    ReflectionObjectiveFunction(const std::vector<double> &augmentedInterfacesIn,
                                const std::vector<double> &augmentedSlownessesIn,
                                const double sourceDepthIn,
                                const double stationDepthIn,
                                const double stationOffsetIn,
                                const int sourceLayerIn,
                                const int stationLayerIn) :
        augmentedInterfaces(augmentedInterfacesIn),
        augmentedSlownesses(augmentedSlownessesIn),
        sourceDepth(sourceDepthIn),
        stationDepth(stationDepthIn),
        stationOffset(stationOffsetIn),
        sourceLayer(sourceLayerIn),
        stationLayer(stationLayerIn)
    {
        lastLayer = static_cast<int> (augmentedInterfaces.size()) - 2;
    }   

    double operator()(const double &takeOffAngle) const override
    {   
        nEvaluations = nEvaluations + 1;
        constexpr bool allowCriticalRefractions{false};
        constexpr bool keepOnlyHits{false};
        constexpr double rayHitTolerance{1};
        auto paths = ::shoot(takeOffAngle,
                             augmentedInterfaces,
                             augmentedSlownesses,
                             sourceLayer,
                             sourceDepth,
                             stationLayer,
                             stationOffset,
                             stationDepth,
                             allowCriticalRefractions,
                             rayHitTolerance,
                             keepOnlyHits,
                             lastLayer);
        if (!paths.empty())
        {
            auto nSegments = paths[0].size();
            const auto lastSegment = paths[0].at(nSegments - 1);
            auto xOffset = stationOffset
                         - lastSegment.getEndPoint().getPositionInX();
            if (takeAbsoluteValue){xOffset = std::abs(xOffset);}
            if (keepRayPath){rayPath = std::move(paths[0]);}
            return xOffset;
        }
        throw std::runtime_error("No paths for reflected wave");
    }
    mutable EikonalXX::Ray::Path2D rayPath;
    std::vector<double> augmentedInterfaces;
    std::vector<double> augmentedSlownesses;
    double sourceDepth;
    double stationDepth;
    double stationOffset;
    int sourceLayer;
    int stationLayer;
    int lastLayer;
    bool takeAbsoluteValue{false};
    bool keepRayPath{false};
    mutable int nEvaluations{0};
};


/// @result 
void computeTakeOffAngles(
    const int n, const double theta0, const double theta1,
    std::vector<double> *takeOffAngles)
{
#ifndef NDEBUG
    assert(n > 0);
#endif
    if (static_cast<int> (takeOffAngles->size()) != n)
    {
        takeOffAngles->resize(n);
    }
    if (n == 1)
    {
        takeOffAngles->at(0) = 0.5*(theta0 + theta1);
    }
    auto dTheta = (theta1 - theta0)/(n - 1);
    auto *__restrict__ takeOffAnglePointer = takeOffAngles->data();
    for (int i = 0; i < n; ++i)
    {
        takeOffAnglePointer[i] = theta0 + dTheta*i;
    }
}

/*
[[nodiscard]]
/// @param[in] shootRay  A function, that given a take-off angle, shoots a
///                      ray and computes the offset.
/// @param[in] fa  The left-bracketing take-off angle and offset.
/// @param[in] fb  An intermediate take-off angle and offset.
/// @param[in] fc  The right-bracketing take-off angle and offset.
/// @param[in] nIterations  The number of iterations in the optimization.
/// @result The optimal take-off angle and ray path.
double quadraticFitSearch(const IObjectiveFunction &objectiveFunction,
                          const std::pair<double, double> &fa,
                          const std::pair<double, double> &fb,
                          const std::pair<double, double> &fc,
                          const int nIterations = 7)
{
    constexpr double third{1./3.};
    auto a  = fa.first;
    auto ya = fa.second;
    auto b  = fb.first;
    auto yb = fb.second;
    auto c  = fc.first;
    auto yc = fc.second;
    auto xStar = third*(a + b + c);
    auto fStar = std::max(ya, std::max(yb, yc));
    for (int k = 0; k < nIterations; ++k)
    {
#ifndef NDEBUG
        assert(a < b);
        assert(b < c);
#endif
        auto a2 = a*a;
        auto b2 = b*b;
        auto c2 = c*c;
        auto xNum = ya*(b2 - c2) + yb*(c2 - a2) + yc*(a2 - b2);
        auto xDen = 2*(ya*(b - c) + yb*(c - a) +  yc*(a - b));
        xStar = third*(a + b + c);
        if (std::abs(xDen) > 1.e-14)
        {
            xStar = xNum/xDen;
        }
        else
        {
            std::cerr << "Division by 0: " << xNum << ","
                      << xDen << std::endl;
            if (ya < yb && ya < yc){return a;}
            if (yb < ya && yb < yc){return b;}
            return c;
        }
        try
        {
            fStar = objectiveFunction(xStar); 
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            if (ya < yb && ya < yc){return a;} 
            if (yb < ya && yb < yc){return b;} 
            return c;
        }
        if (xStar > b)
        {
            if (fStar > yb) 
            {
                c = xStar;
                yc = fStar; 
            }
            else
            {
                a  = b;
                ya = yb; 
                b  = xStar;
                yb = fStar;
            }   
        }
        else if (xStar < b)
        {
            if (fStar > yb) 
            {
                a = xStar;
                ya = fStar;
            }
            else
            {
                c  = b;
                yc = yb; 
                b  = xStar;
                yb = fStar;
            }
        }
        else
        {
            return b;
        }
    }
    return xStar;
}
*/

/// @brief Optimize direct arriving waves.
std::vector<EikonalXX::Ray::Path2D>
    optimizeDirect(
    const std::vector<double> &augmentedInterfaces,
    const std::vector<double> &augmentedSlownesses,
    const int sourceLayer,
    const double sourceDepth,
    const int stationLayer,
    const double stationDepth,
    const double stationOffset,
    const double rayHitTolerance,
    const int nBits = 30)
{
    std::vector<EikonalXX::Ray::Path2D> rayPaths;
    // Source layer = stations layer -> easy
    if (sourceLayer == stationLayer)
    {
        EikonalXX::Ray::Path2D path;
#ifndef NDEBUG
        auto result =
#endif
        ::traceDirectSameLayer(augmentedInterfaces,
                               augmentedSlownesses,
                               sourceLayer,
                               sourceDepth,
                               stationLayer,
                               stationOffset,
                               stationDepth,
                               &path);
#ifndef NDEBUG
        assert(result == ReturnCode::Hit);
#endif
        rayPaths.push_back(path);
        return rayPaths;
    }
    // Flip the problem around
    if (stationDepth > sourceDepth)
    {
        //std::cout << "Station below source for direct" << std::endl;
        rayPaths = ::optimizeDirect(augmentedInterfaces,
                                    augmentedSlownesses,
                                    stationLayer,
                                    stationDepth,
                                    sourceLayer,
                                    sourceDepth,
                                    stationOffset,
                                    rayHitTolerance,
                                    nBits);
        for (auto &rayPath : rayPaths)
        {
            rayPath.reverse();
        }
        return rayPaths;
    }
    // Initialize objective function
    DirectObjectiveFunction objectiveFunction{augmentedInterfaces,
                                              augmentedSlownesses,
                                              sourceDepth,
                                              stationDepth,
                                              stationOffset,
                                              sourceLayer,
                                              stationLayer};
    objectiveFunction.keepRayPath = false;
    // Make a bunch of take-off angles
    constexpr int nInitialAngles{90};
    constexpr double thetaMin{ 90.0001};
    constexpr double thetaMax{180};
    std::vector<double> takeOffAngles;
    ::computeTakeOffAngles(nInitialAngles, thetaMin, thetaMax, &takeOffAngles);
    // Shoot these rays
    std::vector<std::pair<double, double>> angleOffsets(nInitialAngles);
/*
    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<size_t> (0, takeOffAngles.size()),
        [&](const oneapi::tbb::blocked_range<size_t> &range)
    {
        for (size_t i = range.begin(); i != range.end(); ++i)
        {
            auto takeOffAngle = takeOffAngles[i];
            try
            {
                auto offset = objectiveFunction(takeOffAngle);
                angleOffsets[i] = std::pair{takeOffAngle, offset};
            }
            catch (const std::exception &e)
            {
                angleOffsets[i]
                    = std::pair {takeOffAngle,
                                 std::numeric_limits<double>::max()};
                std::cerr << "Initial parallel shooting with take-off angle: "
                          << takeOffAngles[i] << " failed with: " 
                          << e.what() << std::endl;

            }
        }
    });
*/
    for (int i = 0; i < nInitialAngles; ++i)
    {
        auto takeOffAngle = takeOffAngles[i];
        try
        {
            auto offset = objectiveFunction(takeOffAngle);
            angleOffsets[i] = std::pair{takeOffAngle, offset};
        }
        catch (const std::exception &e)
        {
            angleOffsets[i] = std::pair {takeOffAngle,
                                         std::numeric_limits<double>::max()};
            std::cerr << "Initial shooting with take-off angle: " 
                      << takeOffAngles[i] << " failed with: " 
                      << e.what() << std::endl;
        }
    }
    // Remove the failures
    angleOffsets.erase(
        std::remove_if(angleOffsets.begin(), angleOffsets.end(),
                       [](const std::pair<double, double> &o)
                       {
                          return o.second == std::numeric_limits<double>::max();
                       }),
        angleOffsets.end());
    // Nothing hit
    if (angleOffsets.empty()){return rayPaths;}
    // Now we optimize in the bracketed solution(s)
    auto nPaths = static_cast<int> (angleOffsets.size());
    objectiveFunction.takeAbsoluteValue = true; // Switch to minimization
    objectiveFunction.keepRayPath = true; 
    for (int iPath = 0; iPath < nPaths - 1; ++iPath)
    {
        if (std::copysign(1.0, angleOffsets[iPath].second) !=
            std::copysign(1.0, angleOffsets[iPath + 1].second))
        {
            std::uintmax_t maxIterations{1000};
            auto takeOffAngle0 = angleOffsets[iPath].first;
            auto takeOffAngle1 = angleOffsets[iPath + 1].first;
            try
            {
                std::pair<double, double> result
                    = boost::math::tools::brent_find_minima(objectiveFunction,
                                                        takeOffAngle0,
                                                        takeOffAngle1,
                                                        nBits,
                                                        maxIterations);
                 //std::cout << takeOffAngle0 << " " << result.first << " " << takeOffAngle1 << " | " << result.second << " " << maxIterations << std::endl;
                 objectiveFunction(result.first);
                 //std::cout << objectiveFunction.rayPath.getTravelTime() << std::endl;
                 rayPaths.push_back(objectiveFunction.rayPath);
            }
            catch (const std::exception &e)
            {
                 std::cerr << e.what() << std::endl;
            }
        }
    }
    return rayPaths;
}

/// @brief Optimizes the critically refracted arrivals.  
/// @note This assumes the velocities increase with depth.
std::vector<EikonalXX::Ray::Path2D>
    optimizeCriticallyRefracted(
    const std::vector<double> &augmentedInterfaces,
    const std::vector<double> &augmentedSlownesses,
    const int sourceLayer,
    const double sourceDepth,
    const int stationLayer,
    const double stationDepth,
    const double stationOffset)
{
    constexpr double rayHitTolerance{1}; // 1 meter
    auto nLayers = static_cast<int> (augmentedSlownesses.size()) - 1;
    // Check for velocity inversion
    for (int iLayer = 0; iLayer < nLayers - 1; ++iLayer)
    {   
        if (augmentedSlownesses.at(iLayer) < augmentedSlownesses.at(iLayer + 1)) 
        {
            std::cerr << "Velocity inversion exists!" << std::endl;
        }
    }
    std::vector<EikonalXX::Ray::Path2D> rayPaths;
    // Flip the problem around
    if (stationDepth > sourceDepth)
    {
        //std::cout << "Station below source for critical refraction" << std::endl;
        rayPaths = ::optimizeCriticallyRefracted(augmentedInterfaces,
                                                 augmentedSlownesses,
                                                 stationLayer,
                                                 stationDepth,
                                                 sourceLayer,
                                                 sourceDepth,
                                                 stationOffset);
        for (auto &rayPath : rayPaths)
        {
            rayPath.reverse();
        }
        return rayPaths;
    }
    // The idea here is pretty simple.  Basically, the trick is to spend
    // the least amount of time possible in a slower layer.  To do that
    // we use the smallest angle possible - which is the critical angle.
    // Then we just use some algebra to relate that angle for each layer
    // to the other layers.  Since:
    //    sin(i_c) = V_{i}/V_{i+1} -> i_c = asin(V_{i}/V_{i+1})
    // Then what's this angle in the previous layer?
    //    sin(i)/V_{i-1} = sin(i_c)/V_i
    //                   = (V_i/V_{i+1})/V_i
    //                   = 1/V{i+1}
    //    sin(i)/V_{i-1} = 1/V_{i+1}
    // -> i = asin(V_{i-1}/V_{i+1})
    // And following the recursion down we get, in general, for a velocity
    // model where the velocities increase with depth:
    // -> i = asin(V_{sourceLayer}/V_{i+1})
    // However, more generally, I've found it advantageous to simply
    // try every critical angle.  So that's what we do assuming the
    // velocities increase with depth.
    std::vector<double> takeOffAngles;
    auto lastLayer = nLayers - 1;
    takeOffAngles.reserve(nLayers*(nLayers + 1));
/*
    auto startLayer = std::min(sourceLayer, stationLayer);
    for (int layer = startLayer + 1; layer < nLayers; ++layer)
    {
        if (augmentedSlownesses[startLayer] < augmentedSlownesses[layer])
        {
            std::cerr << "Velocity inversion exists!" << std::endl;
        }
        auto angle = ::computeCriticalAngle(augmentedSlownesses[startLayer],
                                            augmentedSlownesses[layer]);
        // Pad just a touch in case we get burned by inexact math
        takeOffAngles.push_back(angle*(180./M_PI) + 1.e-10);
    }
*/
    for (int iLayer = 0; iLayer < nLayers - 1; ++iLayer)
    {
        for (int jLayer = iLayer + 1; jLayer < nLayers; ++jLayer)
        {

            auto angle = ::computeCriticalAngle(augmentedSlownesses[iLayer],
                                                augmentedSlownesses[jLayer]);
            // Pad just a touch in case we get burned by inexact math
            takeOffAngles.push_back(angle*(180./M_PI) + 1.e-10);
        }
    }
/*
    for (int layer = 1; layer < nLayers; ++layer)
    {
        auto angle = ::computeCriticalAngle(augmentedSlownesses[layer - 1],
                                            augmentedSlownesses[layer]);
        // Pad just a touch in case we get burned by inexact math
        takeOffAngles.push_back(angle*(180./M_PI) + 1.e-10);
    }
*/
    std::sort(takeOffAngles.begin(), takeOffAngles.end());
    takeOffAngles.erase(std::unique(takeOffAngles.begin(),
                                    takeOffAngles.end()), takeOffAngles.end());
    // Shoot these rays
    for (const auto &takeOffAngle : takeOffAngles)
    {
        //std::cout << takeOffAngle << std::endl;
        constexpr bool allowCriticalRefractions{true};
        constexpr bool keepOnlyHits{true};
        auto paths = ::shoot(takeOffAngle,
                             augmentedInterfaces,
                             augmentedSlownesses,
                             sourceLayer,
                             sourceDepth,
                             stationLayer,
                             stationOffset,
                             stationDepth,
                             allowCriticalRefractions,
                             rayHitTolerance,
                             keepOnlyHits,
                             lastLayer);
        if (!paths.empty())
        {
            rayPaths.insert(rayPaths.end(), paths.begin(), paths.end());
        }
    }
    // Sort them
    if (!rayPaths.empty())
    {
        std::sort(rayPaths.begin(), rayPaths.end(),
                  [](const EikonalXX::Ray::Path2D &lhs,
                     const EikonalXX::Ray::Path2D &rhs)
                 {
                     return lhs.getTravelTime() < rhs.getTravelTime();
                 });
        //std::cout << rayPaths[0].getTravelTime() << std::endl;
    }
    return rayPaths;
}
    
/// Optimizes downgoing waves that reflect
std::vector<EikonalXX::Ray::Path2D>
    optimizeFirstArrivingDownGoingReflections(
    const std::vector<double> &augmentedInterfaces,
    const std::vector<double> &augmentedSlownesses,
    const int sourceLayer,
    const double sourceDepth,
    const int stationLayer,
    const double stationDepth,
    const double stationOffset,
    const double rayHitTolerance,
    const int nBits = 30)
{
    std::vector<EikonalXX::Ray::Path2D> rayPaths;
    // Flip the problem around
    if (stationDepth > sourceDepth)
    {
        //std::cout << "Station below source for downgoing" << std::endl;
        rayPaths
            = ::optimizeFirstArrivingDownGoingReflections(augmentedInterfaces,
                                                          augmentedSlownesses,
                                                          stationLayer,
                                                          stationDepth,
                                                          sourceLayer,
                                                          sourceDepth,
                                                          stationOffset,
                                                          rayHitTolerance,
                                                          nBits);
        for (auto &rayPath : rayPaths)
        {
            rayPath.reverse();
        }
        return rayPaths;
    }
    // Initialize objective function
    ReflectionObjectiveFunction objectiveFunction{augmentedInterfaces,
                                                  augmentedSlownesses,
                                                  sourceDepth,
                                                  stationDepth,
                                                  stationOffset,
                                                  sourceLayer,
                                                  stationLayer};
    objectiveFunction.keepRayPath = false;
    // Make a bunch of take-off angles
    constexpr int nInitialAngles{90};
    constexpr double thetaMin{0};
    constexpr double thetaMax{90.0001};
    std::vector<double> takeOffAngles;
    ::computeTakeOffAngles(nInitialAngles, thetaMin, thetaMax, &takeOffAngles);
    // Shoot these rays
    std::vector<std::pair<double, double>> angleOffsets;
    angleOffsets.reserve(takeOffAngles.size());
    for (const auto &takeOffAngle : takeOffAngles)
    {
        try
        {
            angleOffsets.push_back(std::pair {takeOffAngle,
                                              objectiveFunction(takeOffAngle)});
        }
        catch (const std::exception &e) 
        {
            std::cerr << e.what() << std::endl;
        }
    }
    // Nothing hit
    if (angleOffsets.empty()){return rayPaths;}
    // Now we optimize in the bracketed solution(s)
    auto nPaths = static_cast<int> (angleOffsets.size());
    objectiveFunction.takeAbsoluteValue = true; // Switch to minimization
    objectiveFunction.keepRayPath = true;
    for (int iPath = 0; iPath < nPaths - 1; ++iPath)
    {
        if (std::copysign(1.0, angleOffsets[iPath].second) !=
            std::copysign(1.0, angleOffsets[iPath + 1].second))
        {
            std::uintmax_t maxIterations{1000};
            auto takeOffAngle0 = angleOffsets[iPath].first;
            auto takeOffAngle1 = angleOffsets[iPath + 1].first;
            try
            {
                std::pair<double, double> result
                    = boost::math::tools::brent_find_minima(objectiveFunction,
                                                            takeOffAngle0,
                                                            takeOffAngle1,
                                                            nBits,
                                                            maxIterations);
                 objectiveFunction(result.first);
                 rayPaths.push_back(objectiveFunction.rayPath);
            }
            catch (const std::exception &e)
            {
                 std::cerr << e.what() << std::endl;
            }
        }
    }
    return rayPaths;
}

}
#endif
