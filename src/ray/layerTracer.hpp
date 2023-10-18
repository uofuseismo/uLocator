#ifndef PRIVATE_LAYER_TRACER_HPP
#define PRIVATE_LAYER_TRACER_HPP
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#ifndef NDEBUG
#include <cassert>
#endif
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"

namespace
{

enum class ReturnCode
{
    UnderShot =-1,
    Hit = 0,
    OverShot =+1,
    RayDoesNotTurn,
    RayTurnsTooEarly
};

struct Segment
{
    double slowness{0};
    double x0{0};
    double z0{0};
    double x1{0};
    double z1{0};
    int layer{0};
    void reverse()
    {
        std::swap(x0, x1);
        std::swap(z0, z1);
    }
    [[nodiscard]] EikonalXX::Ray::Segment2D toSegment() const
    {
        EikonalXX::Ray::Point2D startPoint{x0, z0};
        EikonalXX::Ray::Point2D endPoint{x1, z1};
        EikonalXX::Ray::Segment2D segment;
        segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
        segment.setSlowness(slowness);
        segment.setVelocityModelCellIndex(layer);
        return segment;
    }
    [[nodiscard]] double getTravelTime() const
    {
        return std::hypot(x1 - x0, z1 - z0)*slowness;
    }
};

/// Reverses a segment
void reverseSegments(std::vector<::Segment> *segments)
{
    std::reverse(segments->begin(), segments->end());
    for (auto &segment : *segments)
    {
        segment.reverse();
    }
}

/// @brief Converts the local ray path to the library variant whilst
///        performing any merging.
[[nodiscard]]
EikonalXX::Ray::Path2D
    toRayPath(const std::vector<::Segment> &segments)
{
    EikonalXX::Ray::Path2D rayPath;
    constexpr double tolerance{0.001}; // 1 millimeter
    if (segments.empty()){throw std::runtime_error("Segments is empty");}
    auto nSegments = static_cast<int> (segments.size());
    std::vector<EikonalXX::Ray::Segment2D> raySegments;
    raySegments.reserve(nSegments);
    // Handle edge case
    if (nSegments == 1)
    {
        raySegments.push_back(segments[0].toSegment());
        rayPath.set(std::move(raySegments));
        return rayPath;
    }
    // Merge ray paths in same layer
    for (int i = 0; i < nSegments; ++i)
    {
        ::Segment segment{segments[i]};
        std::array<double, 2> vi{segment.x1 - segment.x0,
                                 segment.z1 - segment.z0};
        auto viLength = std::hypot(vi[0], vi[1]);
        int iNext = 0;
        // Merge chunk of segments
        for (int j = i + 1; j < nSegments; ++j)
        {
            std::array<double, 2> vj{segments[j].x1 - segments[j].x0,
                                     segments[j].z1 - segments[j].z0};
            // 0-length segment
            auto vjLength = std::hypot(vj[0], vj[1]);
            if (segment.layer == segments[j].layer)
            {
#ifndef NDEBUG
                assert(std::abs(segment.slowness - segments[j].slowness) < 
                       1.e-10);
#endif
                /*
                // 0-length segment in same layer - merge it
                if (vjLength < tolerance ||
                    (j == i + 1 && viLength < tolerance))
                {
                    segment.x1 = segments[j].x1;
                    segment.z1 = segments[j].z1;
                    iNext = iNext + 1;
                    continue;
                }
                */
                // Parallel lines in same layer:
                //    u.v/(|u||v|) = 1 -> |u||v| = u.v
                // Note, we don't want -1 b/c that is a vertical refrlection.
                // Additionally, this works for a zero-length segment.
                auto dotProduct = vi[0]*vj[0] + vi[1]*vj[1];
                if (std::abs(viLength*vjLength - dotProduct) < tolerance)
                {
                    segment.x1 = segments[j].x1;
                    segment.z1 = segments[j].z1;
                    iNext = iNext + 1;
                    continue;
                }
                else
                {
                    break;
                }
            }
            break;
        }
        raySegments.push_back(segment.toSegment());
        i = i + iNext;
    }
/*
    raySegments.reserve(segments.size());
    for (const auto &segment : segments)
    {   
        raySegments.push_back(segment.toSegment());
    }
*/
    rayPath.set(std::move(raySegments));
    return rayPath;
}

[[nodiscard]] int getLayer(const double depth,
                           const std::vector<double> &interfaces,
                           const bool isAugmented = true)
{    
    if (interfaces.size() == 1){return 0;}
#ifndef NDEBUG
    if (isAugmented){assert(interfaces.size() > 1);}
#endif
    auto nLayers = static_cast<int> (interfaces.size());
    if (isAugmented){nLayers = nLayers - 1;}
    int layer = 0; 
    if (depth < interfaces.front())
    {
        layer = 0;
    }
    else if (depth >= interfaces.back())
    {
        layer = nLayers - 1;
    }
    else
    {
        layer = std::distance(
                   interfaces.begin(),
                   std::upper_bound(interfaces.begin(),
                                    interfaces.begin() + nLayers,
                                    depth)) - 1; 
#ifndef NDEBUG
        assert(depth >= interfaces.at(layer) &&
               depth < interfaces.at(layer + 1)); 
#endif
    }
#ifndef NDEBUG
    assert(layer >= 0);
    if (nLayers > 1){assert(layer < nLayers);} 
#endif
    return layer;
}   

/// @brief Computes the critical angle from the slownesses.
/// @param[in] s0   The slowness in the current layer in s/m.
/// @param[in] s1   The slowness in the transmission layer in s/m.
/// @result The critical angle in radians.
/// @note Cerveny, Page 121
[[nodiscard]] double computeCriticalAngle(const double s0,
                                          const double s1)
{
    if (s0 <= s1){return M_PI_2;} // Use max transmission angle
    // arcsin(v0/v1) = arcsin( (1/s0)/(1/s1) ) = arcsin(s1/s0)
    return std::asin(s1/s0);
}

/// @brief Given the incidence angle and the slowness in this layer and the
///        next layer, this computes the corresponding transmission angle.
/// @param[in] incidenceAngle  The angle of incidence in radians.
/// @param[in] s0              The slowness in the current layer in s/m.
/// @param[in] s1              The slowness in the transmission laye rin s/m.
/// @result The transmission angle in radians.
/// @note Cerveny pg 42
[[nodiscard]] double computeTransmissionAngle(const double incidenceAngle,
                                              const double s0,
                                              const double s1)
{
#ifndef NDEBUG
    assert(incidenceAngle >= 0 && incidenceAngle <= M_PI);
#endif
    // sin(i0)/v0 = sin(i1)/v1 -> i1 = arcsin((s0/s1)*sin(i0))
    auto asinArgument = (s0/s1)*std::sin(incidenceAngle);
    if (asinArgument >= 0 && asinArgument <= 1)
    {
        return std::asin(asinArgument);
    }
    return -incidenceAngle; // Post critical -> total internal reflection
}

/// @result Defines the ray convergence
[[nodiscard]]
ReturnCode checkRayConvergence(const std::vector<::Segment> &segments,
                               const double stationOffset,
                               const double tolerance = 1)
{
    if (std::abs(segments.back().x1 - stationOffset) < tolerance)
    {
        return ReturnCode::Hit;
    }
    else
    {
        if (segments.back().x1 > stationOffset)
        {
            return ReturnCode::OverShot;
        }
        else
        {
            return ReturnCode::UnderShot;
        }
    }
#ifndef NDEBUG
    assert(false);
#endif
}

/// @brief Traces a ray strictly down through a stack of layers.
ReturnCode traceVerticalReflectionDown(const std::vector<double> &interfaces,
                                       const std::vector<double> &slownesses,
                                       const int sourceLayer,
                                       const int endLayer,
                                       const double sourceDepth,
                                       const int stationLayer,
                                       const double stationDepth,
                                       const double stationOffset,
                                       std::vector<::Segment> *segments,
                                       const double tolerance = 1)
{
    segments->clear();
    auto nLayers = static_cast<int> (interfaces.size()) - 1;
#ifndef NDEBUG
    assert(segments != nullptr);
    assert(interfaces.size() == slownesses.size());
    //if (stationDepth > sourceDepth){assert(stationLayer <= endLayer);}
    //assert(endLayer < nLayers - 1); // Can't bounce from whole-space
    assert(stationOffset >= 0); 
    //assert(sourceDepth >= interfaces[sourceLayer] &&
    //       sourceDepth <  interfaces[sourceLayer + 1]);
    //assert(stationDepth >= interfaces[stationLayer] &&
    //       stationDepth <  interfaces[stationLayer + 1]);
#endif
    if (stationDepth > sourceDepth)
    {
        auto returnCode = ::traceVerticalReflectionDown(interfaces,
                                                        slownesses,
                                                        stationLayer,
                                                        endLayer,
                                                        stationDepth,
                                                        sourceLayer,
                                                        sourceDepth,
                                                        stationOffset,
                                                        segments,
                                                        tolerance);
        std::reverse(segments->begin(), segments->end());
        for (auto &segment : *segments){segment.reverse();}
        return returnCode;
    }
#ifndef NDEBUG
    assert(sourceLayer >= 0 && sourceLayer <= endLayer);
    assert(sourceDepth >= stationDepth);  
    assert(stationLayer >= 0);
    assert(stationLayer <= endLayer);
    assert(endLayer < nLayers - 1); // Can't bounce from whole-space
    assert(stationLayer <= endLayer);
    assert(sourceDepth >= interfaces[sourceLayer] &&
           sourceDepth <  interfaces[sourceLayer + 1]);
    assert(stationDepth >= interfaces[stationLayer] &&
           stationDepth <  interfaces[stationLayer + 1]);
#endif
/*
    // N.B. don't do the direct ray b/c the function name says `reflection'
    if (sourceLayer == stationLayer && sourceLayer == endLayer)
    {
        ::Segment segment{slownesses[sourceLayer],
                          x0, z0, x1, stationDepth,
                          sourceLayer};
        segments->push_back(segment);
        return ::checkRayConvergence(*segments, stationOffset, tolerance);
    } 
*/
    if (segments->capacity() < 2*interfaces.size() + 1)
    {
        segments->reserve(2*interfaces.size() + 1);
    }
    // Trace from source down
    constexpr double x0{0};
    constexpr double x1{0};
    double z0{sourceDepth};
    for (int i = sourceLayer; i <= endLayer; ++i)
    {
        double z1 = interfaces[i + 1]; 
        ::Segment segment{slownesses[i],
                          x0, z0, x1, z1, 
                          i};
        segments->push_back(segment);
        // Update
        z0 = z1;
    }
    // Unwind this thing back to the station 
    for (int i = endLayer; i >= stationLayer; --i)
    {
        double z0 = interfaces[i];
        if (i == stationLayer)
        {
            z0 = stationDepth;
        }
        ::Segment segment{slownesses[i],
                          x0,
                          interfaces[i + 1],
                          x1,
                          z0,
                          i};
        segments->push_back(segment);
    }
    // Check convergence
    return ::checkRayConvergence(*segments, stationOffset, tolerance);
}

/// @brief Traces a symmetric ray path from a source, to a max
///        interface, then back up to a station at the same depth.
/// @param[in] interfaces    The (augmented) interfaces.  This indicates
///                          the depth of the top of each layer in meters
///                          which increases nadir.
/// @param[in] slownesses    The (agumented) slownesses in s/m in each
///                          layer.
/// @param[in] takeOffAngle  The take-off angle in degrees.
/// @param[in] sourceLayer   The layer containing the source.
/// @param[in] endLayer      The layer to which to trace.
/// @param[in] sourceDepth   The source depth in meters.
/// @param[in] stationOffset The source-receiver offset in meters.
/// @param[out] segments     The segments comprising the ray path.
/// @param[in] tolerance     A  hit is defined to be within this tolerance.
/// @result An indicator describing the ray convergence.
[[nodiscard]]
ReturnCode traceDown(const std::vector<double> &interfaces,
                     const std::vector<double> &slownesses,
                     const double takeOffAngle,
                     const int sourceLayer,
                     const int endLayer,
                     const double sourceDepth,
                     const double stationOffset,
                     std::vector<::Segment> *segments,
                     const bool allowCriticalRefractions,
                     const double tolerance)
{
    segments->clear();
    auto nLayers = static_cast<int> (interfaces.size()) - 1;
#ifndef NDEBUG
    assert(segments != nullptr);
    assert(interfaces.size() == slownesses.size());
    assert(sourceLayer >= 0 && sourceLayer <= endLayer);
    assert(endLayer < nLayers);
    assert(takeOffAngle >= 0 && takeOffAngle < 90);
    assert(stationOffset >= 0);
    assert(sourceDepth >= interfaces.at(sourceLayer) &&
           sourceDepth <  interfaces.at(sourceLayer + 1));
#endif
    if (segments->capacity() < 2*interfaces.size() + 1)
    {   
        segments->reserve(2*interfaces.size() + 1);
    } 
    double x0{0};
    double z0{sourceDepth};
    double currentAngle = takeOffAngle*(M_PI/180);
    int iStart =-1;
    for (int i = sourceLayer; i <= endLayer; ++i)
    {
        // Draw the ray
        double z1 = interfaces.at(i + 1);
        double x1 = x0 + (z1 - z0)*std::tan(currentAngle);
        ::Segment segment{slownesses[i],
                          x0, z0, x1, z1,
                          i};
        segments->push_back(segment);
        iStart = iStart + 1;
        // Turning?
        auto criticalAngle 
            = ::computeCriticalAngle(slownesses[i],
                                     slownesses[i + 1]);
        if (currentAngle > criticalAngle)
        {
            if (i < endLayer - 1)
            {
                segments->clear();
                return ReturnCode::RayTurnsTooEarly;
            }
            // Add the last segment by hand.  Compute the distance to the
            // half offset.
            auto dxHalf = 0.5*stationOffset - x1;
            // Now travel the mirror'ing critically refracted half
            // - i.e., 2*dxHalf.  Otherwise, we just bounce off this layer 
            if (allowCriticalRefractions && dxHalf > 0)
            {
                ::Segment segment{slownesses[i + 1],
                                  x1, z1, x1 + 2*dxHalf, z1,
                                  i + 1};
                segments->push_back(segment);
            }
            break; 
        }
        auto transmissionAngle
            = ::computeTransmissionAngle(currentAngle,
                                         slownesses[i],
                                         slownesses[i + 1]);
#ifndef NDEBUG
        assert(transmissionAngle >= 0);
#endif
        // Update
        x0 = x1;
        z0 = z1;
        currentAngle = transmissionAngle;
    }
    // Unwind this thing   
    x0 = segments->back().x1;
    for (int i = iStart; i >= 0; --i)
    {
        auto dx = segments->at(i).x1 - segments->at(i).x0;
        ::Segment segment{segments->at(i).slowness,
                          x0,
                          segments->at(i).z1,
                          x0 + dx,
                          segments->at(i).z0,
                          segments->at(i).layer};
        segments->push_back(segment);
        x0 = x0 + dx;
    }
    // Check convergence
    return ::checkRayConvergence(*segments, stationOffset, tolerance);
}

/// Traces from a direct ray upwards
ReturnCode traceDirect(const std::vector<double> &interfaces,
                       const std::vector<double> &slownesses,
                       const double takeOffAngle,
                       const int sourceLayer,
                       const double sourceDepth,
                       const int stationLayer,
                       const double stationOffset,
                       const double stationDepth,
                       std::vector<::Segment> *segments,
                       const double tolerance)
{
    segments->clear();
    auto nLayers = static_cast<int> (interfaces.size()) - 1; 
#ifndef NDEBUG
    assert(segments != nullptr);
    assert(takeOffAngle > 90 && takeOffAngle <= 180);
    assert(sourceLayer >= 0 && sourceLayer < nLayers);
    assert(stationOffset >= 0);
    assert(sourceDepth >= interfaces[sourceLayer] &&
           sourceDepth <  interfaces[sourceLayer + 1]);
    assert(stationDepth >= interfaces[stationLayer] &&
           stationDepth <  interfaces[stationLayer + 1]);
#endif
    // Flip this around
    if (sourceDepth < stationDepth)
    {
        auto returnCode = ::traceDirect(interfaces,
                                        slownesses,
                                        takeOffAngle,
                                        stationLayer,
                                        stationDepth,
                                        sourceLayer,
                                        stationOffset,
                                        sourceDepth,
                                        segments,
                                        tolerance);
        std::reverse(segments->begin(), segments->end()); 
        for (auto &segment : *segments)
        {
            segment.reverse();
        }
        return returnCode;
    }
#ifndef NDEBUG
    assert(stationLayer <= sourceLayer);
    assert(sourceDepth >= stationDepth);
#endif
    segments->reserve(sourceLayer - stationLayer + 1);
    // Simplify geometry
    double currentAngle = std::abs(180 - takeOffAngle)*(M_PI/180);
    if (takeOffAngle == 180){currentAngle = 0;}
    double x0{0};
    double z0{sourceDepth};
    for (int layer = sourceLayer; layer >= stationLayer; --layer) 
    {
        double z1 = interfaces[layer];
        if (layer == stationLayer){z1 = stationDepth;}
        double x1 = x0 + (z0 - z1)*std::tan(currentAngle);
        if (takeOffAngle == 180){x1 = x0;}
        ::Segment segment{slownesses[layer],
                          x0, z0, x1, z1,
                          layer};
        segments->push_back(segment);
        if (layer == stationLayer){break;}
        double transmissionAngle
            = ::computeTransmissionAngle(currentAngle,
                                         slownesses[layer],
                                         slownesses.at(layer - 1));
#ifndef NDEBUG
        assert(transmissionAngle >= 0); // Should be fine without velocity ivnersions
#endif
        if (transmissionAngle < 0)
        {
            std::cerr << "TraceDirect may have velocity inversion" << std::endl;
            segments->clear();
            return ReturnCode::RayTurnsTooEarly; 
        }
        //double criticalAngle = ::computeCriticalAngle(slownesses[layer],
        //                                              slownesses[layer - 1]);
        // Update
        x0 = x1;
        z0 = z1;
        currentAngle = transmissionAngle;
    }
    // Check convergence
    return checkRayConvergence(*segments, stationOffset, tolerance);
}

/// @brief Traces a ray down from the source, back up to the source depth,
///        then up to the station. 
ReturnCode traceDownThenUp(const std::vector<double> &augmentedInterfaces,
                           const std::vector<double> &augmentedSlownesses,
                           const double takeOffAngle,
                           const int sourceLayer,
                           const double sourceDepth,
                           const int stationLayer,
                           const double stationOffset,
                           const double stationDepth,
                           const int endLayer,
                           std::vector<::Segment> *segments,
                           const bool allowCriticalRefractions,
                           const double tolerance)
{
#ifndef NDEBUG
    assert(takeOffAngle >= 0 && takeOffAngle < 90);
    assert(segments != nullptr);
    assert(sourceDepth >= augmentedInterfaces[sourceLayer] &&
           sourceDepth <  augmentedInterfaces[sourceLayer + 1]);
    assert(stationDepth >= augmentedInterfaces[stationLayer] &&
           stationDepth <  augmentedInterfaces[stationLayer + 1]);
#endif
    segments->clear();
    if (std::abs(stationDepth - sourceDepth) < 1.e-10)
    {
        return ::traceDown(augmentedInterfaces,
                           augmentedSlownesses,
                           takeOffAngle,
                           sourceLayer,
                           endLayer,
                           sourceDepth,
                           stationOffset,
                           segments,
                           allowCriticalRefractions,
                           tolerance);
    }
    // Deal with the station below source by interchanging the values
    // of the source and station and then reversing the result
    if (sourceDepth < stationDepth)
    {
        auto returnCode = traceDownThenUp(augmentedInterfaces,
                                          augmentedSlownesses,
                                          takeOffAngle,
                                          stationLayer, //const int sourceLayer,
                                          stationDepth, //const double sourceDepth,
                                          sourceLayer,  //const int stationLayer,
                                          stationOffset,
                                          sourceDepth,  // const double stationDepth,
                                          endLayer,
                                          segments,
                                          allowCriticalRefractions,
                                          tolerance);
        ::reverseSegments(segments);
        return returnCode;
    }
    // Trace out the top segments
    auto takeOffAngleUp = 180 - takeOffAngle;
    std::vector<::Segment> upSegments;
    auto result = ::traceDirect(augmentedInterfaces,
                                augmentedSlownesses,
                                takeOffAngleUp,
                                sourceLayer,
                                sourceDepth,
                                stationLayer,
                                stationOffset,
                                stationDepth,
                                &upSegments,
                                tolerance);
    if (result != ReturnCode::Hit &&
        result != ReturnCode::UnderShot &&
        result != ReturnCode::OverShot)
    {
        return result;
    }
    // Trace the bottom segments from source to `source'
#ifndef NDEBUG
    assert(upSegments.front().x0 == 0);
#endif
    auto offsetUp = upSegments.back().x1;
    std::vector<::Segment> downSegments;
    result = ::traceDownThenUp(augmentedInterfaces,
                               augmentedSlownesses,
                               takeOffAngle,
                               sourceLayer,
                               sourceDepth,
                               sourceLayer, 
                               std::max(0.0, stationOffset - offsetUp),
                               sourceDepth,
                               endLayer,
                               &downSegments,
                               allowCriticalRefractions,
                               tolerance);
    if (result != ReturnCode::Hit &&
        result != ReturnCode::UnderShot &&
        result != ReturnCode::OverShot)
    {
        return result;
    }
    auto downOffset = downSegments.back().x1;
    // Put the ray path together 
    segments->resize(downSegments.size() + upSegments.size());
    std::copy(downSegments.begin(), downSegments.end(),
              segments->begin());
    int j = static_cast<int> (downSegments.size());
    for (auto &upSegment : upSegments)
    {
        upSegment.x0 = upSegment.x0 + downOffset;
        upSegment.x1 = upSegment.x1 + downOffset;
        segments->at(j) = std::move(upSegment);
        j = j + 1;
    }
    return ::checkRayConvergence(*segments, stationOffset, tolerance);
}

/// @brief Converts the 1D velocity stack to a slowness stack. 
[[nodiscard]]
std::vector<double> toSlownessVector(const std::vector<double> &velocities)
{
    std::vector<double> slownesses(velocities.size());
    std::transform(velocities.begin(), velocities.end(), slownesses.begin(),
                   [&](const double velocity)
                   {
#ifndef NDEBUG
                       assert(velocity >= 0);
#endif
                       return 1./velocity;
                   });
    return slownesses;
}

/// @brief Modifies the velocity vector for sfae use with the above routines.
[[nodiscard]] 
std::vector<double> augmentVelocityVector(const std::vector<double> &x)
{
    auto y = x;
    if (!y.empty()){y.push_back(x.back());}
    return y;
}

/// @brief Modifies the interfaces vector for safe use with the above routines.
[[nodiscard]]
std::vector<double> augmentInterfacesVector(const std::vector<double> &x)
{
    auto y = x;
    y.push_back(std::numeric_limits<double>::max());
    return y;
}

/// Traces a direct wave when the source and receiver are in the smae layer.
[[nodiscard]]
ReturnCode
    traceDirectSameLayer(const std::vector<double> &interfaces,
                         const std::vector<double> &slownesses,
                         const int sourceLayer,
                         const double sourceDepth,
                         const int stationLayer,
                         const double stationOffset,
                         const double stationDepth,
                         EikonalXX::Ray::Path2D *path)
{
    path->clear();
#ifndef NDEBUG
    assert(sourceLayer == stationLayer);
    assert(!slownesses.empty());
    assert(slownesses.size() == interfaces.size());
    assert(sourceLayer < interfaces.size() - 1);
    assert(path != nullptr);
#endif
    EikonalXX::Ray::Point2D startPoint{0, sourceDepth};
    EikonalXX::Ray::Point2D endPoint{stationOffset, stationDepth};
    EikonalXX::Ray::Segment2D segment;
    segment.setStartAndEndPoint(std::pair {startPoint, endPoint});
    segment.setSlowness(slownesses.at(sourceLayer));
    segment.setVelocityModelCellIndex(sourceLayer);
    path->open();
    path->append(segment);
    path->close();
    return ReturnCode::Hit;
}

/// @brief Traces a direct ray in a whole space.
ReturnCode
    traceWholeSpace(const double slowness,
                    const double sourceDepth,
                    const double stationOffset,
                    const double stationDepth,
                    EikonalXX::Ray::Path2D *path)
{
#ifndef NDEBUG
    assert(slowness >= 0);
    assert(path != nullptr);
#endif
    const std::vector<double> interfaces{std::numeric_limits<double>::min(),
                                         std::numeric_limits<double>::max()};
    const std::vector<double> slownesses{slowness, slowness};
    return ::traceDirectSameLayer(interfaces, slownesses,
                                  0, sourceDepth,
                                  0, stationOffset, stationDepth,
                                  path);
}

[[nodiscard]] std::vector<EikonalXX::Ray::Path2D>
    computeVerticalRayPaths(const std::vector<double> &augmentedInterfaces,
                            const std::vector<double> &augmentedSlownesses,
                            const int sourceLayer,
                            const double sourceDepth,
                            const int stationLayer,
                            const double stationDepth,
                            const double stationOffset,
                            const double rayHitTolerance = 1)
{
    std::vector<EikonalXX::Ray::Path2D> rayPaths;
    // Trace the direct ray path upwards 
    if (sourceDepth >= stationDepth)
    {
        constexpr double upTakeOffAngle{180};
        std::vector<::Segment> segments;
#ifndef NDEBUG
        auto returnCode =
#endif
        ::traceDirect(augmentedInterfaces,
                      augmentedSlownesses,
                      upTakeOffAngle,
                      sourceLayer,
                      sourceDepth,
                      stationLayer,
                      stationOffset,
                      stationDepth,
                      &segments,
                      rayHitTolerance);
#ifndef NDEBUG 
        assert(returnCode == ReturnCode::Hit ||
               returnCode == ReturnCode::UnderShot);
#endif
        if (!segments.empty()){rayPaths.push_back(::toRayPath(segments));}
    }
    else // Same layer but station below
    {
        if (sourceLayer == stationLayer)
        {
            EikonalXX::Ray::Path2D path;     
#ifndef NDEBUG
            auto returnCode =
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
            assert(returnCode == ReturnCode::Hit);
#endif
            rayPaths.push_back(path);
        }
    }
    // Loop through stack and bounce rays
    auto nLayers = static_cast<int> (augmentedInterfaces.size()) - 1;
    for (int layer = sourceLayer; layer < nLayers - 1; ++layer)
    {
        std::vector<::Segment> segments;
#ifndef NDEBUG
        auto returnCode = 
#endif
        ::traceVerticalReflectionDown(augmentedInterfaces,
                                      augmentedSlownesses,
                                      sourceLayer,
                                      layer,
                                      sourceDepth,
                                      stationLayer,
                                      stationDepth,
                                      stationOffset,
                                      &segments,
                                      rayHitTolerance);
#ifndef NDEBUG
        assert(returnCode == ReturnCode::Hit ||
               returnCode == ReturnCode::UnderShot);
#endif
        if (!segments.empty()){rayPaths.push_back(::toRayPath(segments));}
    }
    std::sort(rayPaths.begin(), rayPaths.end(),
              [](const EikonalXX::Ray::Path2D &lhs,
                 const EikonalXX::Ray::Path2D &rhs)
              {
                  return lhs.getTravelTime() < rhs.getTravelTime();
              });
    return rayPaths;
}

/// @brief Shoots a ray with the given take-off angle.
/// @note This cannot handle velocity inversions. 
[[nodiscard]] [[maybe_unused]]
std::vector<EikonalXX::Ray::Path2D> 
    shoot(const double takeOffAngle,
          const std::vector<double> &augmentedInterfaces,
          const std::vector<double> &augmentedSlownesses,
          const int sourceLayer,
          const double sourceDepth,
          const int stationLayer,
          const double stationOffset,
          const double stationDepth,
          const bool allowCriticalRefractions,
          const double rayHitTolerance, 
          const bool keepOnlyHits,
          const int lastLayer)
{
    std::vector<EikonalXX::Ray::Path2D> result;
    // Edge case - straight down or up
    if (takeOffAngle < 1.e-7 || std::abs(180 - takeOffAngle) < 1.e-7)
    {
        if (keepOnlyHits && stationOffset < rayHitTolerance)
        {
            std::vector<EikonalXX::Ray::Path2D> temporaryRays;  
            try
            {
                if (takeOffAngle < 1.e-7)
                {
                    temporaryRays
                       = ::computeVerticalRayPaths(augmentedInterfaces,
                                                   augmentedSlownesses,
                                                   sourceLayer,
                                                   sourceDepth,
                                                   stationLayer,
                                                   stationDepth,
                                                   stationOffset,
                                                   rayHitTolerance);
                }
                else
                {
                    std::vector<::Segment> segments;
                    ::traceDirect(augmentedInterfaces,
                                  augmentedSlownesses,
                                  180,
                                  sourceLayer,
                                  sourceDepth,
                                  stationLayer,
                                  stationOffset,
                                  stationDepth,
                                  &segments,
                                  rayHitTolerance); 
                    temporaryRays.push_back(::toRayPath(segments));
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "shoot error in computeVerticalRayPaths 1: "
                          << e.what() << std::endl;
            }
            result.reserve(temporaryRays.size());
            for (auto &ray : temporaryRays)
            {
                auto nSegments = ray.size();
                if (nSegments > 0)
                {
                    auto offset
                         = ray.at(nSegments - 1).getEndPoint().getPositionInX();
                    if (std::abs(offset - stationOffset) < rayHitTolerance)
                    {
                        result.push_back(std::move(ray));
                    }
                }
            }
        }
        else
        {
            try
            {
                if (takeOffAngle < 1.e-7)
                {
                    result = ::computeVerticalRayPaths(augmentedInterfaces,
                                                       augmentedSlownesses,
                                                       sourceLayer,
                                                       sourceDepth,
                                                       stationLayer,
                                                       stationDepth,
                                                       stationOffset,
                                                       rayHitTolerance);
                }
                else
                {
                    std::vector<::Segment> segments;
                    ::traceDirect(augmentedInterfaces,
                                  augmentedSlownesses,
                                  180,
                                  sourceLayer,
                                  sourceDepth,
                                  stationLayer,
                                  stationOffset,
                                  stationDepth,
                                  &segments,
                                  rayHitTolerance); 
                    result.push_back(::toRayPath(segments));
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "shoot error in computeVerticalRayPaths 2: "
                          << e.what() << std::endl;
            }
        } 
        return result;
    }
    // 90 degree take-off angles won't converge
    if (std::abs(takeOffAngle - 90) < 1.e-7){return result;}
    // Up-going rays
    if (takeOffAngle > 90)
    {
        // Not bouncing off above layers so just quit now
        if (stationDepth > sourceDepth){return result;}
        std::vector<::Segment> segments;
        try
        {
            auto returnCode = ::traceDirect(augmentedInterfaces,
                                            augmentedSlownesses,
                                            takeOffAngle,
                                            sourceLayer,
                                            sourceDepth,
                                            stationLayer,
                                            stationOffset,
                                            stationDepth,
                                            &segments,
                                            rayHitTolerance); 
            if (returnCode == ReturnCode::Hit ||
                (!keepOnlyHits && returnCode == ReturnCode::UnderShot) ||
                (!keepOnlyHits && returnCode == ReturnCode::OverShot))
            {
                auto rayPath = ::toRayPath(segments);
                //auto descriptor = "Hit";
                //if (returnCode == ReturnCode::UnderShot)
                //{
                //    descriptor = "Under-shot";
                //}
                //else if (returnCode == ReturnCode::OverShot)
                //{
                //    descriptor = "Over-shot";
                //}
                //std::cout << descriptor << " offset, takeoff angle, travel time: " << stationOffset << "," << takeOffAngle << "," << rayPath.getTravelTime() << std::endl;
                result.push_back(rayPath);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Problem in traceDirect called from shoot: "
                      << e.what() << std::endl;
        }
    }
    else // Down-going rays
    {
        // Trace down through the velocity model stack and tabulate the
        // ray paths
        auto nLayers = static_cast<int> (augmentedInterfaces.size()) - 1;
        for (int endLayer = sourceLayer;
             endLayer < std::min(lastLayer, nLayers - 1);
             ++endLayer)
        {
            try
            {
                std::vector<::Segment> segments;
                auto returnCode = ::traceDownThenUp(augmentedInterfaces,
                                                    augmentedSlownesses,
                                                    takeOffAngle,
                                                    sourceLayer,
                                                    sourceDepth,
                                                    stationLayer,
                                                    stationOffset,
                                                    stationDepth,
                                                    endLayer,
                                                    &segments,
                                                    allowCriticalRefractions,
                                                    rayHitTolerance);
                if (returnCode == ReturnCode::Hit ||
                    (!keepOnlyHits && returnCode == ReturnCode::UnderShot) ||
                    (!keepOnlyHits && returnCode == ReturnCode::OverShot))
                {
                    auto rayPath = ::toRayPath(segments);
                    //std::cout << std::endl;
                    //for (const auto &segment : path)
                    //{
                    //    std::cout << std::setprecision(12)
                    //              << segment.getStartPoint().getPositionInX() << ","
                    //              << segment.getStartPoint().getPositionInZ() << ","
                    //              << segment.getEndPoint().getPositionInX() << ","
                    //              << segment.getEndPoint().getPositionInZ() << ","
                    //              << 1./segment.getSlowness() << "," 
                    //              << 1./mSlownessModel[segment.getVelocityModelCellIndex()] << std::endl;
                    //}
                    //std::cout << "Offset, takeoff angle, travel time: " << stationOffset << "," << takeOffAngle << "," << rayPath.getTravelTime() << std::endl;
                    result.push_back(rayPath);
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error in traceDownThenUp called from shoot: "
                          << e.what() << std::endl;
            }
        } // Loop on stack of layers
    } // End check on upgoing vs downgoing
    // Sort the rays in increasing order based on the travel times
    if (!result.empty())
    {
        std::sort(result.begin(), result.end(), 
                  [](const EikonalXX::Ray::Path2D &lhs,
                     const EikonalXX::Ray::Path2D &rhs)
                  {
                     return lhs.getTravelTime() < rhs.getTravelTime();
                  });
    }
    return result;
}

}

#endif
