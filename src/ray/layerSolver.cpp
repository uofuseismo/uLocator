#include <iostream>
#include <iomanip>
#include <functional>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#ifndef NDEBUG
#include <cassert>
#endif
#include "layerSolver.hpp"
#include "segment2d.hpp"
#include "path2d.hpp"
#include "point2d.hpp"
#include "layerTracer.hpp"
#include "optimize.hpp"

/// 1.e-4 changes about 1 m for every 500 km
#define MIN_DOWNGOING_ANGLE 0.0
//1.e-4
#define MAX_DOWNGOING_ANGLE 89.9999
#define MIN_UPGOING_ANGLE 90.0001
#define MAX_UPGOING_ANGLE 180


using namespace EikonalXX::Ray;

class LayerSolver::LayerSolverImpl
{
public:
    /// For when the source and station are in the same layer this computes
    /// the direct ray. 
    [[nodiscard]] Path2D computeDirectRaySameLayer() const
    {
#ifndef NDEBUG
        assert(mSourceLayer == mStationLayer);
#endif
        Path2D path;
#ifndef NDEBUG
        auto returnCode =
#endif
        ::traceDirectSameLayer(mAugmentedInterfaces,
                               mAugmentedSlownesses,
                               mSourceLayer,
                               mSourceDepth,
                               mStationLayer,
                               mStationOffset,
                               mStationDepth,
                               &path);
#ifndef NDEBUG 
        assert(returnCode == ReturnCode::Hit);
#endif
        return path;
    }
    /// For the whole or half space problem
    [[nodiscard]] Path2D computeWholeSpace() const
    {
#ifndef NDEBUG
        assert(mSourceLayer == 0);
        assert(mStationLayer == mSourceLayer);
#endif
        Path2D rayPath;
#ifndef NDEBUG
        auto returnCode =
#endif
        ::traceWholeSpace(mSlownessModel[mSourceLayer],
                          mSourceDepth,
                          mStationOffset,
                          mStationDepth,
                          &rayPath);
#ifndef NDEBUG
        assert(returnCode == ReturnCode::Hit);
#endif
        return rayPath;
    }
    /// For when the source and station are vertically aligned this computes
    /// the direct ray (assuming the source is below the station) and the
    /// reflected ray paths.
    [[nodiscard]] std::vector<Path2D> computeVerticalRayPaths() const
    {
        std::vector<Path2D> rayPaths;
        if (mSourceLayer == mStationLayer)
        {
            rayPaths.push_back(computeDirectRaySameLayer());
        }
        else
        {
            constexpr double upTakeOffAngle{180};
            std::vector<::Segment> segments;
#ifndef NDEBUG
            auto returnCode =
#endif
            ::traceDirect(mAugmentedInterfaces,
                          mAugmentedSlownesses,
                          upTakeOffAngle,
                          mSourceLayer,
                          mSourceDepth,
                          mStationLayer,
                          mStationOffset,
                          mStationDepth,
                          &segments,
                          mRayHitTolerance);
#ifndef NDEBUG 
            assert(returnCode == ReturnCode::Hit ||
                  returnCode == ReturnCode::UnderShot);
#endif
            if (!segments.empty()){rayPaths.push_back(::toRayPath(segments));}
        }
        // Loop through stack and bounce rays
        auto nLayers = static_cast<int> (mInterfaces.size());
        for (int layer = std::max(mSourceLayer, mStationLayer);
             layer < nLayers - 1;
             ++layer)
        {
            std::vector<::Segment> segments;
#ifndef NDEBUG
            auto returnCode = 
#endif
                 ::traceVerticalReflectionDown(mAugmentedInterfaces,
                                               mAugmentedSlownesses,
                                               mSourceLayer,
                                               layer,
                                               mSourceDepth,
                                               mStationLayer,
                                               mStationDepth,
                                               mStationOffset,
                                               &segments,
                                               mRayHitTolerance);
#ifndef NDEBUG
            assert(returnCode == ReturnCode::Hit ||
                   returnCode == ReturnCode::UnderShot);
#endif
            if (!segments.empty()){rayPaths.push_back(::toRayPath(segments));}
        }
        std::sort(rayPaths.begin(), rayPaths.end(),
                  [](const Path2D &lhs, const Path2D &rhs)
                  {
                      return lhs.getTravelTime() < rhs.getTravelTime();
                  });
        return rayPaths;
    }
    /// @brief Shoots a ray with the given take-off angle.
    /// @note This cannot handle velocity inversions. 
    std::vector<Path2D> shoot(const double takeOffAngle)
    {
        auto lastLayer = static_cast<int> (mAugmentedInterfaces.size()) - 2;
        constexpr bool allowCriticalRefractions{true};
        return ::shoot(takeOffAngle,
                       mAugmentedInterfaces,
                       mAugmentedSlownesses,
                       mSourceLayer,
                       mSourceDepth,
                       mStationLayer,
                       mStationOffset,
                       mStationDepth,
                       allowCriticalRefractions,
                       mRayHitTolerance,
                       mKeepOnlyHits,
                       lastLayer);
    }
    std::vector<Path2D> mRayPaths;
    std::vector<double> mSlownessModel;
    std::vector<double> mInterfaces;
    std::vector<double> mAugmentedSlownesses;
    std::vector<double> mAugmentedInterfaces;
    std::function<std::vector<Path2D> (const double )> 
        mShootRay{std::bind(&LayerSolverImpl::shoot,
                            this,
                            std::placeholders::_1)};
    double mRayHitTolerance{150};
    double mMaximumOffset{500000};
    double mStationDepth{0};
    double mSourceDepth{0};
    double mStationOffset{0};
    int mSourceLayer{0};
    int mStationLayer{0};
    bool mHaveStationOffsetAndDepth{false};
    bool mHaveSourceDepth{false};
    bool mHaveRayPaths{false};
    bool mAllowCriticalRefractions{true};
    bool mKeepOnlyHits{true};
};

/// Constructor
LayerSolver::LayerSolver() :
    pImpl(std::make_unique<LayerSolverImpl> ())
{
}

/// Move constructor
LayerSolver::LayerSolver(LayerSolver &&solver) noexcept
{
    *this = std::move(solver);
}

/// Move assignment
LayerSolver& LayerSolver::operator=(LayerSolver &&solver) noexcept
{
    if (&solver == this){return *this;}
    pImpl = std::move(solver.pImpl);
    return *this;
}

/// Reset class
void LayerSolver::clear() noexcept
{
    pImpl = std::make_unique<LayerSolverImpl> ();
}

/// Destructor
LayerSolver::~LayerSolver() = default;

/// Set the velocity model
void LayerSolver::setVelocityModel(const std::vector<double> &interfaces,
                                   const std::vector<double> &velocityModel)
{
    if (velocityModel.empty())
    {
        throw std::invalid_argument("Velocity model is empty");
    }
    if (interfaces.size() != velocityModel.size())
    {
        throw std::invalid_argument(
            "Velocity model size must equal interfaces size"
        );
    }
    for (const auto &velocity : velocityModel)
    {
        if (velocity <= 0)
        {
            throw std::invalid_argument("All velocities must be positive");
        }
    }
    for (int i = 0; i < static_cast<int> (velocityModel.size()) - 1; ++i)
    {
        if (velocityModel[i + 1] <= velocityModel[i])
        {
            throw std::invalid_argument("Velocity inversions not yet done");
        }
    }
    if (!std::is_sorted(interfaces.begin(), interfaces.end(),
                        [=](const double lhs, const double rhs)
                        {
                            return lhs < rhs;
                        }))
    {
        throw std::invalid_argument("Interfaces must increase with depth");
    }
    pImpl->mSlownessModel = ::toSlownessVector(velocityModel);
    pImpl->mInterfaces = interfaces;
    pImpl->mAugmentedSlownesses
         = ::toSlownessVector(::augmentVelocityVector(velocityModel));
    pImpl->mAugmentedInterfaces
         = ::augmentInterfacesVector(pImpl->mInterfaces);
}

bool LayerSolver::haveVelocityModel() const noexcept
{
    return !pImpl->mSlownessModel.empty();
}

/// Source depth
void LayerSolver::setSourceDepth(const double depth)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (depth < pImpl->mInterfaces[0])
    {
        throw std::invalid_argument("Source is in the air");
    }
    pImpl->mSourceDepth = depth; 
    pImpl->mSourceLayer = ::getLayer(pImpl->mSourceDepth,
                                     pImpl->mAugmentedInterfaces,
                                     true);
    pImpl->mHaveSourceDepth = true;
}

double LayerSolver::getSourceDepth() const
{
    if (!haveSourceDepth()){throw std::runtime_error("Source depth not set");}
    return pImpl->mSourceDepth;
}

bool LayerSolver::haveSourceDepth() const noexcept
{
    return pImpl->mHaveSourceDepth;
}

/// Station
void LayerSolver::setStationOffset(const double offset)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (offset < 0){throw std::invalid_argument("Offset must be non-negative");}
    pImpl->mStationOffset = offset;
    pImpl->mStationDepth = pImpl->mInterfaces.at(0);
    pImpl->mStationLayer = ::getLayer(pImpl->mStationDepth,
                                      pImpl->mAugmentedInterfaces,
                                      true);
    pImpl->mHaveStationOffsetAndDepth = true;
}

void LayerSolver::setStationOffsetAndDepth(const double offset,
                                           const double depth)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (depth < pImpl->mInterfaces[0])
    {   
        throw std::invalid_argument("Station is in the air");
    }
    if (offset < 0){throw std::invalid_argument("Offset must be non-negative");}
    pImpl->mStationOffset = offset;
    pImpl->mStationDepth = depth; 
    pImpl->mStationLayer = ::getLayer(pImpl->mStationDepth,
                                      pImpl->mAugmentedInterfaces,
                                      true);
    pImpl->mHaveStationOffsetAndDepth = true;
#ifndef NDEBUG
    assert(pImpl->mStationLayer >= 0 &&
           pImpl->mStationLayer < static_cast<int> (pImpl->mInterfaces.size()));
#endif
}

bool LayerSolver::haveStationOffsetAndDepth() const noexcept
{
    return pImpl->mHaveStationOffsetAndDepth;
}

/// Get the ray paths
std::vector<Path2D> LayerSolver::getRayPaths() const
{
    if (!haveRayPaths()){throw std::runtime_error("Ray paths not computed");}
    return pImpl->mRayPaths;
}

const std::vector<Path2D> &LayerSolver::getRayPathsReference() const
{
    return *&pImpl->mRayPaths;
}

/// Have ray path?
bool LayerSolver::haveRayPaths() const noexcept
{
    return pImpl->mHaveRayPaths;
}

/// Shoot
std::vector<Path2D> LayerSolver::shoot(const double takeOffAngle,
                                       const bool keepOnlyHits)
{
    if (!haveVelocityModel())
    {   
        throw std::runtime_error("Velocity model not set");
    }   
    if (!haveSourceDepth()){throw std::runtime_error("Source depth not set");}
    if (!haveStationOffsetAndDepth())
    {   
        throw std::runtime_error("Station offset not set");
    }
    if (takeOffAngle < 0 || takeOffAngle > 180)
    {
        throw std::invalid_argument("Take-off angle must be in range [0,180]");
    }
    pImpl->mKeepOnlyHits = keepOnlyHits;
    return pImpl->shoot(takeOffAngle);
}

/// Solve
void LayerSolver::solve(const bool doReflections)
{
    if (!haveVelocityModel())
    {
        throw std::runtime_error("Velocity model not set");
    }
    if (!haveSourceDepth()){throw std::runtime_error("Source depth not set");}
    if (!haveStationOffsetAndDepth())
    {
        throw std::runtime_error("Station offset not set");
    }
    Path2D rayPath;
    pImpl->mHaveRayPaths = false;
    pImpl->mRayPaths.clear();
    // Special case of just one layer
    if (pImpl->mSlownessModel.size() == 1)
    {
        pImpl->mRayPaths.push_back(pImpl->computeWholeSpace());
        pImpl->mHaveRayPaths = true;
        return;
    }
    // Special case station directly above/below station
    if (pImpl->mStationOffset < std::numeric_limits<double>::epsilon())
    {
        pImpl->mRayPaths = pImpl->computeVerticalRayPaths();
        pImpl->mHaveRayPaths = true;
        return;
    }
    // General case where we have to optimize
    auto directRayPath = ::optimizeDirect(pImpl->mAugmentedInterfaces,
                                          pImpl->mAugmentedSlownesses,
                                          pImpl->mSourceLayer,
                                          pImpl->mSourceDepth,
                                          pImpl->mStationLayer,
                                          pImpl->mStationDepth,
                                          pImpl->mStationOffset,
                                          pImpl->mRayHitTolerance);
    if (!directRayPath.empty())
    {
        pImpl->mRayPaths.insert(pImpl->mRayPaths.end(),
                                directRayPath.begin(),
                                directRayPath.end());
    }
    // Optimize the critically refracted rays
    auto refractedRayPaths = optimizeCriticallyRefracted(
        pImpl->mAugmentedInterfaces,
        pImpl->mAugmentedSlownesses,
        pImpl->mSourceLayer,
        pImpl->mSourceDepth,
        pImpl->mStationLayer,
        pImpl->mStationDepth,
        pImpl->mStationOffset);
    if (!refractedRayPaths.empty())
    {
        pImpl->mRayPaths.insert(pImpl->mRayPaths.end(),
                                refractedRayPaths.begin(),
                                refractedRayPaths.end());
    }
    // Optimize reflections
    if (doReflections)
    {
        auto reflectedRayPaths
            = optimizeFirstArrivingDownGoingReflections(
                 pImpl->mAugmentedInterfaces,
                 pImpl->mAugmentedSlownesses,
                 pImpl->mSourceLayer,
                 pImpl->mSourceDepth,
                 pImpl->mStationLayer,
                 pImpl->mStationDepth,
                 pImpl->mStationOffset,
                 pImpl->mRayHitTolerance);
        if (!reflectedRayPaths.empty())
        {   
            pImpl->mRayPaths.insert(pImpl->mRayPaths.end(),
                                    reflectedRayPaths.begin(),
                                    reflectedRayPaths.end());
        }   
    }
    // Sort these based on travel time
    std::sort(pImpl->mRayPaths.begin(), pImpl->mRayPaths.end(), 
              [&](const Path2D &lhs, const Path2D &rhs)
              {
                  return lhs.getTravelTime() < rhs.getTravelTime();
              });
    pImpl->mHaveRayPaths = true;
}
