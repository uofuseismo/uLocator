#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
#ifndef NDEBUG
#include <cassert>
#endif
#include "eikonalxx/ray/path2d.hpp"
#include "eikonalxx/ray/segment2d.hpp"
#include "eikonalxx/ray/point2d.hpp"
#include "length.hpp"

using namespace EikonalXX::Ray;

namespace
{
void checkSegment(const Segment2D &segment)
{
    if (!segment.haveStartAndEndPoint())
    {
        throw std::invalid_argument("Segment start/end point not set");
    }
    if (!segment.haveVelocity())
    {
        throw std::invalid_argument("Segment velocity not set");
    }
}
}

class Path2D::Path2DImpl
{
public:
    void reverse()
    {
        if (mSegments.empty()){return;}
        std::reverse(mSegments.begin(), mSegments.end());
        for (auto &segment : mSegments)
        {
            segment.reverse();
        }
#ifndef NDEBUG
        for (int i = 0; i < static_cast<int> (mSegments.size()) - 1; ++i)
        {
            auto point0 = mSegments[i].getEndPoint(); 
            auto point1 = mSegments[i + 1].getStartPoint();
            assert(std::abs(point0.getPositionInX() - point1.getPositionInX())
                   < 1.e-8);
            assert(std::abs(point0.getPositionInZ() - point1.getPositionInZ())
                   < 1.e-8);
        }
#endif
    }
    void integrate()
    {
        double length = 0;
        double travelTime = 0;
        for (const auto &segment : mSegments)
        {
            travelTime = travelTime + segment.getTravelTime();
            length = length + segment.getLength();
        }
        mTravelTime = travelTime;
        mLength = length;
    }
    std::vector<Segment2D> mSegments;
    std::vector<Segment2D> mSegmentsConstruction;
    double mTravelTime{0};
    double mLength{0};
    bool mOpened{false};
};

/// Constructor
Path2D::Path2D() :
    pImpl(std::make_unique<Path2DImpl> ())
{
}

/// Copy constructor
Path2D::Path2D(const Path2D &path)
{
    *this = path;
}

/// Move constructor
Path2D::Path2D(Path2D &&path) noexcept
{
    *this = std::move(path);
}

/// Reset class
void Path2D::clear() noexcept
{
    pImpl = std::make_unique<Path2DImpl> ();
}

/// Destructor
Path2D::~Path2D() = default;

/// Copy assignment
Path2D& Path2D::operator=(const Path2D &path)
{
    if (&path == this){return *this;}
    pImpl = std::make_unique<Path2DImpl> (*path.pImpl);
    return *this;
}

/// Move constructor
Path2D& Path2D::operator=(Path2D &&path) noexcept
{
    if (&path == this){return *this;}
    pImpl = std::move(path.pImpl);
    return *this;
}

/// Open the ray path
void Path2D::open()
{
    pImpl->mSegmentsConstruction.clear();
    pImpl->mSegmentsConstruction.reserve(
        std::max(128, static_cast<int> (pImpl->mSegments.size())));
    pImpl->mOpened = true;
}

bool Path2D::isOpen() const noexcept
{
    return pImpl->mOpened;
}

/// Append to the ray path
void Path2D::append(const Segment2D &segment)
{
    auto temporarySegment = segment;
    append(std::move(temporarySegment));
}

void Path2D::append(Segment2D &&segment)
{
    ::checkSegment(segment);
    if (pImpl->mSegmentsConstruction.empty())
    {
        pImpl->mSegmentsConstruction.push_back(std::move(segment));
    }
    else
    {
        auto point0 = pImpl->mSegmentsConstruction.back().getEndPoint();
        auto point1 = segment.getStartPoint();
        auto distance = ::computeLength(point0, point1);
        if (distance > std::numeric_limits<double>::epsilon()*100)
        {
            throw std::invalid_argument(
               "Segment does not start at last segments's end point"); 
        }
        pImpl->mSegmentsConstruction.push_back(std::move(segment));
    }
}

void Path2D::set(const std::vector<Segment2D> &segments)
{
    auto temporarySegments = segments;
    set(std::move(temporarySegments));
}

void Path2D::set(std::vector<Segment2D> &&segments)
{
    for (int i = 0; i < static_cast<int> (segments.size()); ++i)
    {
        ::checkSegment(segments[i]);
        if (i > 0)
        {
            auto distance = ::computeLength(segments[i-1].getEndPoint(),
                                            segments[i].getStartPoint());
            if (distance > std::numeric_limits<double>::epsilon()*100)
            {
                throw std::invalid_argument("Segment " + std::to_string(i)
                                          + " does no start at last segments "
                                          + std::to_string(i - 1)
                                          + " end point");
            }
        }
    }
    pImpl->mSegments = std::move(segments);
    pImpl->mSegmentsConstruction.clear();
    pImpl->mOpened = false;             
    // sum along the path
    pImpl->integrate();
}

/// Close the ray path
void Path2D::close()
{
    if (!isOpen()){throw std::runtime_error("Path not open");}
    // Do some tap dancing and release memory
    pImpl->mSegments = std::move(pImpl->mSegmentsConstruction);
    pImpl->mSegmentsConstruction = std::vector<Segment2D> ();
    pImpl->mOpened = false;
    // Sum along the path
    pImpl->integrate();
}

size_t Path2D::size() const noexcept
{
    return static_cast<int> (pImpl->mSegments.size());
}

bool Path2D::empty() const noexcept
{
    return pImpl->mSegments.empty();
}

/// Take-off angle
double Path2D::getTakeOffAngle() const
{
    if (empty()){throw std::runtime_error("No segments in path");}
    auto nSegments = pImpl->mSegments.size();
    double dx = 0;
    double dz = 0;
    bool foundOne = false;
    // TODO i might want to average the first few
    double oneTakeOffAngle = 0;
    double averageTakeOffAngle = 0;
    int nAverage = 1;
    int iAverage = 0;
    for (int i = 0; i < nSegments; ++i)
    {
        if (pImpl->mSegments.at(i).getLength() > 0)
        {
            auto point0 = pImpl->mSegments.at(i).getStartPoint();
            auto point1 = pImpl->mSegments.at(i).getEndPoint();
            dx = point1.getPositionInX() - point0.getPositionInX();
            dz = point1.getPositionInZ() - point0.getPositionInZ();
            auto angle = std::atan2(std::abs(dx), dz)*(180/M_PI);
            if (!foundOne)
            {
                oneTakeOffAngle = angle;
            }
            averageTakeOffAngle = averageTakeOffAngle + angle;
            iAverage = iAverage + 1;
            foundOne = true;
            if (iAverage == nAverage){break;}
         }
    }
    if (!foundOne)
    {
        std::cerr << "There may exist a zero-length segment at start"
                  << std::endl;
    }
//std::cout << std::setprecision(16) << dx << "," << dz << std::endl;
    // +z is down so this is like atan2's +x
    // +x is is like atan2's +y.  however, i don't care about the sign
    // function is : atan2(y, x) -> atan2(dx, dz)
    if (iAverage == nAverage)
    {
        return averageTakeOffAngle/nAverage;
    }
    return oneTakeOffAngle;
}

/// Travel time
double Path2D::getTravelTime() const
{
    return pImpl->mTravelTime;
}

/// Length
double Path2D::getLength() const
{
    return pImpl->mLength;
}

/// Iterators
Path2D::iterator Path2D::begin()
{
    return pImpl->mSegments.begin();
}

Path2D::const_iterator Path2D::begin() const noexcept
{
    return pImpl->mSegments.begin();
}

Path2D::const_iterator Path2D::cbegin() const noexcept
{
    return pImpl->mSegments.cbegin();
}

Path2D::iterator Path2D::end()
{
    return pImpl->mSegments.end();
}

Path2D::const_iterator Path2D::end() const noexcept
{
    return pImpl->mSegments.cend();
}

Path2D::const_iterator Path2D::cend() const noexcept
{
    return pImpl->mSegments.cend();
}

Segment2D& Path2D::at(const size_t index)
{
    return pImpl->mSegments.at(index);
}

const Segment2D& Path2D::at(const size_t index) const
{
    return pImpl->mSegments.at(index);
}

void Path2D::reverse()
{
    if (isOpen()){throw std::runtime_error("Construction is open");}
    pImpl->reverse();
}
