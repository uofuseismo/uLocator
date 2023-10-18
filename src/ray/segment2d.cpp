#include <string>
#include <cmath>
#include "segment2d.hpp"
#include "point2d.hpp"
#include "length.hpp"

using namespace EikonalXX::Ray;

class Segment2D::Segment2DImpl
{
public:
    Point2D mStartPoint;
    Point2D mEndPoint;
    double mLength{0};
    double mSlowness{0};
    int mVelocityModelCellIndex{-1};
    bool mHaveStartAndEndPoint{false};
};

/// Constructor
Segment2D::Segment2D() :
    pImpl(std::make_unique<Segment2DImpl> ())
{
}

/// Copy constructor
Segment2D::Segment2D(const Segment2D &segment)
{
    *this = segment;
}

/// Move constructor
Segment2D::Segment2D(Segment2D &&segment) noexcept
{
    *this = std::move(segment);
}

/// Copy assignment
Segment2D& Segment2D::operator=(const Segment2D &segment)
{
    if (&segment == this){return *this;}
    pImpl = std::make_unique<Segment2DImpl> (*segment.pImpl);
    return *this;
}

/// Move assignment
Segment2D& Segment2D::operator=(Segment2D &&segment) noexcept
{
    if (&segment == this){return *this;}
    pImpl = std::move(segment.pImpl);
    return *this;
}

/// Reset class
void Segment2D::clear() noexcept
{
    pImpl = std::make_unique<Segment2DImpl> ();
}

/// Destructor
Segment2D::~Segment2D() = default;

/// Start and end point
void Segment2D::setStartAndEndPoint(
    const std::pair<Point2D, Point2D> &startAndEndPoint)
{
    if (!startAndEndPoint.first.havePositionInX())
    {
        throw std::invalid_argument("Start x position not set");
    }
    if (!startAndEndPoint.first.havePositionInZ())
    {
        throw std::invalid_argument("Start z position not set");
    }
    if (!startAndEndPoint.second.havePositionInX())
    {
        throw std::invalid_argument("End x position not set");
    }
    if (!startAndEndPoint.second.havePositionInZ())
    {
        throw std::invalid_argument("End z position not set");
    }
    pImpl->mStartPoint = startAndEndPoint.first;
    pImpl->mEndPoint = startAndEndPoint.second;
    pImpl->mLength = ::computeLength(pImpl->mStartPoint, pImpl->mEndPoint);
    pImpl->mHaveStartAndEndPoint = true;
}

Point2D Segment2D::getStartPoint() const
{
    if (!haveStartAndEndPoint())
    {   
        throw std::runtime_error("Start point not set");
    }   
    return pImpl->mStartPoint;
}

Point2D Segment2D::getEndPoint() const
{
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("End point not set");
    }
    return pImpl->mEndPoint;
}

double Segment2D::getLength() const
{
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("Start and end point not set");
    }
    return pImpl->mLength;
}

bool Segment2D::haveStartAndEndPoint() const noexcept
{
    return pImpl->mHaveStartAndEndPoint;
}

/// Velocity
void Segment2D::setVelocity(const double velocity)
{
    if (velocity <= 0)
    {
        throw std::runtime_error("Velocity = " + std::to_string(velocity)
                               + " must be positive");
    }
    setSlowness(1./velocity);
}

void Segment2D::setSlowness(const double slowness)
{
    if (slowness <= 0)
    {
        throw std::runtime_error("Slowness = " + std::to_string(slowness)
                              + " must be positive");}
    pImpl->mSlowness = slowness;
}

double Segment2D::getVelocity() const
{
    if (!haveVelocity()){throw std::runtime_error("Velocity not set");}
    return 1./getSlowness();
}

double Segment2D::getSlowness() const
{
    return pImpl->mSlowness;
}

bool Segment2D::haveVelocity() const noexcept
{
    return (pImpl->mSlowness > 0);
}

/// Travel time
double Segment2D::getTravelTime() const
{
    return getSlowness()*getLength();
}

/// Cell index
void Segment2D::setVelocityModelCellIndex(const int cellIndex)
{
    if (cellIndex < 0)
    {
        throw std::invalid_argument("Cell index must be positive");
    }
    pImpl->mVelocityModelCellIndex = cellIndex;
}

int Segment2D::getVelocityModelCellIndex() const
{
    if (!haveVelocityModelCellIndex())
    {
        throw std::runtime_error("Velocity model cell index not set");
    }
    return pImpl->mVelocityModelCellIndex;
}

bool Segment2D::haveVelocityModelCellIndex() const noexcept
{
    return (pImpl->mVelocityModelCellIndex >= 0);
}

void EikonalXX::Ray::swap(Segment2D &lhs, Segment2D &rhs) noexcept
{
    std::swap(lhs.pImpl, rhs.pImpl);
}

/// Reverses the segments
void Segment2D::reverse()
{   
    if (!haveStartAndEndPoint())
    {
        throw std::runtime_error("Start and end point not set");
    }
    std::swap(pImpl->mStartPoint, pImpl->mEndPoint);
} 

std::ostream& EikonalXX::Ray::operator<<(std::ostream &os,
                                         const Segment2D &segment)
{
    std::string result;
    if (segment.haveStartAndEndPoint())
    {   
        auto point0 = segment.getStartPoint();
        auto point1 = segment.getEndPoint();
        result = "StartPoint : ("
               + std::to_string(point0.getPositionInX()) + ","
               + std::to_string(point0.getPositionInZ()) + ")\n"
               + "EndPoint : ("
               + std::to_string(point1.getPositionInX()) + ","
               + std::to_string(point1.getPositionInZ()) + ")\n"
               + "Length : " + std::to_string(segment.getLength()); 
    }
    if (segment.haveVelocity())
    {
        if (!result.empty()){result = result + "\n";}
        result = result 
               + "Velocity : " + std::to_string(segment.getVelocity());
        if (segment.haveStartAndEndPoint())
        {
            result = result + "\nTravelTime : "
                   + std::to_string(segment.getTravelTime());
        }
    }
    if (segment.haveVelocityModelCellIndex())
    {
        if (!result.empty()){result = result + "\n";}
        result = result + "VelocityModelIndex : "
               + std::to_string(segment.getVelocityModelCellIndex());
    }
    return os << result;
}

