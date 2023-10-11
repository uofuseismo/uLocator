#include <string>
#ifndef NDEBUG
#include <cassert>
#endif
#include "uLocator/arrival.hpp"
#include "uLocator/station.hpp"

using namespace ULocator;

class Arrival::ArrivalImpl
{
public:
    Station mStation;
    std::string mPhase;
    double mTime{0};
    double mResidual{0};
    double mStandardError{0.1};
    double mDistance{0};
    double mAzimuth{0};
    double mBackAzimuth{0};
    int64_t mIdentifier{0};
    bool mHaveDistance{false};
    bool mHaveTime{false};
    bool mHaveResidual{false};
    bool mHaveStandardError{false};
    bool mHaveStation{false};
    bool mHaveAzimuth{false};
    bool mHaveBackAzimuth{false};
};

/// Constructor
Arrival::Arrival() :
    pImpl(std::make_unique<ArrivalImpl> ())
{
}

/// Copy constructor
Arrival::Arrival(const Arrival &arrival)
{
    *this = arrival;
}

/// Move constructor
Arrival::Arrival(Arrival &&arrival) noexcept
{
    *this = std::move(arrival);
}

/// Copy assignment
Arrival& Arrival::operator=(const Arrival &arrival)
{
    if (&arrival == this){return *this;}
    pImpl = std::make_unique<ArrivalImpl> (*arrival.pImpl);
    return *this;
}

/// Move assignment
Arrival& Arrival::operator=(Arrival &&arrival) noexcept
{
    if (&arrival == this){return *this;}
    pImpl = std::move(arrival.pImpl);
    return *this;
}

/// Destructor
Arrival::~Arrival() = default;

/// Reset class
void Arrival::clear() noexcept
{
    pImpl = std::make_unique<ArrivalImpl> ();
}

/// Phase
void Arrival::setPhase(const PhaseType type) noexcept
{
    if (type == PhaseType::P)
    {
        setPhase("P");
    }
    else if (type == PhaseType::S)
    {
        setPhase("S");
    }
#ifndef NDEBUG
    else
    {
        assert(false);
    }
#endif
}

void Arrival::setPhase(const std::string &phase)
{
    if (phase.empty()){throw std::invalid_argument("Phase is empty");}
    pImpl->mPhase = phase;
}

std::string Arrival::getPhase() const
{
    if (!havePhase()){throw std::runtime_error("Phase not set");}
    return pImpl->mPhase;
}

bool Arrival::havePhase() const noexcept
{
    return !pImpl->mPhase.empty();
}

/// Time
void Arrival::setTime(const double time) noexcept
{
    pImpl->mTime = time;
    pImpl->mHaveTime = true;
}

double Arrival::getTime() const
{
    if (!haveTime()){throw std::runtime_error("Time not set");}
    return pImpl->mTime;
}

bool Arrival::haveTime() const noexcept
{
    return pImpl->mHaveTime;
}

/// Identifier
void Arrival::setIdentifier(const int64_t identifier) noexcept
{
    pImpl->mIdentifier = identifier;
}

int64_t Arrival::getIdentifier() const noexcept
{
    return pImpl->mIdentifier;
}

/// Station
void Arrival::setStation(const Station &station)
{
    if (!station.haveGeographicPosition())
    {
        throw std::invalid_argument("Geographic position not set");
    }
    if (!station.haveNetwork())
    {
        throw std::invalid_argument("Network code not set");
    }
    if (!station.haveName())
    {
        throw std::invalid_argument("Station code not set");
    }
    pImpl->mStation = station;
    pImpl->mHaveStation = true;
}

Station Arrival::getStation() const
{
    if (!haveStation()){throw std::runtime_error("Station not set");}
    return pImpl->mStation;
}

const Station& Arrival::getStationReference() const
{
    if (!haveStation()){throw std::runtime_error("Station not set");}
    return *&pImpl->mStation;
}

bool Arrival::haveStation() const noexcept
{
    return pImpl->mHaveStation;
}

/// Standard error
void Arrival::setStandardError(const double standardError)
{
    if (standardError <= 0)
    {
        throw std::invalid_argument("Standard error must be positive");
    }
    pImpl->mStandardError = standardError;
    pImpl->mHaveStandardError = true;
}

double Arrival::getStandardError() const
{
    if (!haveStandardError())
    {
        throw std::runtime_error("Standard error not set");
    }
    return pImpl->mStandardError;
}

bool Arrival::haveStandardError() const noexcept
{
    return pImpl->mHaveStandardError;
}

/// Travel time residual
void Arrival::setResidual(const double residual) noexcept
{
    pImpl->mResidual = residual;
    pImpl->mHaveResidual = true;
}

double Arrival::getResidual() const
{
    if (!haveResidual()){throw std::runtime_error("Residual not set");}
    return pImpl->mResidual;
}

bool Arrival::haveResidual() const noexcept
{
    return pImpl->mHaveResidual;
}

/// Source-receiver distance
void Arrival::setDistance(const double distance)
{
    if (distance < 0)
    {
        throw std::invalid_argument("Source receiver distance must positive");
    }
    pImpl->mDistance = distance;
    pImpl->mHaveDistance = true;
}

double Arrival::getDistance() const
{
    if (!haveDistance())
    {
        throw std::runtime_error("Source-receiver distance not set");
    }
    return pImpl->mDistance;
}

bool Arrival::haveDistance() const noexcept
{
    return pImpl->mHaveDistance;
} 

/// Source-to-receiver azimuth
void Arrival::setAzimuth(const double azimuth)
{
    if (azimuth < 0 || azimuth >= 360)
    {   
        throw std::invalid_argument(
           "Source-receiver azimuth not in range [0,360)");
    }   
    pImpl->mAzimuth = azimuth;
    pImpl->mHaveAzimuth = true;
}

double Arrival::getAzimuth() const
{
    if (!haveAzimuth())
    {
        throw std::runtime_error("Source-receiver azimuth not set");
    }
    return pImpl->mAzimuth;
}

bool Arrival::haveAzimuth() const noexcept
{
    return pImpl->mHaveAzimuth;
} 

/// Receiver-to-source azimuth
void Arrival::setBackAzimuth(const double backAzimuth)
{
    if (backAzimuth < 0 || backAzimuth >= 360)
    {
        throw std::invalid_argument(
           "Receiver-to-source azimuth of " + std::to_string(backAzimuth)
         + " not in range [0,360)");
    }
    pImpl->mBackAzimuth = backAzimuth;
    pImpl->mHaveBackAzimuth = true;
}

double Arrival::getBackAzimuth() const
{
    if (!haveBackAzimuth())
    {
        throw std::runtime_error("Receiver-to-source azimuth not set");
    }
    return pImpl->mBackAzimuth;
}

bool Arrival::haveBackAzimuth() const noexcept
{
    return pImpl->mHaveBackAzimuth;
}

