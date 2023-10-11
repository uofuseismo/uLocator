#include <string>
#include <algorithm>
#include "uLocator/station.hpp"
#include "uLocator/position/wgs84.hpp"

using namespace ULocator;

namespace
{
[[nodiscard]] std::string removeBlanksAndCapitalize(const std::string &s) 
{
    auto result = s;
    result.erase(std::remove_if(result.begin(), result.end(), ::isspace),
                 result.end());
    std::transform(result.begin(), result.end(), result.begin(), ::toupper);
    return result;
}
}

class Station::StationImpl
{
public:
    void updateHash()
    {
        mHash = std::hash<std::string> {}(mNetwork + "." + mStation);
    }
    std::string mNetwork;
    std::string mStation;
    Position::WGS84 mGeographicPosition; 
    size_t mHash{0};
    double mElevation{0};
    bool mHaveElevation{false};
};

/// Constructor
Station::Station() :
    pImpl(std::make_unique<StationImpl> ())
{
}

/// Copy constructor
Station::Station(const Station &station)
{
    *this = station;
}

/// Move constructor
Station::Station(Station &&station) noexcept
{
    *this = std::move(station);
}

/// Copy assignment
Station& Station::operator=(const Station &station)
{
    if (&station == this){return *this;}
    pImpl = std::make_unique<StationImpl> (*station.pImpl);
    return *this;
}

/// Move assignment
Station& Station::operator=(Station &&station) noexcept
{
    if (&station == this){return *this;}
    pImpl = std::move(station.pImpl);
    return *this;
}

/// Reset class
void Station::clear() noexcept
{
    pImpl = std::make_unique<StationImpl> ();
}

/// Destructor
Station::~Station() = default;

/// Network code
void Station::setNetwork(const std::string &networkIn)
{
    auto network = ::removeBlanksAndCapitalize(networkIn);
    if (network.empty()){throw std::invalid_argument("Network is empty");}
    pImpl->mNetwork = network;
    pImpl->updateHash();
}

std::string Station::getNetwork() const
{
    if (!haveNetwork()){throw std::runtime_error("Network not set");}
    return pImpl->mNetwork;
}

bool Station::haveNetwork() const noexcept
{
   return !pImpl->mNetwork.empty(); 
}

/// Station name
void Station::setName(const std::string &nameIn)
{
    auto name = ::removeBlanksAndCapitalize(nameIn);
    if (name.empty()){throw std::invalid_argument("Name is empty");}
    pImpl->mStation = name;
    pImpl->updateHash();
}

std::string Station::getName() const
{
    if (!haveName()){throw std::runtime_error("Station name not set");}
    return pImpl->mStation;
}

bool Station::haveName() const noexcept
{
    return !pImpl->mStation.empty();
}

/// Hash
size_t Station::getHash() const
{
    if (haveName() && haveNetwork()){return pImpl->mHash;}
    throw std::runtime_error("Both network code and station name must be set");
}

/// Station position
void Station::setGeographicPosition(const Position::WGS84 &position)
{
    if (!position.havePosition())
    {
        throw std::invalid_argument("Position not set");
    }
    pImpl->mGeographicPosition = position;
}

Position::WGS84 Station::getGeographicPosition() const
{
    if (!haveGeographicPosition()){std::runtime_error("Position not set");}
    return pImpl->mGeographicPosition;
}

const Position::WGS84& Station::getGeographicPositionReference() const
{
    if (!haveGeographicPosition()){std::runtime_error("Position not set");}
    return *&pImpl->mGeographicPosition;
}

bool Station::haveGeographicPosition() const noexcept
{
    return pImpl->mGeographicPosition.havePosition();
}

/// Elevation
void Station::setElevation(const double elevation)
{
    if (elevation < -6400000)
    {
        throw std::invalid_argument("Elevation cannot be in center of earth");
    }
    if (elevation > 8860)
    {
        throw std::invalid_argument(
            "Elevation cannot exceed Mt Everest summit");
    }
    pImpl->mElevation = elevation; 
    pImpl->mHaveElevation = true;
}

double Station::getElevation() const
{
    if (!haveElevation()){throw std::runtime_error("Elevation not set");}
    return pImpl->mElevation;
}

bool Station::haveElevation() const noexcept
{
    return pImpl->mHaveElevation;
}
