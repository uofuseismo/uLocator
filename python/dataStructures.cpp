#include <uLocator/arrival.hpp>
#include <uLocator/origin.hpp>
#include "dataStructures.hpp"
#include "position.hpp"

using namespace ULocator::Python;

Arrival::Arrival() :
    mArrival(std::make_unique<ULocator::Arrival> ())
{
}

Arrival::Arrival(const Arrival &arrival)
{
    *this = arrival;
}

Arrival::Arrival(Arrival &&arrival) noexcept
{
    *this = std::move(arrival);
}

Arrival& Arrival::operator=(const Arrival &arrival)
{
    if (&arrival == this){return *this;}
    mArrival = std::make_unique<ULocator::Arrival> (*arrival.mArrival);
    return *this;
}

Arrival& Arrival::operator=(Arrival &&arrival) noexcept
{
    if (&arrival == this){return *this;}
    mArrival = std::move(arrival.mArrival);
    return *this;
}

void Arrival::setIdentifier(const int64_t identifier) noexcept
{
    mArrival->setIdentifier(identifier);
}

int64_t Arrival::getIdentifier() const
{
    return mArrival->getIdentifier();
}

void Arrival::setTime(const double time) noexcept
{
    mArrival->setTime(time);
}

double Arrival::getTime() const
{
    return mArrival->getTime();
}

void Arrival::setStandardError(const double error)
{
    mArrival->setStandardError(error);
}

double Arrival::getStandardError() const
{
    return mArrival->getStandardError();
}

void Arrival::setPhase(const std::string &phase)
{
    mArrival->setPhase(phase);
}

std::string Arrival::getPhase() const
{
    return mArrival->getPhase();
}

void Arrival::clear() noexcept
{
    mArrival->clear();
}

double Arrival::getAzimuth() const
{
    return mArrival->getAzimuth();
}

double Arrival::getBackAzimuth() const
{
    return mArrival->getBackAzimuth();
}

double Arrival::getDistance() const
{
    return mArrival->getDistance();
}

Arrival::~Arrival() = default;

//----------------------------------------------------------------------------//

Station::Station() :
    pImpl(std::make_unique<ULocator::Station> ()) 
{
}

Station::Station(const Station &station)
{
    *this = station;
}

Station::Station(Station &&station) noexcept
{
    *this = std::move(station);
}

void Station::clear() noexcept
{
    pImpl = std::make_unique<ULocator::Station> ();
}

Station::~Station() = default;

Station& Station::operator=(const Station &station)
{
    if (&station == this){return *this;}
    pImpl = std::make_unique<ULocator::Station> (*station.pImpl);
    return *this;
}

Station& Station::operator=(Station &&station) noexcept
{
    if (&station == this){return *this;}
    pImpl = std::move(station.pImpl);
    return *this;
}

void Station::setNetwork(const std::string &network)
{
    pImpl->setNetwork(network);
}

std::string Station::getNetwork() const
{
    return pImpl->getNetwork();
}

void Station::setName(const std::string &name)
{
    pImpl->setName(name);
}

std::string Station::getName() const
{
    return pImpl->getName();
}

void Station::setElevation(const double elevation)
{
    pImpl->setElevation(elevation);
}

double Station::getElevation() const
{
    return pImpl->getElevation();
}

void Station::setGeographicPosition(const ULocator::Python::Position::WGS84 &position,
                                    const ULocator::Python::Position::IGeographicRegion &region)
{
    pImpl->setGeographicPosition(
        position.getBaseClassReference(),
        *region.getBaseClass());
}

std::pair<double, double> Station::getLocalCoordinates() const
{
    return pImpl->getLocalCoordinates();
}

ULocator::Python::Position::WGS84 Station::getGeographicPosition() const
{
    ULocator::Python::Position::WGS84
       wgs84{pImpl->getGeographicPosition()};
    return std::move(wgs84);
}

const ULocator::Station &Station::getBaseClassReference() const noexcept
{
    return *pImpl;
}


//----------------------------------------------------------------------------//

void ULocator::Python::initialize(pybind11::module &m)
{
    
    ///-----------------------------Arrival-----------------------------------///
    pybind11::enum_<ULocator::Arrival::PhaseType> (m, "PhaseType")
        .value("P", ULocator::Arrival::PhaseType::P,
               "A direct compressional wave arrival.")
        .value("S", ULocator::Arrival::PhaseType::S,
               "A direct shear wave arrival.");
    pybind11::class_<ULocator::Python::Arrival> arrival(m, "Arrival");
    arrival.def(pybind11::init<> ());
    arrival.doc() = R""""(
This defines a seismic phase arrival.

Required Read-Write Properties:
   time : float
      The time (UTC) since the epoch of the phase arrival. 
   standard_error : float
      The standard error of the phase arrival in seconds.
   station : Station
      Has the station's network code and name on which the arrival was observed.
   phase : str
      The arrival's phase type.

Optional Read-Write Properties:
   identifier : int
      The unique arrival identifier.

Read-Only Properties:
   azimuth : float
      The source-to-receiver azimuth in degrees measured positive east of north.
   back_azimuth : float
      The receiver-to-source azimuth in degrees measured positive east of north.
   distance : float
      The source-receiver distance in meters.
)"""";
    arrival.def("__copy__", [](const Arrival &self)
    {
        return Arrival(self);
    });
    arrival.def("clear",
                &Arrival::clear,
                "Resets the class.");
    arrival.def_property_readonly("azimuth",
                                   &Arrival::getAzimuth);
    arrival.def_property_readonly("back_azimuth",
                                   &Arrival::getBackAzimuth);
    arrival.def_property_readonly("distance",
                                   &Arrival::getDistance);
    arrival.def_property("time",
                         &Arrival::getTime, &Arrival::setTime);
    arrival.def_property("standard_error",
                         &Arrival::getStandardError, &Arrival::setStandardError);
    arrival.def_property("phase",
                         &Arrival::getPhase, &Arrival::setPhase);
    arrival.def_property("identifier",
                         &Arrival::getIdentifier, &Arrival::setIdentifier);

    ///-----------------------------Station----------------------------------///
    pybind11::class_<ULocator::Python::Station> station(m, "Station");
    station.def(pybind11::init<> ());
    station.def("__copy__", [](const Station &self)
    {
        return Station(self);
    });
    station.def("clear",
                &Station::clear,
                "Resets the class.");
    station.def_property("network",
                         &Station::getNetwork, &Station::setNetwork);
    station.def_property("name",
                         &Station::getName, &Station::setName);
    station.def_property("elevation",
                         &Station::getElevation, &Station::setElevation);
    station.def("set_geographic_position",
                &Station::setGeographicPosition,
                "Sets the geographic position of the stationl");
}
