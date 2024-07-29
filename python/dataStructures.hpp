#ifndef ULOCATOR_PYTHON_DATA_STRUCTURES_HPP
#define ULOCATOR_PYTHON_DATA_STRUCTURES_HPP
#include <pybind11/pybind11.h>
#include <uLocator/arrival.hpp>
#include <uLocator/origin.hpp>
#include <uLocator/station.hpp>
namespace ULocator::Python::Position
{
 class WGS84;
 class IGeographicRegion;
}
namespace ULocator::Python
{
class Arrival
{
public:
    Arrival();
    Arrival(const Arrival &arrival);
    Arrival(Arrival &&arrival) noexcept;
    Arrival& operator=(const Arrival &arrival);
    Arrival& operator=(Arrival &&arrival) noexcept;

    void setTime(double time) noexcept;
    double getTime() const;
    void setStandardError(double error);
    double getStandardError() const;
    void setPhase(const std::string &phase);
    std::string getPhase() const;
    void setIdentifier(int64_t identifier) noexcept;
    int64_t getIdentifier() const;

    double getAzimuth() const;
    double getBackAzimuth() const;
    double getDistance() const;

    void clear() noexcept;
    ~Arrival();
private:
    std::unique_ptr<ULocator::Arrival> mArrival;
};

class Station
{
public:
    Station();
    Station(const Station &station);
    Station(Station &&station) noexcept;

    void setNetwork(const std::string &network);
    [[nodiscard]] std::string getNetwork() const;

    void setName(const std::string &name);
    [[nodiscard]] std::string getName() const;

    void setGeographicPosition(const ULocator::Python::Position::WGS84 &position,
                               const ULocator::Python::Position::IGeographicRegion &region);
    [[nodiscard]] std::pair<double, double> getLocalCoordinates() const;
    [[nodiscard]] ULocator::Python::Position::WGS84 getGeographicPosition() const;
    [[nodiscard]] const ULocator::Python::Position::WGS84& getGeographicPositionReference() const;
    [[nodiscard]] const ULocator::Station &getBaseClassReference() const noexcept;

    void setElevation(double elevation);
    [[nodiscard]] double getElevation() const;


    void clear() noexcept;
    ~Station();

    Station& operator=(const Station &); 
    Station& operator=(Station &&) noexcept;
private:
    std::unique_ptr<ULocator::Station> pImpl;
};

void initialize(pybind11::module &m);
}
#endif
