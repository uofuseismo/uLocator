#include <fstream>
#include <nlohmann/json.hpp>
#include "uLocator/station.hpp"
#include "uLocator/arrival.hpp"
#include "uLocator/position/wgs84.hpp"
#include <gtest/gtest.h>

namespace
{

using namespace ULocator;

TEST(ULocator, 60086357)
{
    constexpr int utmZone{12};
    std::ifstream f{"60086357.json"};
    auto jsonData = nlohmann::json::parse(f);
    f.close();
    for (const auto &element : jsonData)
    {
        //std::cout << element << std::endl;
        Station station;
        auto arrivalIdentifier = element["arrival_id"].get<int64_t> ();
        auto networkCode = element["network"].get<std::string> ();
        auto stationName = element["station"].get<std::string> ();
        auto phase = element["phase"].get<std::string> ();
        auto time = element["arrival_time"].get<double> ();
        auto stationLatitude = element["station_latitude"].get<double> ();
        auto stationLongitude = element["station_longitude"].get<double> ();
        //std::cout << networkCode << " " << stationName << " " << stationLatitude << " " << stationLongitude << std::endl;
        Arrival arrival;
        station.setNetwork(networkCode);
        station.setName(stationName);
        station.setGeographicPosition(
            Position::WGS84 {stationLatitude, stationLongitude, utmZone});
        arrival.setTime(time);
        arrival.setStation(station);
        arrival.setIdentifier(arrivalIdentifier); 
        if (phase == "P")
        {
            arrival.setPhase(Arrival::PhaseType::P);
            arrival.setStandardError(0.1);
        }
        else if (phase == "S")
        {
            arrival.setPhase(Arrival::PhaseType::S);
            arrival.setStandardError(0.2);
        }
        else
        {
            std::cerr << "Unhandled phase: " << phase << std::endl;
        }
//{"arrival_id":678909,"arrival_time":1411640134.6037428,"network":"US","phase":"P","station":"HWUT","station_latitude":41.6069,"station_longitude":-111.5652}

    }

}

}
