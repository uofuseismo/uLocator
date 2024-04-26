#ifndef LOAD_CATALOG_HPP
#define LOAD_CATALOG_HPP
#include <cmath>
#include <boost/tokenizer.hpp>
#include <uLocator/arrival.hpp>
#include <uLocator/origin.hpp>
#include <uLocator/station.hpp>
#include <uLocator/position/geographicRegion.hpp>
#include <uLocator/position/wgs84.hpp>
#ifdef WITH_UMPS
#include <umps/logging/standardOut.hpp>
#else
#include "logging/standardOut.hpp"
#endif
namespace
{

/// @result Converts an event type to a string
[[nodiscard]] std::string eventTypeToString(
    const ULocator::Origin::EventType eventType)
{
    if (eventType == ULocator::Origin::EventType::Earthquake)
    {
        return "eq";
    }
    else if (eventType == ULocator::Origin::EventType::QuarryBlast)
    {
        return "qb";
    }
    throw std::runtime_error("Unhandled event type");
}

[[nodiscard]] ULocator::Origin::EventType 
stringToEventType(const std::string &eventType)
{
    if (eventType == "eq"){return ULocator::Origin::EventType::Earthquake;}
    if (eventType == "qb"){return ULocator::Origin::EventType::QuarryBlast;}
    throw std::runtime_error("Unhandled event type of " + eventType);
}

/// @brief Utility routine to convert from a Jiggle quality to a standard error.
/// @param[in] quality   The quality reported in the AQMS database.
double qualityToStandardError(const double quality)
{
    if (std::abs(quality - 1) < 1.e-4)
    {
        return 0.03;
    }
    else if (std::abs(quality - 0.75) < 1.e-4)
    {
        return 0.06;
    }
    else if (std::abs(quality - 0.5) < 1.e-4)
    {
        return 0.15;
    }
    else if (std::abs(quality - 0.25) < 1.e-4)
    {
        return 0.30;
    }
#ifndef NDEBUG
    assert(false);
#endif
    return 1;
}

/// @result The origins in the catalog.
/// @param[in] fileName          The name of the CSV file with the entire 
///                              catalog (earthquakes and picks).
/// @param[in] eventsFileName    The name of the file with the list of
///                              events to locate.
/// @param[in] region            The geographic region.
/// @param[in] logger            A logger to report errors.
/// @param[in] catalogVersion    Indicates the catalog format version.
/// @result The origins in the 
std::vector<ULocator::Origin>
    loadCatalog(const std::string &fileName,
                const std::string &eventsFileName,
                const ULocator::Position::IGeographicRegion &region,
                std::shared_ptr<UMPS::Logging::ILog> logger,
                const int catalogVersion = 3,
                const int utmZone = 12)
{
    // Get events to locate
    std::string line;
    std::vector<int64_t> eventIdentifiers;
    if (std::filesystem::exists(eventsFileName))
    {
        std::ifstream eventsFile(eventsFileName);
        getline(eventsFile, line); // Header
        while (getline(eventsFile, line))
        {
            std::vector<std::string> splitLine;
            splitLine.reserve(64);
            boost::tokenizer<boost::escaped_list_separator<char>>
               tokenizer(line,
                         boost::escaped_list_separator<char>('\\', ',', '\"'));
            for (auto i = tokenizer.begin(); i != tokenizer.end(); ++i) 
            {
                splitLine.push_back(*i);
            }
            eventIdentifiers.push_back(std::stol(splitLine[0]));
        }
    } 
    else
    {
        logger->warn("Will locate all events in catalog...");
    }
    // Get entire catalog
    std::vector<ULocator::Origin> origins;
    origins.reserve(15000);
    std::ifstream csvFile(fileName); 
    if (!csvFile.is_open())
    {
        throw std::runtime_error("Couldn't open: " + fileName);
    }
    ULocator::Origin origin;
    std::vector<ULocator::Arrival> arrivals;
    int64_t evidOld{-1};
    getline(csvFile, line); // Header
    while (getline(csvFile, line))
    {
        std::vector<std::string> splitLine;
        splitLine.reserve(64);
        boost::tokenizer<boost::escaped_list_separator<char>>
           tokenizer(line,
                     boost::escaped_list_separator<char>('\\', ',', '\"'));
        for (auto i = tokenizer.begin(); i != tokenizer.end(); ++i) 
        {
            splitLine.push_back(*i);
        }
        int64_t evid;
        std::string network;
        std::string stationName;
        std::string phase;
        int64_t arid;
        double time;
        double residual;
        double uncertainty;
        double stationLat;
        double stationLon;
        double stationElev;
        double eventLat;
        double eventLon;
        double eventDepth;
        double originTime;
        double sourceReceiverDistance{0};
        std::string eventType;
        if (catalogVersion == 1)
        {
            evid = std::stol(splitLine[0]);
            network = splitLine[1];
            stationName = splitLine[2];
            phase = splitLine[5];
            arid = std::stol(splitLine[6]);
            time = std::stod(splitLine[7]);
            residual = std::stod(splitLine[13]);
            uncertainty = ::qualityToStandardError(std::stod(splitLine[8]));
            stationLat = std::stod(splitLine[14]);
            stationLon = std::stod(splitLine[15]);
            stationElev = std::stod(splitLine[16]);
            eventLat = std::stod(splitLine[17]); 
            eventLon = std::stod(splitLine[18]);
            eventDepth = 1000*std::stod(splitLine[19]);
            originTime = std::stod(splitLine[20]);
            eventType = splitLine[24];
        }
        else if (catalogVersion == 2)
        {
            evid = std::stol(splitLine[0]);
            network = splitLine[1];
            stationName = splitLine[2];
            phase = splitLine[7];
            arid = std::stol(splitLine[12]);
            time = std::stod(splitLine[5]);
            residual = std::stod(splitLine[6]);
            uncertainty = std::stod(splitLine[11]);
            stationLat = std::stod(splitLine[8]);
            stationLon = std::stod(splitLine[9]);
            stationElev = std::stod(splitLine[10]);
            eventLat = std::stod(splitLine[14]); 
            eventLon = std::stod(splitLine[15]);
            eventDepth = 1000*std::stod(splitLine[16]);
            originTime = std::stod(splitLine[17]);
            eventType = splitLine[21];
        }
        else if (catalogVersion == 3 || catalogVersion == 4)
        {
//   0             1          2           3        4         5     6                  7       8       9       10            11    12           13           14      15       16               17                 18
//event_identifier,event_type,origin_time,latitude,longitude,depth,arrival_identifier,network,station,channel,location_code,phase,arrival_time,first_motion,quality,residual,station_latitude,station_longitude,station_elevation
            evid = std::stol(splitLine[0]);
            eventType = splitLine[1];
            originTime = std::stod(splitLine[2]);
            eventLat = std::stod(splitLine[3]);
            eventLon = std::stod(splitLine[4]);
            eventDepth = 1000*std::stod(splitLine[5]);
            arid = std::stol(splitLine[6]);
            network = splitLine[7];
            stationName = splitLine[8];
            phase = splitLine[11];
            time = std::stod(splitLine[12]);
            if (catalogVersion == 3)
            {
                uncertainty = ::qualityToStandardError(std::stod(splitLine[14]));
            }
            else
            {
                uncertainty = std::stod(splitLine[14]);
            }
            residual = std::stod(splitLine[15]);
            stationLat = std::stod(splitLine[16]);
            stationLon = std::stod(splitLine[17]);
            stationElev = std::stod(splitLine[18]);
            // Compute distance
            if (catalogVersion == 4)
            {
                ULocator::Position::WGS84
                     sourceLocation{eventLat, eventLon, utmZone};
                ULocator::Position::WGS84
                     receiverLocation{stationLat, stationLon, utmZone};
                ULocator::Position::computeDistanceAzimuth(sourceLocation,
                                                           receiverLocation,
                                                           nullptr,
                                                           &sourceReceiverDistance,
                                                           nullptr,
                                                           nullptr);
            }
        }
        else
        {
            throw std::runtime_error("Undefined catalog format");
        }
        if (evid != evidOld)
        {
            // When we get to a new event finish out this origin
            if (evidOld !=-1)
            {
                // Save the event if it's in our list or, if there's no list,
                // save the event b/c we'll locate the entire catalog
                bool saveIt{false};
                if (!eventIdentifiers.empty())
                {
                    if (std::find(eventIdentifiers.begin(),
                                  eventIdentifiers.end(),
                                  origin.getIdentifier())
                       != eventIdentifiers.end())
                    {
                       saveIt = true;
                    }
                }
                else
                {
                    saveIt = true;
                }
                if (saveIt && arrivals.size() < 4)
                {
                    logger->warn("Too few arrivals for " + std::to_string(evid)
                               + "- skipping");
                    saveIt = false;
                }
                if (saveIt)
                {
                    origin.setArrivals(arrivals);
                    origins.push_back(origin);
                }
            }
            // Update
            evidOld = evid;
            arrivals.clear();
            origin.clear();
            // Set initial origin information
            ULocator::Position::WGS84 epicenter{eventLat, eventLon, utmZone};
            origin.setIdentifier(evid);
            origin.setEpicenter(epicenter);
            origin.setDepth(eventDepth);
            origin.setTime(originTime);
            origin.setEventType(::stringToEventType(eventType));
        }
        ULocator::Position::WGS84 stationPosition{stationLat,
                                                  stationLon,
                                                  utmZone};
        ULocator::Arrival arrival;
        ULocator::Station station;
        station.setNetwork(network);
        station.setName(stationName);
        station.setGeographicPosition(stationPosition, region);
        station.setElevation(stationElev);
        arrival.setTime(time);
        arrival.setStandardError(uncertainty);
        arrival.setResidual(residual);
        arrival.setPhase(phase);
        arrival.setStation(station);
        arrival.setIdentifier(arid);
        // Assert the arrival does not already exist
        bool exists{false};
        for (const auto &a : arrivals)
        {
            const auto &si = a.getStation();
            if (si.getNetwork() == network &&
                si.getName()    == stationName &&
                a.getPhase()    == phase)
            {
                exists = true;
                logger->warn("Duplicate phase arrival for: "
                            + std::to_string(evid)
                            + " " + network + "." + stationName + "." + phase
                            + "; skipping...");
            }
        }
        if (!exists && sourceReceiverDistance < 300000)
        {
            arrivals.push_back(std::move(arrival));
        }
    }
    return origins;    
}

/// @param[in] origins  The list of origins.
/// @result The unique list of stations given all the picks defining the
///         origins.
std::vector<ULocator::Station> uniqueStationsFromOrigins(
    const std::vector<ULocator::Origin> &origins)
{
    std::vector<ULocator::Station> uniqueStations;
    std::vector<std::string> uniqueStationNames;
    for (const auto &origin : origins)
    {
        const auto arrivals = origin.getArrivalsReference();
        for (const auto &arrival : arrivals)
        {
            const auto station = arrival.getStationReference();
            auto name = station.getNetwork() + "." + station.getName();
            bool exists{false};
            for (const auto &uniqueStationName : uniqueStationNames)
            {
                if (name == uniqueStationName)
                {
                    exists = true;
                    break;
                }
            }
            if (!exists)
            {
                uniqueStationNames.push_back(name);
                uniqueStations.push_back(station);
            }
        }
    }
    return uniqueStations;
}

/// @result The weighted root-mean-squared error from the arrivals. 
[[nodiscard]]
double computeWeightedRMSE(const std::vector<ULocator::Arrival> &arrivals)
{
    if (arrivals.empty())
    {
        throw std::invalid_argument("No arrivals");
    }
    double numerator{0};
    double denominator{0};
    for (const auto &arrival : arrivals)
    {
        auto weight = 1./arrival.getStandardError();
        auto residual = arrival.getResidual();
        numerator = numerator + weight*(residual*residual);
        denominator = denominator + weight; 
    }
    if (denominator == 0)
    {
        throw std::runtime_error("All weigths are zero");
    }
    return std::sqrt(numerator/denominator);
} 

}
#endif
