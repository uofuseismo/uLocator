#include "uLocator/rayTracer.hpp"
#include "uLocator/uussRayTracer.hpp"
#include "uLocator/station.hpp"
#include "uLocator/origin.hpp"
#include "uLocator/position/utahRegion.hpp"
#include "uLocator/position/ynpRegion.hpp"
#include "uLocator/position/wgs84.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace
{

ULocator::Station toStationYNP(const std::string &network,
                               const std::string &name,
                               const double latitude,
                               const double longitude,
                               const double elevation,
                               const int utmZone = 12) 
{
    ULocator::Position::YNPRegion ynp;
    ULocator::Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        ULocator::Position::WGS84 {latitude, longitude, utmZone}, ynp);
    station.setElevation(elevation);
    return station;
}

ULocator::Station toStation(const std::string &network,
                            const std::string &name,
                            const double latitude,
                            const double longitude,
                            const double elevation,
                            const int utmZone = 12)
{
    ULocator::Position::UtahRegion utah;
    ULocator::Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        ULocator::Position::WGS84 {latitude, longitude, utmZone}, utah);
    station.setElevation(elevation);
    return station;
}

}

TEST_CASE("ULocator::RayTracer", "[uncorrectedTravelTimes]")
{
    ULocator::Origin referenceOrigin;
    referenceOrigin.setTime(1646479475.880001);
    referenceOrigin.setEpicenter(ULocator::Position::WGS84 {38.5615, -112.2196667} );
    referenceOrigin.setDepth(2000);
    ULocator::Position::UtahRegion utah;
    auto [xSource, ySource] = utah.geographicToLocalCoordinates(38.5615, -112.2196667);

    std::vector<ULocator::Station> stations;
    stations.push_back(::toStation("UU", "ASU7", 38.552595, -112.220772, 1609));
    stations.push_back(::toStation("LB", "MVU" , 38.5037,   -112.212303, 2239));
    stations.push_back(::toStation("UU", "TCRU", 38.6095,   -112.4472,   2293));
    stations.push_back(::toStation("UU", "WCU",  38.96467,  -112.09067,  2673));
    stations.push_back(::toStation("UU", "NMU",  38.5165,   -112.85,     1853));
    stations.push_back(::toStation("UU", "ECUT", 39.171697, -112.133201, 2136));
    stations.push_back(::toStation("UU", "DWU",  38.10534,  -112.9975,   2270));
    stations.push_back(::toStation("UU", "CVRU", 38.9176,   -111.1716,   1912));
    stations.push_back(::toStation("UU", "SWUT", 39.3286,   -113.1954,   1644));
    const std::vector<double> referenceTimes{0.8452878124706655,
                                             1.704260847820826,
                                             4.057303365184186,
                                             8.478570990177197,
                                             9.812868545764282,
                                             12.07825821309153,
                                             14.92522075012022,
                                             17.30671308165072,
                                             20.75111804581613};

    const std::string phase{"P"};
    const std::vector<double> utahInterfaces{-4500,  40,  15600, 26500, 40500};
    const std::vector<double> utahPVelocities{3400, 5900,  6400,  7500,  7900};
    for (int i = 0; i < static_cast<int> (stations.size()); ++i)
    {
        ULocator::RayTracer rayTracer(stations[i],
                                      phase,
                                      utahInterfaces,
                                      utahPVelocities);
        constexpr bool applyCorrection{false};
        auto travelTime = rayTracer.evaluate(0, xSource, ySource, 
                                             referenceOrigin.getDepth(),
                                             nullptr, nullptr, nullptr, nullptr,
                                             applyCorrection);
        REQUIRE_THAT(travelTime,
                     Catch::Matchers::WithinAbs(referenceTimes[i], 1.e-4));
        travelTime = rayTracer.evaluate(1, xSource, ySource,
                                        referenceOrigin.getDepth(),
                                        nullptr, nullptr, nullptr, nullptr,
                                        applyCorrection);
        REQUIRE_THAT(travelTime,
                     Catch::Matchers::WithinAbs(1 + referenceTimes[i], 1.e-4));
        //std::cout << std::setprecision(16) << travelTime << std::endl;
    }
}

TEST_CASE("ULocator::RayTracer", "[uncorrectedUtahToLayerTracerTimes]")
{
    ULocator::Position::UtahRegion utah;
    auto [xSource, ySource] = utah.geographicToLocalCoordinates(38.5615, -112.2196667);

    std::vector<ULocator::Station> stations;
    stations.push_back(::toStation("UU", "ASU7", 38.552595, -112.220772, 1609));
    stations.push_back(::toStation("LB", "MVU" , 38.5037,   -112.212303, 2239));
    stations.push_back(::toStation("UU", "TCRU", 38.6095,   -112.4472,   2293));
    stations.push_back(::toStation("UU", "WCU",  38.96467,  -112.09067,  2673));
    stations.push_back(::toStation("UU", "NMU",  38.5165,   -112.85,     1853));
    stations.push_back(::toStation("UU", "ECUT", 39.171697, -112.133201, 2136));
    stations.push_back(::toStation("UU", "DWU",  38.10534,  -112.9975,   2270));
    stations.push_back(::toStation("UU", "CVRU", 38.9176,   -111.1716,   1912));
    stations.push_back(::toStation("UU", "SWUT", 39.3286,   -113.1954,   1644));

    const std::vector<double> depths{-4200, 20, 25000, 32000, 55000};
    const std::vector<double> utahInterfaces{-4500,  40,  15600, 26500, 40500};
    const std::vector<double> utahPVelocities{3664, 5965, 6566, 7100, 7856};
    const std::vector<double> utahSVelocities{2004, 3480, 3858, 4168, 5006};
    for (int i = 0; i < static_cast<int> (stations.size()); ++i)
    {   
        ULocator::RayTracer rayTracer(stations[i],
                                      "P",
                                      utahInterfaces,
                                      utahPVelocities);
        ULocator::UUSSRayTracer pUtah(stations[i],
                                      ULocator::UUSSRayTracer::Phase::P,
                                      ULocator::UUSSRayTracer::Region::Utah);

        constexpr bool applyCorrection{false};
        for (const auto &depth : depths)
        {
            double dtdt0Ref, dtdxRef, dtdyRef, dtdzRef;
            auto travelTimeRef
                = rayTracer.evaluate(1, xSource, ySource, depth,
                                     &dtdt0Ref, &dtdxRef, &dtdyRef, &dtdzRef,
                                     applyCorrection);
            double dtdt0, dtdx, dtdy, dtdz; 
            auto travelTime
                = pUtah.evaluate(1, xSource, ySource, depth,
                                 &dtdt0, &dtdx, &dtdy, &dtdz,
                                 applyCorrection);
            REQUIRE_THAT(travelTimeRef,
                         Catch::Matchers::WithinAbs( travelTime, 1.e-8));
            REQUIRE_THAT(dtdt0Ref,
                         Catch::Matchers::WithinAbs( dtdt0,      1.e-8));
            REQUIRE_THAT(dtdxRef,
                         Catch::Matchers::WithinAbs( dtdx,       1.e-8));
            REQUIRE_THAT(dtdyRef,
                         Catch::Matchers::WithinAbs( dtdy,       1.e-8));
            REQUIRE_THAT(dtdzRef,
                         Catch::Matchers::WithinAbs( dtdz,       1.e-8));
        }
    }

    for (int i = 0; i < static_cast<int> (stations.size()); ++i)
    {   
        ULocator::RayTracer rayTracer(stations[i],
                                      "S",
                                      utahInterfaces,
                                      utahSVelocities);
        ULocator::UUSSRayTracer sUtah(stations[i],
                                      ULocator::UUSSRayTracer::Phase::S,
                                      ULocator::UUSSRayTracer::Region::Utah);

        constexpr bool applyCorrection{false};
        for (const auto &depth : depths)
        {
            double dtdt0Ref, dtdxRef, dtdyRef, dtdzRef;
            auto travelTimeRef
                = rayTracer.evaluate(1, xSource, ySource, depth,
                                     &dtdt0Ref, &dtdxRef, &dtdyRef, &dtdzRef,
                                     applyCorrection);
            double dtdt0, dtdx, dtdy, dtdz; 
            auto travelTime
                = sUtah.evaluate(1, xSource, ySource, depth,
                                 &dtdt0, &dtdx, &dtdy, &dtdz,
                                 applyCorrection);
            REQUIRE_THAT(travelTimeRef,
                         Catch::Matchers::WithinAbs( travelTime, 1.e-8));
            REQUIRE_THAT(dtdt0Ref,
                         Catch::Matchers::WithinAbs( dtdt0,      1.e-8));
            REQUIRE_THAT(dtdxRef, 
                         Catch::Matchers::WithinAbs( dtdx,       1.e-8));
            REQUIRE_THAT(dtdyRef,
                         Catch::Matchers::WithinAbs( dtdy,       1.e-8));
            REQUIRE_THAT(dtdzRef,
                         Catch::Matchers::WithinAbs( dtdz,       1.e-8));
        }
    }   
}

TEST_CASE("ULocator::RayTracer", "[uncorrectedYNPToLayerTracerTimes]")
{
    ULocator::Position::UtahRegion ynp;
    auto [xSource, ySource] = ynp.geographicToLocalCoordinates(44.7433333, -111.0675);

    std::vector<ULocator::Station> stations;
    stations.push_back(::toStationYNP("WY", "YHB", 44.7508, -111.1962,  2167.0));
    stations.push_back(::toStationYNP("WY", "YMR", 44.66867,-110.965,   2149.0));
    stations.push_back(::toStationYNP("WY", "YDC", 44.7095, -111.23967, 2025.0));
    stations.push_back(::toStationYNP("WY", "YHH", 44.78833,-110.8505,  2717.0));

    const std::vector<double> depths{-4200, 20, 3000, 7000, 13000, 18000, 24000, 55000};
    const std::vector<double> ynpInterfaces{-4500, -1000,  2000,  5000,  8000,
                                            12000, 16000, 21000, 50000};
    const std::vector<double> ynpPVelocities{2713, 2787, 5207, 5565, 5812,
                                             6150, 6290, 6611, 7990};
    const std::vector<double> ynpSVelocities{1745, 2316, 3066, 3375, 3529,
                                             3650, 3713, 4014, 5054};
    for (int i = 0; i < static_cast<int> (stations.size()); ++i)
    {
        ULocator::RayTracer rayTracer(stations[i],
                                      "P",
                                      ynpInterfaces,
                                      ynpPVelocities);
        ULocator::UUSSRayTracer pYNP(stations[i],
                                     ULocator::UUSSRayTracer::Phase::P,
                                     ULocator::UUSSRayTracer::Region::YNP);

        constexpr bool applyCorrection{false};
        for (const auto &depth : depths)
        {
            double dtdt0Ref, dtdxRef, dtdyRef, dtdzRef;
            auto travelTimeRef
                = rayTracer.evaluate(1, xSource, ySource, depth,
                                     &dtdt0Ref, &dtdxRef, &dtdyRef, &dtdzRef,
                                     applyCorrection);
            double dtdt0, dtdx, dtdy, dtdz;
            auto travelTime
                = pYNP.evaluate(1, xSource, ySource, depth,
                                &dtdt0, &dtdx, &dtdy, &dtdz,
                                applyCorrection);
            REQUIRE_THAT(travelTimeRef,
                         Catch::Matchers::WithinAbs( travelTime, 1.e-8));
            REQUIRE_THAT(dtdt0Ref,
                         Catch::Matchers::WithinAbs( dtdt0,      1.e-8));
            REQUIRE_THAT(dtdxRef,
                         Catch::Matchers::WithinAbs( dtdx,       1.e-8));
            REQUIRE_THAT(dtdyRef,
                         Catch::Matchers::WithinAbs( dtdy,       1.e-8));
            REQUIRE_THAT(dtdzRef,
                         Catch::Matchers::WithinAbs( dtdz,       1.e-8));
        }
    }

    for (int i = 0; i < static_cast<int> (stations.size()); ++i)
    {
        ULocator::RayTracer rayTracer(stations[i],
                                      "S",
                                      ynpInterfaces,
                                      ynpSVelocities);
        ULocator::UUSSRayTracer sYNP(stations[i],
                                     ULocator::UUSSRayTracer::Phase::S,
                                     ULocator::UUSSRayTracer::Region::YNP);

        constexpr bool applyCorrection{false};
        for (const auto &depth : depths)
        {
            double dtdt0Ref, dtdxRef, dtdyRef, dtdzRef;
            auto travelTimeRef
                = rayTracer.evaluate(1, xSource, ySource, depth,
                                     &dtdt0Ref, &dtdxRef, &dtdyRef, &dtdzRef,
                                     applyCorrection);
            double dtdt0, dtdx, dtdy, dtdz;
            auto travelTime
                = sYNP.evaluate(1, xSource, ySource, depth,
                                &dtdt0, &dtdx, &dtdy, &dtdz,
                                applyCorrection);
            REQUIRE_THAT(travelTimeRef,
                         Catch::Matchers::WithinAbs(travelTime, 1.e-8));
            REQUIRE_THAT(dtdt0Ref,
                         Catch::Matchers::WithinAbs(dtdt0,      1.e-8));
            REQUIRE_THAT(dtdxRef,
                         Catch::Matchers::WithinAbs(dtdx,       1.e-8));
            REQUIRE_THAT(dtdyRef,
                         Catch::Matchers::WithinAbs(dtdy,       1.e-8));
            REQUIRE_THAT(dtdzRef,
                         Catch::Matchers::WithinAbs(dtdz,       1.e-8));
        }
    }
}
