#include <string>
#include <uLocator/origin.hpp>
#include <uLocator/position/wgs84.hpp>
#include <uLocator/position/utahRegion.hpp>
#include <uLocator/station.hpp>
#include <uLocator/arrival.hpp>
#include <uLocator/uussRayTracer.hpp>
#include <uLocator/travelTimeCalculatorMap.hpp>
#include <uLocator/topography/constant.hpp>
#include <uLocator/optimizers/nlopt/dividedRectangles.hpp>
#include <uLocator/optimizers/nlopt/stochasticGradientOptimization.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

namespace
{

ULocator::Arrival toArrival(const std::string &network,
                            const std::string &name,
                            const std::string &phase,
                            const double latitude,
                            const double longitude,
                            const double elevation,
                            const double arrivalTime,
                            const double standardError, 
                            const int utmZone = 12) 
{
    ULocator::Position::UtahRegion utah;
    ULocator::Station station;
    station.setNetwork(network);
    station.setName(name);
    station.setGeographicPosition(
        ULocator::Position::WGS84 {latitude, longitude, utmZone}, utah);
    station.setElevation(elevation);

    ULocator::Arrival arrival;
    arrival.setStation(station);
    arrival.setPhase(phase);
    arrival.setTime(arrivalTime);
    arrival.setStandardError(standardError);
    return arrival;
}

/// Magna event uu60484462
std::vector<ULocator::Arrival> createEarthquakeArrivals()
{
    std::vector<ULocator::Arrival> arrivals;
    arrivals.push_back(::toArrival("UU", "ASU7", "P", 38.552595, -112.220772, 1609, 1646479476.6776488,   0.03));
    arrivals.push_back(::toArrival("LB", "MVU" , "P", 38.5037,   -112.212303, 2239, 1646479477.4650624,   0.03));
    arrivals.push_back(::toArrival("UU", "ASU7", "S", 38.552595, -112.220772, 1609, 1646479477.5781794,   0.15));
    arrivals.push_back(::toArrival("UU", "TCRU", "P", 38.6095,   -112.4472,   2293, 1646479479.779858,    0.03));
    arrivals.push_back(::toArrival("UU", "TCRU", "S", 38.6095,   -112.4472,   2293, 1646479483.044294,    0.15));
    arrivals.push_back(::toArrival("UU", "WCU",  "P", 38.96467,  -112.09067,  2673, 1646479483.8797174,   0.06));
    arrivals.push_back(::toArrival("UU", "NMU",  "P", 38.5165,   -112.85,     1853, 1646479485.68,        0.06));
    arrivals.push_back(::toArrival("UU", "ECUT", "P", 39.171697, -112.133201, 2136, 1646479487.7271104,   0.15));
    arrivals.push_back(::toArrival("UU", "DWU",  "P", 38.10534,  -112.9975,   2270, 1646479490.553661,    0.15));
    arrivals.push_back(::toArrival("UU", "NMU",  "S", 38.5165,   -112.85,     1853, 1646479492.619098,    0.06));
    arrivals.push_back(::toArrival("UU", "CVRU", "P", 38.9176,  -111.1716,    1912, 1646479492.699132,    0.15));
    arrivals.push_back(::toArrival("UU", "SWUT", "P", 39.3286,  -113.1954,    1644, 1646479496.371336,    0.15));
    return arrivals;
}

/// Quarry blast uu60484742
std::vector<ULocator::Arrival> createQuarryBlastArrivals()
{
    std::vector<ULocator::Arrival> arrivals;
    arrivals.push_back(::toArrival("UU", "CWU",  "P", 40.44584,-112.10217, 1945.0, 1646684879.079617,   0.06));
    arrivals.push_back(::toArrival("UU", "MID",  "P", 40.51733,-112.25467, 1722.0, 1646684879.3381984,  0.03));
    arrivals.push_back(::toArrival("UU", "NOQ",  "P", 40.6525, -112.12033, 1622.0, 1646684880.4258294,  0.03));
    arrivals.push_back(::toArrival("UU", "WTU",  "P", 40.45483,-111.9535,  1552.0, 1646684880.785694,   0.15));
    arrivals.push_back(::toArrival("UU", "GMU",  "P", 40.5755, -111.76317, 1829.0, 1646684883.5544589,  0.15));
    arrivals.push_back(::toArrival("UU", "SAIU", "P", 40.85483,-112.1815,  1384.0, 1646684883.9708252,  0.03));
    arrivals.push_back(::toArrival("UU", "RBU",  "P", 40.78083,-111.80833, 1676.0, 1646684885.032984,   0.15));
    arrivals.push_back(::toArrival("UU", "SNUT", "P", 40.88567,-112.509,   1652.0, 1646684886.6219738,  0.15));
    arrivals.push_back(::toArrival("UU", "NAIU", "P", 41.01617,-112.228,   1472.0, 1646684887.0673397,  0.15));
    arrivals.push_back(::toArrival("UU", "NLU",  "P", 39.95483,-112.075,   2036.0, 1646684888.7377942,  0.15));
    arrivals.push_back(::toArrival("UU", "MOUT", "P", 41.19894,-111.87953, 2748.0, 1646684891.405937,   0.15));
    arrivals.push_back(::toArrival("UU", "SPU",  "P", 41.30867,-112.44917, 2086.0, 1646684893.4537792,  0.15));
    arrivals.push_back(::toArrival("UU", "HVU",  "P", 41.77967,-112.775,   1609.0, 1646684902.3164327,  0.06));
    return arrivals;
}

}

TEST(ULocatorNLOptDividedRectangles, FixedDepthLp)
{
// TODO redo this with the original UUSS velocity models
    ULocator::Origin referenceOrigin;
    referenceOrigin.setTime(1646479475.880001);
    referenceOrigin.setEpicenter(ULocator::Position::WGS84 {38.5615, -112.2196667} );
    referenceOrigin.setDepth(2000);

    auto arrivals = createEarthquakeArrivals();
    auto calculatorMap = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    for (const auto &arrival : arrivals)
    {
        std::unique_ptr<ULocator::ITravelTimeCalculator> calculator;
        if (arrival.getPhase() == "P")
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::P, 
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        else
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::S,
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        calculatorMap->insert(arrival.getStation(), arrival.getPhase(),
                              std::move(calculator));
    }
    auto topography = std::make_unique<ULocator::Topography::Constant> ();
    topography->set(4000);
 
    ULocator::Optimizers::NLOpt::DividedRectangles direct;
    direct.setGeographicRegion(ULocator::Position::UtahRegion {});
    direct.setTravelTimeCalculatorMap(std::move(calculatorMap));
    direct.setTopography(std::move(topography));
    direct.setArrivals(arrivals);
    EXPECT_TRUE(direct.locallyBias());
    EXPECT_TRUE(direct.normalize());
    EXPECT_NEAR(direct.getLocationTolerance(), 1, 1.e-8);
    EXPECT_NEAR(direct.getOriginTimeTolerance(), 0.01, 1.e-8);
    EXPECT_EQ(direct.getMaximumNumberOfObjectiveFunctionEvaluations(), 5000);

    ULocator::Origin noGuess;
    EXPECT_NO_THROW(direct.setMaximumNumberOfObjectiveFunctionEvaluations(2500));
    EXPECT_NO_THROW(
        direct.locateEventWithFixedDepth(
            referenceOrigin.getDepth(),
            ULocator::Optimizers::IOptimizer::Norm::LeastSquares )
    );
    EXPECT_TRUE(direct.haveOrigin());
    auto origin = direct.getOrigin();
    auto dLatitude = referenceOrigin.getEpicenter().getLatitude()
                   - origin.getEpicenter().getLatitude();
    auto dLongitude = referenceOrigin.getEpicenter().getLongitude()
                    - origin.getEpicenter().getLongitude();
    auto dDepth = referenceOrigin.getDepth() - origin.getDepth();
    auto dTime = referenceOrigin.getTime() - origin.getTime();
    EXPECT_NEAR(dLatitude,  0, 1.e-2);
    EXPECT_NEAR(dLongitude, 0, 1.e-2);
    EXPECT_NEAR(dDepth,     0, 1.e-14);
    EXPECT_NEAR(dTime,      0, 0.2);
    //std::cout << dLat << " " << dLon << " " << dDep << " " << dTime << std::endl;
    //std::cout << std::setprecision(12) << origin.getEpicenter().getLatitude() << " " << origin.getEpicenter().getLongitude() << " " << " " << origin.getDepth() << " " << origin.getTime() << std::endl;

}

TEST(ULocatorNLOptDividedRectangles, BlastLp)
{
// TODO redo this with the original UUSS velocity models
    ULocator::Origin referenceOrigin;
    referenceOrigin.setTime(1646684876.8099976);
    referenceOrigin.setEpicenter(ULocator::Position::WGS84 {40.5151667, -112.1455} );
    referenceOrigin.setDepth(-2000);

    auto arrivals = createQuarryBlastArrivals();
    auto calculatorMap = std::make_unique<ULocator::TravelTimeCalculatorMap> (); 
    for (const auto &arrival : arrivals)
    {   
        std::unique_ptr<ULocator::ITravelTimeCalculator> calculator;
        if (arrival.getPhase() == "P")
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::P, 
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        else
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::S,
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        calculatorMap->insert(arrival.getStation(), arrival.getPhase(),
                              std::move(calculator));
    }
    auto topography = std::make_unique<ULocator::Topography::Constant> ();
    topography->set(2000);

    ULocator::Optimizers::NLOpt::DividedRectangles direct;
    direct.setGeographicRegion(ULocator::Position::UtahRegion {});
    direct.setTravelTimeCalculatorMap(std::move(calculatorMap));
    direct.setTopography(std::move(topography));
    direct.setArrivals(arrivals);
    EXPECT_TRUE(direct.locallyBias());
    EXPECT_TRUE(direct.normalize());
    EXPECT_NEAR(direct.getLocationTolerance(), 1, 1.e-8);
    EXPECT_NEAR(direct.getOriginTimeTolerance(), 0.01, 1.e-8);
    EXPECT_EQ(direct.getMaximumNumberOfObjectiveFunctionEvaluations(), 5000);

    ULocator::Origin noGuess;
    EXPECT_NO_THROW(direct.setMaximumNumberOfObjectiveFunctionEvaluations(2500));

    EXPECT_NO_THROW(
        direct.locateEventAtFreeSurface(
            ULocator::Optimizers::IOptimizer::Norm::LeastSquares )
    );
    EXPECT_TRUE(direct.haveOrigin());
    EXPECT_EQ(direct.getNumberOfGradientEvaluations(), 0);

    auto origin = direct.getOrigin();
    auto dLatitude = referenceOrigin.getEpicenter().getLatitude()
                   - origin.getEpicenter().getLatitude();
    auto dLongitude = referenceOrigin.getEpicenter().getLongitude()
                    - origin.getEpicenter().getLongitude();
    auto dDepth = referenceOrigin.getDepth() - origin.getDepth();
    auto dTime = referenceOrigin.getTime() - origin.getTime();
    EXPECT_NEAR(dLatitude,  0, 1.e-2);
    EXPECT_NEAR(dLongitude, 0, 1.e-2);
    EXPECT_NEAR(dDepth,     0, 1.e-14);
    EXPECT_NEAR(dTime,      0, 0.2);
std::cout << direct.getOptimalObjectiveFunction() << std::endl;
ULocator::Origin newOrigin;
newOrigin.setEpicenter(origin.getEpicenter());
newOrigin.setDepth(origin.getDepth());
std::cout << direct.evaluateObjectiveFunction(newOrigin, ULocator::Optimizers::IOptimizer::Norm::LeastSquares) << std::endl;
std::cout << direct.evaluateObjectiveFunction(newOrigin, ULocator::Optimizers::IOptimizer::Norm::L1) << std::endl;
std::cout << direct.evaluateObjectiveFunction(newOrigin, ULocator::Optimizers::IOptimizer::Norm::Lp) << std::endl;
}

TEST(ULocatorNLOptDividedRectangles, StoGoLp)
{
// TODO redo this with the original UUSS velocity models
    ULocator::Origin referenceOrigin;
    referenceOrigin.setTime(1646479475.880001);
    referenceOrigin.setEpicenter(ULocator::Position::WGS84 {38.5615, -112.2196667} );
    referenceOrigin.setDepth(2000);

    auto arrivals = createEarthquakeArrivals();
    auto calculatorMap = std::make_unique<ULocator::TravelTimeCalculatorMap> ();
    for (const auto &arrival : arrivals)
    {
        std::unique_ptr<ULocator::ITravelTimeCalculator> calculator;
        if (arrival.getPhase() == "P")
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::P,
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        else
        {
            calculator
                = std::make_unique<ULocator::UUSSRayTracer>
                  (arrival.getStation(),
                   ULocator::UUSSRayTracer::Phase::S,
                   ULocator::UUSSRayTracer::Region::Utah);
        }
        calculatorMap->insert(arrival.getStation(), arrival.getPhase(),
                              std::move(calculator));
    }
    auto topography = std::make_unique<ULocator::Topography::Constant> ();
    topography->set(4000);

    ULocator::Optimizers::NLOpt::StochasticGradientOptimization stogo;
    stogo.setGeographicRegion(ULocator::Position::UtahRegion {});
    stogo.setTravelTimeCalculatorMap(std::move(calculatorMap));
    stogo.setTopography(std::move(topography));
    stogo.setArrivals(arrivals);
    EXPECT_TRUE(stogo.randomize());
    EXPECT_NEAR(stogo.getLocationTolerance(), 1, 1.e-8);
    EXPECT_NEAR(stogo.getOriginTimeTolerance(), 0.01, 1.e-8);
    EXPECT_EQ(stogo.getMaximumNumberOfObjectiveFunctionEvaluations(), 5000);

    ULocator::Origin noGuess;
    noGuess.setEpicenter(referenceOrigin.getEpicenter());
    noGuess.setTime(referenceOrigin.getTime());
    noGuess.setDepth(2000);
    stogo.locate(referenceOrigin,
                 ULocator::Optimizers::IOptimizer::LocationProblem::ThreeDimensionsAndTime,
                 ULocator::Optimizers::IOptimizer::Norm::Lp);
    EXPECT_EQ(stogo.getNumberOfObjectiveFunctionEvaluations(), 0);
    EXPECT_TRUE(stogo.getNumberOfGradientEvaluations() > 0);
 
    EXPECT_TRUE(stogo.haveOrigin());
    auto origin = stogo.getOrigin();
    auto dLatitude = referenceOrigin.getEpicenter().getLatitude()
                   - origin.getEpicenter().getLatitude();
    auto dLongitude = referenceOrigin.getEpicenter().getLongitude()
                    - origin.getEpicenter().getLongitude();
    auto dDepth = referenceOrigin.getDepth() - origin.getDepth();
    auto dTime = referenceOrigin.getTime() - origin.getTime();
std::cout << dLatitude << " " << dLongitude << " " << dDepth << " " << dTime << std::endl;


}

