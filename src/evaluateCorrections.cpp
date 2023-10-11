#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <boost/tokenizer.hpp>
#include "uLocator/sourceSpecificStationCorrection.hpp"
#include "uLocator/staticCorrection.hpp"

using namespace ULocator;

struct Data
{
    std::string network;
    std::string station;
    std::string phase;
    double arrivalTime;
    double uncorrectedTime;
    double residual;
    double sourceLat;
    double sourceLon;
    double sourceDep;
};

struct NSP
{
    NSP(const std::string &n, const std::string &s, const std::string &p) :
       network(n), station(s), phase(p)
    {
    }
    std::string network;
    std::string station;
    std::string phase;
};

std::vector<Data> loadCatalog(const std::string &fileName)
{
    std::vector<Data> result;
    std::string line;
    std::ifstream csvFile(fileName);
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
        //event_identifier,latitude,longitude,depth,origin_time,network,station,phase,arrival_time,standard_error,residual,uncorrected_travel_time,source_receiver_distance,event_type,n_objective_function_evaluations,weightedRMS
        //60000071,39.79217290175797,-113.7219946522614,6400,1349058838.948746,UU,NOQ,P,1349058865.790473,0.03,0.04147601127624512,26.80025029182434,166461.9340866167,eq,206,1.004027800049852
        Data data;
        data.sourceLat = std::stod(splitLine[1]);
        data.sourceLon = std::stod(splitLine[2]);
        data.sourceDep = std::stod(splitLine[3]);
        data.network = splitLine[5];
        data.station = splitLine[6];
        data.phase = splitLine[7];
        data.arrivalTime = std::stod(splitLine[8]);
        data.residual = std::stod(splitLine[10]);
        data.uncorrectedTime = std::stod(splitLine[4]) + std::stod(splitLine[11]);
        result.push_back(data);
    }
    csvFile.close();
    return result;
}


int main(int argc, char *argv[])
{
    std::string catalogFile{"../examples/uuss/utahResultsL2/utahTestingCatalog.csv"};
    std::string correctionsFile{"../examples/uuss/utahResultsL2/utahCorrections.3.h5"};
    auto arrivals = loadCatalog(catalogFile);
    std::vector<NSP> uniqueNSPs;
    std::map<std::string, std::unique_ptr<StaticCorrection>> statics;
    std::map<std::string, std::unique_ptr<SourceSpecificStationCorrection>> ssscs;
    for (const auto &a : arrivals)
    {
        bool found = false;
        for (const auto &nsp : uniqueNSPs)
        {
            if (nsp.network == a.network  &&
                nsp.station == a.station &&
                nsp.phase   == a.phase)
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            NSP nsp(a.network, a.station, a.phase);
            auto name = a.network + "." + a.station + "." + a.phase;
            uniqueNSPs.push_back(nsp);
            try
            {
                auto staticCorrection = std::make_unique<StaticCorrection> ();
                staticCorrection->initialize(a.network, a.station, a.phase, 0.0);
                staticCorrection->load(correctionsFile);
                statics.insert(std::pair {name, std::move(staticCorrection)});
            }
            catch (const std::exception &e)
            {

            }
            try
            {
                auto sssc = std::make_unique<SourceSpecificStationCorrection> ();
                sssc->initialize(a.network, a.station, a.phase);
                sssc->load(correctionsFile);
                ssscs.insert(std::pair {name, std::move(sssc)});
            }
            catch (const std::exception &e)
            {

            }
        }
    }
    std::ofstream ofl;
    ofl.open("trainingCorrectionResiduals.csv");
    ofl << "uncorrected_residual,static_residual,sssc_residual,static_sssc_residual" << std::endl;
    for (const auto &a : arrivals)
    {
if (a.phase == "P"){continue;}
        auto name = a.network + "." + a.station + "." + a.phase;
        double c1 = 0;
        double c2 = 0;
        double residual = a.arrivalTime - a.uncorrectedTime;
        double res1 = residual;
        double sgn =+1;
        if (statics.contains(name))
        {
            res1 = a.arrivalTime - (a.uncorrectedTime + sgn*statics[name]->evaluate(0));
        }
        else
        {
            continue;
        } 
        double res2 = residual;
        double res3 = residual;
        if (ssscs.contains(name))
        {
            double t1 = statics[name]->evaluate(0);
            res2 = a.arrivalTime - (a.uncorrectedTime + sgn*ssscs[name]->evaluate(a.sourceLat, a.sourceLon, a.sourceDep));
            res3 = a.arrivalTime - (a.uncorrectedTime + sgn*(t1 + ssscs[name]->evaluate(a.sourceLat, a.sourceLon, a.sourceDep)));
        }
        else
        {
            continue;
        }
        ofl << residual << "," << res1 << "," << res2 << "," << res3 << std::endl;
    }
    ofl.close();
    return EXIT_SUCCESS;
}
