#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include "simulation.hpp"
#include "temp.hpp"
#include "random.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;

std::string SimulationSummaryHeader() {
  return "Run #, MutationProb, TempDriftSD, SimThresh, Nspecies, ECDF, AvgOtempdegC, Nalive, EndTime, AvgElevation";
}

void printSimulationSummary(ofstream &out, int r, const Simulation &sim){
  out<<r;
  out<<", " << sim.mutation_probability;
  out<<", " << sim.temperature_drift_sd;
  out<<", " << sim.species_sim_thresh;
  out<<", " << sim.nspecies;
  out<<", " << sim.ecdf;
  out<<", " << sim.avg_otempdegC;
  out<<", " << sim.salive;
  out<<", " << sim.endtime;
  out<<", " << sim.avg_elevation;
  out<<endl;
}

int main(int argc, char **argv){
  if(argc<4 || argc>5){
    cout<<"Syntax: "<<argv[0]<<" <Summary Stats> <Persistance Graph Output Base> <Phylogeny Output Base> [once]"<<endl;
    cout<<"'once' indicates the program should be run once, for testing."<<endl;
    cout<<"Names marked base will have file extensions automatically appended."<<endl;
    return -1;
  }

  Temperature::getInstance().init("data/temp_series_degreesC_0_65MYA_by_0.001MY.csv");
  if(!vary_temp) {
    Temperature::getInstance().testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE
  }

  string out_summary   = argv[1];
  string out_persist   = argv[2];
  string out_phylogeny = argv[3];

  //seed_rand(); //TODO: Uncomment this line before production

  if(argc==5 && std::string(argv[4])==std::string("once")){
    Simulation sim(0.001, 0.01, 0.96, 1, 0.5); //The last argument sets the timestep
    sim.runSimulation();

    out_persist   += ".csv";
    out_phylogeny += ".tre";

    std::ofstream f_persist  (out_persist.c_str());
    std::ofstream f_phylogeny(out_phylogeny.c_str());
    std::ofstream f_summary  (out_summary.c_str());
    f_summary<<SimulationSummaryHeader()<<endl;
    printSimulationSummary(f_summary, 0, sim);
    sim.phylos.persistGraph(f_persist);
    f_phylogeny<<sim.phylos.printNewick()<<endl;
    return 0;
  }

  //Create a vector to hold the runs
  std::vector<Simulation> runs;

  //Set up the runs
  for(int iterationnumber=0; iterationnumber<1; iterationnumber++)
  for(double mutation_probability=1e-3; mutation_probability<1e-1; mutation_probability+=5e-3)
  for(double temperature_drift_sd=0.1; temperature_drift_sd<10; temperature_drift_sd+=5e-1)
  for(double sim_thresh=0.90; sim_thresh<1; sim_thresh+=0.02)
  for(double tempdeathfactor=1; tempdeathfactor<=1; tempdeathfactor+=0.1){
    Simulation temp(mutation_probability,temperature_drift_sd,sim_thresh,tempdeathfactor,0.5); //The last argument sets the timestep
    runs.push_back(temp);
  }

  //Run the runs
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    runs[i].runSimulation();
  }

  //Print out the summary statistics of all of the runs
  std::ofstream f_summary(out_summary.c_str());
  f_summary<<SimulationSummaryHeader()<<endl;
  for(unsigned int r=0;r<runs.size();++r)
    printSimulationSummary(f_summary, r, runs[r]);

  //Print out the phylogenies and persistance graphs of the runs which approximate
  //the phylogeny of Kozak and Wiens (2010)
  //for(unsigned int i=0;i<runs.size();++i){
    //Output persistance table for each run within the boundries
    //string outputname_persist=std::string(out_persist)+"_run_"+std::to_string(i)+".csv";
    //std::ofstream f_persist(outputname_persist);
    //runs[i].phylos.persistGraph(f_persist);

    //Output phylogeny for each run within the boundries
    //string outputname_phylo=std::string(out_phylogeny)+"_run_"+std::to_string(i)+".tre";
    //std::ofstream f_phylogeny(outputname_phylo);
    //f_phylogeny   <<runs[i].phylos.printNewick() <<endl;
  //}

  return 0;
}
