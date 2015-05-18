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
  bool vary_temp;
  bool vary_height;
  bool run_once;
  string out_summary;
  string out_persist;
  string out_phylogeny;
  string out_species_stats;
  int maxiter;
  double mprob;
  double tdrift;
  double simthresh;

  seed_rand(); //TODO: Uncomment this line before production

  if(argc!=11){
    cout<<"Syntax: "<<argv[0]<<" <Summary Stats> <Persistance Graph Output Base> <Phylogeny Output Base> <SpeciesStats Output Base> <(No)VaryHeight> <(No)VaryTemp> <RunOnce/RunMany> <Maximum Iterations> <Mortality Prob (~1e-3)> <Temperature Drift Rate (~0.1)> <Species Similarity Threshold (~0.95)>"<<endl;
    cout<<"'once' indicates the program should be run once, for testing."<<endl;
    cout<<"Names marked base will have file extensions automatically appended."<<endl;
    return -1;
  } else {
    out_summary              = argv[1];
    out_persist              = argv[2];
    out_phylogeny            = argv[3];
    out_species_stats        = argv[4];
    string vary_height_or_no = argv[5];
    string vary_temp_or_no   = argv[6];
    string run_once_or_many  = argv[7];
    maxiter                  = stoi(std::string(argv[8]));
    mprob                    = stod(std::string(argv[9]));
    tdrift                   = stod(std::string(argv[10]));
    simthresh                = stod(std::string(argv[11]));

    if( !(vary_height_or_no=="NoVaryHeight" || vary_height_or_no=="VaryHeight")){
      cerr<<"Unrecognised vary height directive."<<endl;
      return -1;
    }
    if( !(vary_temp_or_no=="NoVaryTemp" || vary_temp_or_no=="VaryTemp")){
      cerr<<"Unrecognised vary temp directive."<<endl;
      return -1;
    }
    if( !(run_once_or_many=="RunOnce" || run_once_or_many=="RunMany")){
      cerr<<"Unrecognised vary temp directive."<<endl;
      return -1;
    }
    vary_height = (vary_height_or_no == "VaryHeight");
    vary_temp   = (vary_temp_or_no   == "VaryTemp"  );
    run_once    = (run_once_or_many  == "RunOnce"   );
  }

  Temperature::getInstance().init("data/temp_series_degreesC_0_65MYA_by_0.001MY.csv");
  if(!vary_temp) {
    Temperature::getInstance().testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE
  }

  if(argc==5 && run_once){
    Simulation sim(0.001, 0.01, 0.96, 1, 0.5, vary_height); //The last argument sets the timestep
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

  //Set up the runs NOTE: OpenMP cannot perform parallel looping with floating-
  //point numbers. So don't try a "pragam omp parallel collapse (5)" here, or
  //some such nonesense.
  //Set up the runs
  for(int iterationnumber=0; iterationnumber<maxiter; iterationnumber++)
  for(double mutation_probability=mprob; mutation_probability<=mprob; mutation_probability+=5e-3)
  for(double temperature_drift_sd=tdrift; temperature_drift_sd<=tdrift; temperature_drift_sd++)
  for(double sim_thresh=simthresh; sim_thresh<=simthresh; sim_thresh++)
  for(double tempdeathfactor=1; tempdeathfactor<=1; tempdeathfactor++){
    Simulation temp(mutation_probability,temperature_drift_sd,sim_thresh,tempdeathfactor,0.5,vary_height); //The last argument sets the timestep
    runs.push_back(temp);
  }

  //Run the runs
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    runs[i].runSimulation();
    //runs[i].dumpPhylogeny();
  }

  //Print out the summary statistics of all of the runs
  cerr<<"Printing summaries to: "<<out_summary<<endl;
  std::ofstream f_summary(out_summary.c_str());
  f_summary<<SimulationSummaryHeader()<<endl;
  for(unsigned int r=0;r<runs.size();++r)
    printSimulationSummary(f_summary, r, runs[r]);

  //Print out the phylogenies and persistance graphs of the runs which approximate
  //the phylogeny of Kozak and Wiens (2010)
  for(unsigned int i=0;i<runs.size();++i){
    //Output persistance table for each run within the boundries
    string outputname_persist=std::string(out_persist)+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_persist(outputname_persist);
    runs[i].phylos.persistGraph(f_persist);

    //Output phylogeny for each run within the boundries
    string outputname_phylo=std::string(out_phylogeny)+"_run_"+std::to_string(i)+".tre";
    std::ofstream f_phylogeny(outputname_phylo);
    f_phylogeny   <<runs[i].phylos.printNewick() <<endl;

    string outputname_species_stats = std::string(out_species_stats)+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_species_stats(outputname_species_stats);
    runs[i].phylos.speciesSummaries(f_species_stats);
  }

  return 0;
}
