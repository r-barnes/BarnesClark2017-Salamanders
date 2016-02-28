//This file reads in information used to set the parameters of a simulation or
//ensemble of simulations, builds a list of simulations, executes them in
//parallel, and prints the results.
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
  return "Run #, MutationProb, TempDriftSD, SimThresh, Nspecies, ECDF, "
         "AvgOtempdegC, Nalive, EndTime, AvgElevation";
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

  //The simulation timestep. This choice is somewhat arbitrary.
  const double timestep = 0.5; //My

  //NOTE: This line should be uncommented before running simulations for
  //purposes other than testing.
  seed_rand();

  if(argc!=12){
    cout<<"Syntax: "<<argv[0]<<" <Summary Stats> <Persistence Graph Output Base> ";
    cout<<"<Phylogeny Output Base> <SpeciesStats Output Base> <(No)VaryHeight> ";
    cout<<"<(No)VaryTemp> <RunOnce/RunMany> <Maximum Iterations> <Mutation Prob (~1e-3)> ";
    cout<<"<Temperature Drift Rate (~0.1)> <Species Similarity Threshold (~0.95)>"<<endl;
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
    Simulation sim(0.001, 0.01, 0.96, 1, timestep, vary_height);
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
  //Generate a list of simulations to run

  //NOTE: These are set right now to not change the any of the parameters - this
  //means that the output is a series of 'maxiter' runs all of which have
  //identical parameters, but will vary because of random factors. The bounds of
  //the loops below could be altered to sweep the parameter space.
  for(int iterationnumber=0; iterationnumber<maxiter; iterationnumber++)
  for(double mutation_probability=mprob; mutation_probability<=mprob; mutation_probability+=5e-3)
  for(double temperature_drift_sd=tdrift; temperature_drift_sd<=tdrift; temperature_drift_sd++)
  for(double sim_thresh=simthresh; sim_thresh<=simthresh; sim_thresh++)
  for(double tempdeathfactor=1; tempdeathfactor<=1; tempdeathfactor++){
    Simulation temp(
      mutation_probability,
      temperature_drift_sd,
      sim_thresh,
      tempdeathfactor,
      timestep,
      vary_height
    );
    runs.push_back(temp);
  }

  //Run the simulations in parallel using OpenMP
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

  //Print out the phylogenies and persistence graphs of the runs which approximate
  //the phylogeny of Kozak and Wiens (2010)
  for(unsigned int i=0;i<runs.size();++i){
    //Output persistence table for each run within the boundaries
    string fname_persist=std::string(out_persist)+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_persist(fname_persist);
    runs[i].phylos.persistGraph(f_persist);

    //Output phylogeny for each run within the boundaries
    string fname_phylo=std::string(out_phylogeny)+"_run_"+std::to_string(i)+".tre";
    std::ofstream f_phylogeny(fname_phylo);
    f_phylogeny   <<runs[i].phylos.printNewick() <<endl;

    //Output summaries of the distribution of species properties at each point
    //in time
    string fname_species_stats = std::string(out_species_stats)+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_species_stats(fname_species_stats);
    runs[i].phylos.speciesSummaries(f_species_stats);
  }

  return 0;
}
