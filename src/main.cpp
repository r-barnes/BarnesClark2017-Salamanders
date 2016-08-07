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
#include <stdexcept>
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

void Input_CheckParamName(std::ifstream &fparam, const std::string &param_name){
  std::string in_param_name;
  fparam>>in_param_name;
  if(in_param_name!=param_name){
    std::cerr<<"Unexpected parameter name! Found '"<<in_param_name<<"' expected '"<<param_name<<"'"<<endl;
    throw std::runtime_error("Unexpected parameter name!");
  }
}

std::string Input_Filename(std::ifstream &fparam, const std::string &param_name){
  Input_CheckParamName(fparam, param_name);
  std::string temp;
  fparam>>temp;
  return temp;
}

bool Input_YesNo(std::ifstream &fparam, const std::string &param_name){
  Input_CheckParamName(fparam, param_name);
  std::string temp;
  fparam>>temp;
  if(! (temp=="YES" || temp=="NO")){
    cerr<<"Paremeter '"<<param_name<<"' must be YES or NO!"<<endl;
    throw std::runtime_error("Bad parameter value!");
  }
  return (temp=="YES");
}

double Input_Double(std::ifstream &fparam, const std::string &param_name){
  Input_CheckParamName(fparam, param_name);
  double temp;
  fparam>>temp;
  return temp;
}

int Input_Integer(std::ifstream &fparam, const std::string &param_name){
  Input_CheckParamName(fparam, param_name);
  int temp;
  fparam>>temp;
  return temp;
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

  std::string temp_series_filename = "data/temp_series_degreesC_0_65MYA_by_0.001MY.csv";

  //The simulation timestep.
  double timestep = 0.5; //My

  if(argc!=2){
    cout<<"Syntax: "<<argv[0]<<" <Parameters File>\n";
    cout<<"Parmeters file can contain:\n";
    cout<<"\tPARAM                     TYPE        DESCRIPTION\n";
    cout<<"\tSummaryStatsFilename      Filename               \n";
    cout<<"\tPersistenceGraphFilename  Filename               \n";
    cout<<"\tPhylogenyFilename         Filename               \n";
    cout<<"\tSpeciesStatsFilename      Filename               \n";
    cout<<"\tVaryHeight                YES/NO                 \n";
    cout<<"\tVaryTemp                  YES/NO                 \n";
    cout<<"\tRunOnce                   YES/NO                 \n";
    cout<<"\tMutationProb              Double                 \n";
    cout<<"\tTemperatureDrift          Double                 \n";
    cout<<"\tSpeciesSimilarity         Double                 \n";
    cout<<"\ttimestep                  Double      Units are in My.\n";
    cout<<"\tmaxiter                   Integer     Numer of times to run each parameter combination.\n";
    cout<<"\tPRNGseed                  Integer     Use 0 for true entropy; otherwise, seed engine to given value.\n";
    cout<<"\tTempSeries                Filename    If not specified, a default time series is used. NO SPACES ALLOWED!\n";
    return -1;
  }

  std::ifstream fparam(argv[1]);
  out_summary       = Input_Filename(fparam,"SummaryStatsFilename");
  out_persist       = Input_Filename(fparam,"PersistenceGraphFilename");
  out_phylogeny     = Input_Filename(fparam,"PhylogenyFilename");
  out_species_stats = Input_Filename(fparam,"SpeciesStatsFilename");

  vary_height = Input_YesNo(fparam,"VaryHeight");
  vary_temp   = Input_YesNo(fparam,"VaryTemp");
  run_once    = Input_YesNo(fparam,"RunOnce");

  mprob     = Input_Double(fparam, "MutationProb");
  tdrift    = Input_Double(fparam, "TemperatureDrift");
  simthresh = Input_Double(fparam, "SpeciesSimilarity");
  timestep  = Input_Double(fparam, "timestep");

  maxiter = Input_Integer(fparam,"maxiter");
  seed_rand(Input_Integer(fparam,"PRNGseed"));

  temp_series_filename = Input_Filename(fparam,"TempSeries");

  {
    std::string test_if_file_is_empty;
    if(fparam>>test_if_file_is_empty){
      std::cerr<<"File contained unexpected parameters at the end!"<<std::endl;
      return -1;
    }
  }

  Temperature::getInstance().init(temp_series_filename);
  if(!vary_temp) {
    Temperature::getInstance().testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE
  }

  if(run_once){
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