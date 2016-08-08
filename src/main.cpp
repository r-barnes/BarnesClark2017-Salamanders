//This file reads in information used to set the parameters of a simulation or
//ensemble of simulations, builds a list of simulations, executes them in
//parallel, and prints the results.
#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include "simulation.hpp"
#include "temp.hpp"
#include "random.hpp"
#include "params.hpp"
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

int main(int argc, char **argv){
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
    cout<<"\tDispersalProb             Double      Prob of salamander trying to move between bins.\n";
    cout<<"\tDispersalType             String      Must be: Global, Better, MaybeWorse\n";
    cout<<"\tmaxiter                   Integer     Numer of times to run each parameter combination.\n";
    cout<<"\tPRNGseed                  Integer     Use 0 for true entropy; otherwise, seed engine to given value.\n";
    cout<<"\tTempSeries                Filename    If not specified, a default time series is used. NO SPACES ALLOWED!\n";
    cout<<"\tInitialAltitude           Integer     Altitude of progenitor species. Must be a valid bin number.\n";
    return -1;
  }

  Params params(argv[1]);

  seed_rand(params.getRandomSeed());

  //TODO: Incorporate dispersal_prob, dispersal_type, initial_altitude

  Temperature::getInstance().init(params.getTempSeriesFilename());
  if(!params.pVaryTemp()) {
    Temperature::getInstance().testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE
  }

  if(params.pRunOnce()){
    Simulation sim(0.001, 0.01, 0.96, 1, params.getTimestep(), params.pVaryHeight());
    sim.runSimulation();

    std::ofstream f_persist  ((params.getOutPersistFilename()+".csv"));
    std::ofstream f_phylogeny((params.getOutPhylogenyFilename()+".tre"));
    std::ofstream f_summary  (params.getOutSummaryFilename());
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
  for(int iterationnumber=0; iterationnumber<params.getMaxiter(); iterationnumber++)
  for(double mutation_probability=params.getMutationProb(); mutation_probability<=params.getMutationProb(); mutation_probability+=5e-3)
  for(double temperature_drift_sd=params.getTempDrift(); temperature_drift_sd<=params.getTempDrift(); temperature_drift_sd++)
  for(double sim_thresh=params.getSimthresh(); sim_thresh<=params.getSimthresh(); sim_thresh++)
  for(double tempdeathfactor=1; tempdeathfactor<=1; tempdeathfactor++){
    Simulation temp(
      mutation_probability,
      temperature_drift_sd,
      sim_thresh,
      tempdeathfactor,
      params.getTimestep(),
      params.pVaryHeight()
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
  cerr<<"Printing summaries to: "<<params.getOutSummaryFilename()<<endl;
  std::ofstream f_summary(params.getOutSummaryFilename());
  f_summary<<SimulationSummaryHeader()<<endl;
  for(unsigned int r=0;r<runs.size();++r)
    printSimulationSummary(f_summary, r, runs[r]);

  //Print out the phylogenies and persistence graphs of the runs which approximate
  //the phylogeny of Kozak and Wiens (2010)
  for(unsigned int i=0;i<runs.size();++i){
    //Output persistence table for each run within the boundaries
    string fname_persist=params.getOutPersistFilename()+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_persist(fname_persist);
    runs[i].phylos.persistGraph(f_persist);

    //Output phylogeny for each run within the boundaries
    string fname_phylo=params.getOutPhylogenyFilename()+"_run_"+std::to_string(i)+".tre";
    std::ofstream f_phylogeny(fname_phylo);
    f_phylogeny   <<runs[i].phylos.printNewick() <<endl;

    //Output summaries of the distribution of species properties at each point
    //in time
    string fname_species_stats = params.getOutSpeciesStatsFilename()+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_species_stats(fname_species_stats);
    runs[i].phylos.speciesSummaries(f_species_stats);
  }

  return 0;
}