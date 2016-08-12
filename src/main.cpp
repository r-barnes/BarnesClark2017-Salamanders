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
#include "omp.h"
using namespace std;

std::string SimulationSummaryHeader() {
  return "RunNum, MutationProb, TempDriftSD, SimThresh, Nspecies, ECDF, "
         "AvgOtempdegC, Nalive, EndTime, AvgElevation";
}

void printSimulationSummary(ofstream &out, int r, const Simulation &sim){
  out<<r;
  out<<", " << TheParams.mutationProb();
  out<<", " << TheParams.tempDrift();
  out<<", " << TheParams.speciesSimthresh();
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
    cout<<"\tDebug                     YES/NO      Print real-time stats about simulation.\n";
    cout<<"\tVaryHeight                YES/NO                 \n";
    cout<<"\tVaryTemp                  YES/NO                 \n";
    cout<<"\tRunOnce                   YES/NO                 \n";
    cout<<"\tNumBins                   Integer     Number of elevation bins which comprise the mountain.\n";
    cout<<"\tMutationProb              Double                 \n";
    cout<<"\tTemperatureDrift          Double                 \n";
    cout<<"\tSpeciesSimilarity         Integer                \n";
    cout<<"\ttimestep                  Double      Units are in My.\n";
    cout<<"\tDispersalProb             Double      Prob of salamander trying to move between bins.\n";
    cout<<"\tDispersalType             String      Must be: Global, Better, MaybeWorse\n";
    cout<<"\tmaxiter                   Integer     Numer of times to run each parameter combination.\n";
    cout<<"\tPRNGseed                  Integer     Use 0 for true entropy; otherwise, seed engine to given value.\n";
    cout<<"\tTempSeries                Filename    If not specified, a default time series is used. NO SPACES ALLOWED!\n";
    cout<<"\tInitialAltitude           Integer     Altitude of progenitor species. Must be a valid bin number.\n";
    cout<<"\tInitialPopSize            Integer     How many salamanders to put in the initial bin.\n";
    cout<<"\tLogitTempWeight           Double      Adjusts how salamander mortality relates to temperature.\n";
    cout<<"\tLogitOffset               Double      \n";
    cout<<"\tLogitCAweight             Double      Conspecific abundance weighting in mortality.\n";
    cout<<"\tLogitHAweight             Double      Heterospecific abundance weighting in mortality.\n";
    cout<<"\tMaxOffspringPerBinPerDt   Integer     \n";
    cout<<"\tMaxTriesToBreed           Integer     \n";
    cout<<"\tToLowlandsProb            Double      \n";
    cout<<"\tFromLowlandsProb          Double      \n";

    return -1;
  }

  TheParams.load(argv[1]);

  seed_rand(TheParams.randomSeed());

  Temperature.init(TheParams.tempSeriesFilename());
  if(!TheParams.pVaryTemp()) {
    Temperature.testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE
  }

  if(TheParams.pRunOnce()){
    std::cerr<<"Feature deprecated!"<<std::endl;
    // Simulation sim();
    // sim.runSimulation();

    // std::ofstream f_persist  ((TheParams.outPersistFilename()+".csv"));
    // std::ofstream f_phylogeny((TheParams.outPhylogenyFilename()+".tre"));
    // std::ofstream f_summary  (TheParams.outSummaryFilename());
    // f_summary<<SimulationSummaryHeader()<<endl;
    // printSimulationSummary(f_summary, 0, sim);
    // sim.phylos.persistGraph(f_persist);
    // f_phylogeny<<sim.phylos.printNewick()<<endl;
    return -1;
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
  for(int i=0;i<TheParams.maxiter();i++)
    runs.emplace_back();

  //Used to show more detailed, real-time info about simulation
  if(TheParams.debug()){
    omp_set_num_threads(1);
    runs.clear();
    runs.emplace_back();
  }

  //Run the simulations in parallel using OpenMP
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    //#pragma omp critical
    //  cout<<"Run #"<<i<<endl;
    runs[i].runSimulation();
    //runs[i].dumpPhylogeny();
  }

  //cerr<<"Printing output information..."<<endl;
  
  //Print out the summary statistics of all of the runs
  std::ofstream f_summary(TheParams.outSummaryFilename());
  f_summary<<SimulationSummaryHeader()<<endl;
  for(unsigned int r=0;r<runs.size();++r)
    printSimulationSummary(f_summary, r, runs[r]);

  //Output persistence table for each run within the boundaries
  {
    std::ofstream f_persist(TheParams.outPersistFilename());
    for(unsigned int i=0;i<runs.size();i++)
      runs[i].phylos.persistGraph(i, f_persist);
  }

  //Output phylogeny for each run within the boundaries
  {
    std::ofstream f_phylogeny(TheParams.outPhylogenyFilename());
    for(unsigned int i=0;i<runs.size();i++)
      f_phylogeny<<i<<" "<<runs[i].phylos.printNewick()<<endl;
  }

  //Output summaries of the distribution of species properties at each point
  //in time
  {
    std::ofstream f_species_stats(TheParams.outSpeciesStatsFilename());
    for(unsigned int i=0;i<runs.size();i++)
      runs[i].phylos.speciesSummaries(i, f_species_stats);
  }

  return 0;
}