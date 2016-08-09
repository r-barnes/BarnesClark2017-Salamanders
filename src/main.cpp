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
  out<<", " << TheParams::get().mutationProb();
  out<<", " << TheParams::get().tempDrift();
  out<<", " << TheParams::get().speciesSimthresh();
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
    cout<<"\tNumBins                   Integer     Number of elevation bins which comprise the mountain.\n";
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
    cout<<"\tLogitTempWeight           Double      Adjusts how salamander mortality relates to temperature.\n";
    cout<<"\tLogitOffset               Double      \n";
    cout<<"\tLogitCAweight             Double      Conspecific abundance weighting in mortality.\n";
    cout<<"\tLogitHAweight             Double      Heterospecific abundance weighting in mortality.\n";
    cout<<"\tMaxOffspringPerBinPerDt   Integer     \n";
    cout<<"\tMaxTriesToBreed           Integer     \n";

    return -1;
  }

  TheParams::get().load(argv[1]);

  seed_rand(TheParams::get().randomSeed());

  Temperature::getInstance().init(TheParams::get().tempSeriesFilename());
  if(!TheParams::get().pVaryTemp()) {
    Temperature::getInstance().testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE
  }

  if(TheParams::get().pRunOnce()){
    std::cerr<<"Feature deprecated!"<<std::endl;
    // Simulation sim();
    // sim.runSimulation();

    // std::ofstream f_persist  ((TheParams::get().outPersistFilename()+".csv"));
    // std::ofstream f_phylogeny((TheParams::get().outPhylogenyFilename()+".tre"));
    // std::ofstream f_summary  (TheParams::get().outSummaryFilename());
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
  for(int i=0;i<TheParams::get().maxiter();i++)
    runs.emplace_back();

  //Run the simulations in parallel using OpenMP
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    runs[i].runSimulation();
    //runs[i].dumpPhylogeny();
  }

  //Print out the summary statistics of all of the runs
  cerr<<"Printing summaries to: "<<TheParams::get().outSummaryFilename()<<endl;
  std::ofstream f_summary(TheParams::get().outSummaryFilename());
  f_summary<<SimulationSummaryHeader()<<endl;
  for(unsigned int r=0;r<runs.size();++r)
    printSimulationSummary(f_summary, r, runs[r]);

  //Print out the phylogenies and persistence graphs of the runs which approximate
  //the phylogeny of Kozak and Wiens (2010)
  for(unsigned int i=0;i<runs.size();++i){
    //Output persistence table for each run within the boundaries
    string fname_persist=TheParams::get().outPersistFilename()+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_persist(fname_persist);
    runs[i].phylos.persistGraph(f_persist);

    //Output phylogeny for each run within the boundaries
    string fname_phylo=TheParams::get().outPhylogenyFilename()+"_run_"+std::to_string(i)+".tre";
    std::ofstream f_phylogeny(fname_phylo);
    f_phylogeny   <<runs[i].phylos.printNewick() <<endl;

    //Output summaries of the distribution of species properties at each point
    //in time
    string fname_species_stats = TheParams::get().outSpeciesStatsFilename()+"_run_"+std::to_string(i)+".csv";
    std::ofstream f_species_stats(fname_species_stats);
    runs[i].phylos.speciesSummaries(f_species_stats);
  }

  return 0;
}