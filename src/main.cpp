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
#include "timer.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include "omp.h"
using namespace std;

string SimulationSummaryHeader() {
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
  Timer timer_overall, timer_calc, timer_io;

  timer_overall.start();

  if(argc!=2){
    cout<<"Syntax: "<<argv[0]<<" <Parameters File>\n";
    cout<<"Parmeters file can contain:\n";
    cout<<"\tPARAM                     TYPE        DESCRIPTION\n";
    cout<<"\tSummaryStatsFilename      Filename               \n";
    cout<<"\tPersistenceGraphFilename  Filename               \n";
    cout<<"\tPhylogenyFilename         Filename               \n";
    cout<<"\tSpeciesStatsFilename      Filename               \n";
    cout<<"\tDebug                     YES/NO      ";
      cout<<"Print real-time stats about simulation.\n";
    cout<<"\tVaryHeight                YES/NO                 \n";
    cout<<"\tVaryTemp                  YES/NO                 \n";
    cout<<"\tRunOnce                   YES/NO                 \n";
    cout<<"\tNumBins                   Integer     ";
      cout<<"Number of elevation bins which comprise the mountain.\n";
    cout<<"\tMutationProb              Double                 \n";
    cout<<"\tTemperatureDrift          Double                 \n";
    cout<<"\tSpeciesSimilarity         Integer                \n";
    cout<<"\ttimestep                  Double      ";
      cout<<"Units are in My.\n";
    cout<<"\tDispersalProb             Double      ";
      cout<<"Prob of salamander trying to move between bins.\n";
    cout<<"\tDispersalType             String      ";
      cout<<"Must be: Global, Better, MaybeWorse\n";
    cout<<"\tmaxiter                   Integer     ";
      cout<<"Numer of times to run each parameter combination.\n";
    cout<<"\tPRNGseed                  Integer     ";
      cout<<"Use 0 for true entropy; otherwise, seed engine to given value.\n";
    cout<<"\tTempSeries                Filename    ";
      cout<<"If not specified, a default time series is used. NO SPACES ALLOWED!\n";
    cout<<"\tInitialAltitude           Integer     ";
      cout<<"Altitude of progenitor species. Must be a valid bin number.\n";
    cout<<"\tInitialPopSize            Integer     ";
      cout<<"How many salamanders to put in the initial bin.\n";
    cout<<"\tLogitTempWeight           Double      ";
      cout<<"Adjusts how salamander mortality relates to temperature.\n";
    cout<<"\tLogitOffset               Double      \n";
    cout<<"\tLogitCAweight             Double      ";
      cout<<"Conspecific abundance weighting in mortality.\n";
    cout<<"\tLogitHAweight             Double      ";
      cout<<"Heterospecific abundance weighting in mortality.\n";
    cout<<"\tMaxOffspringPerBinPerDt   Integer     \n";
    cout<<"\tMaxTriesToBreed           Integer     \n";
    cout<<"\tToLowlandsProb            Double      \n";
    cout<<"\tFromLowlandsProb          Double      \n";

    return -1;
  }

  TheParams.load(argv[1]);

  //Seed random number generator for each thread
  timer_calc.start();
  seed_rand(TheParams.randomSeed());
  timer_calc.stop();

  //Load temperature data into global object
  timer_io.start();
  Temperature.init(TheParams.tempSeriesFilename());
  timer_io.stop();

  //If the temperature is set not to vary, then set it to 34 here, which was the
  //sea level temperature 65Mya. It doesn't really matter, though, because Eve
  //is initialized with an optimum temperature equal to that of her starting
  //bin. Since the temperature never changes, the absolute values are therefore
  //unimportant.
  if(!TheParams.pVaryTemp())
    Temperature.testOn(34); //km CHANGE ADDED TO TEST NO TEMP CHANGE

  //TODO: Cut
  if(TheParams.pRunOnce()){
    cerr<<"Feature deprecated!"<<std::endl;
    return -1;
  }

  //Create a vector to hold the runs
  vector<Simulation> runs;

  //Load a number of runs into the vector. Each run has the same parameters, but
  //the runs will differ due to random factors.
  for(int i=0;i<TheParams.maxiter();i++)
    runs.emplace_back();

  //Used to show more detailed, real-time info about simulation
  if(TheParams.debug()){
    omp_set_num_threads(1);
    runs.clear();
    runs.emplace_back();
  }

  //Run the simulations in parallel using OpenMP
  timer_calc.start();
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    //#pragma omp critical
    //  cout<<"Run #"<<i<<endl;
    runs[i].runSimulation();
    //runs[i].dumpPhylogeny();
  }
  timer_calc.stop();

  //cerr<<"Printing output information..."<<endl;
  
  timer_io.start();
  //Print out the summary statistics of all of the runs
  ofstream f_summary(TheParams.outSummaryFilename());
  f_summary<<SimulationSummaryHeader()<<endl;
  for(unsigned int r=0;r<runs.size();++r)
    printSimulationSummary(f_summary, r, runs[r]);

  //Output persistence table for each run within the boundaries
  {
    ofstream f_persist(TheParams.outPersistFilename());
    for(unsigned int i=0;i<runs.size();i++)
      runs[i].phylos.persistGraph(i, f_persist);
  }

  //Output phylogeny for each run within the boundaries
  {
    ofstream f_phylogeny(TheParams.outPhylogenyFilename());
    for(unsigned int i=0;i<runs.size();i++)
      f_phylogeny<<i<<" "<<runs[i].phylos.printNewick()<<endl;
  }

  //Output summaries of the distribution of species properties at each point
  //in time
  {
    ofstream f_species_stats(TheParams.outSpeciesStatsFilename());
    for(unsigned int i=0;i<runs.size();i++)
      runs[i].phylos.speciesSummaries(i, f_species_stats);
  }
  timer_io.stop();
  timer_overall.stop();

  cerr<<"Time overall: "<<timer_overall.accumulated() <<endl;
  cerr<<"Time calc:    "<<timer_calc.accumulated()    <<endl;
  cerr<<"Time IO:      "<<timer_io.accumulated()      <<endl;

  return 0;
}