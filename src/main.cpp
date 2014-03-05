#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include "simulation.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
using namespace std;

int main(int argc, char **argv){
  if(argc<3 || argc>4){
    cout<<"Syntax: "<<argv[0]<<" <Persistance Graph Output> <Phylogeny Output> [once]"<<endl;
    cout<<"'once' indicates the program should be run once, for testing."<<endl;
    return -1;
  }

  //srand (time(NULL)); //TODO: Uncomment this line before production

  if(argc==4 && std::string(argv[3])==std::string("once")){
    std::ofstream out_persistgraph(argv[1]);
    std::ofstream out_phylogeny   (argv[2]);
    Simulation sim(0.001, 0.01, 0.96);
    sim.runSimulation();
    sim.phylos.persistGraph(out_persistgraph);
    out_phylogeny <<sim.phylos.printNewick() <<endl;
    cout<<"Avg temp: "                <<sim.avg_otempdegC<<endl;
    cout<<"Alive salamanders at end: "<<sim.salive       <<endl;
    return 0;
  }

  //Create a vector to hold the runs
  std::vector<Simulation> runs;

  //Set up the runs
  for(double mutation_probability=1e-4; mutation_probability<1e-3; mutation_probability+=2e-4)
  for(double temperature_drift_sd=1; temperature_drift_sd<10; temperature_drift_sd+=2)
  for(double sim_thresh=0.95; sim_thresh<1; sim_thresh+=0.01){
    Simulation temp(mutation_probability,temperature_drift_sd,sim_thresh);
    runs.push_back(temp);
  }

  //Run the runs and stash the results.
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    runs[i].runSimulation();
  }

  //Print out the final parameters of the runs
  cout<<"Run #, MutationProb, TempDriftSD, SimThresh, Nspecies, ECDF, AvgOtempdegC, Nalive"<<endl;
  for(unsigned int r=0;r<runs.size();++r){
    cout<<r<<",";
    cout<<", " << runs[r].mutation_probability;
    cout<<", " << runs[r].temperature_drift_sd;
    cout<<", " << runs[r].species_sim_thresh;
    cout<<", " << runs[r].nspecies;
    cout<<", " << runs[r].ecdf;
    cout<<", " << runs[r].avg_otempdegC;
    cout<<", " << runs[r].salive;
    cout<<endl;
  }

  int min=0;
  for(unsigned int i=1;i<runs.size();++i)
    //Mark as best match if there are 80-120 species alive at the end of the
    //simulation, and if the ECDF is more similar to the Kozak and Wiens data
    //based on ECDF of mean branch distances than previous results.
    //otemp limits are from Gifford (&Martin) 1970: Mean +/- 2SD for salamander
    //feeding optimum
    if(runs[i].ecdf<runs[min].ecdf && (80<runs[i].nspecies && runs[i].nspecies<120) &&
         (4.629879<runs[i].avg_otempdegC && runs[i].avg_otempdegC<20.97868)) {
      min=i;

      Phylogeny currentphylos=runs[i].phylos;

      //Output persistance table for each run within the boundries
      string outputname_persist=std::string(argv[1])+"_run_"+std::to_string(i);
      std::ofstream out_persistgraph(outputname_persist);
      currentphylos.persistGraph(out_persistgraph);

      //Output phylogeny for each run within the boundries
      string outputname_phylo=std::string(argv[2])+"_run_"+std::to_string(i);
      std::ofstream out_phylogeny(outputname_phylo);
      out_phylogeny   <<currentphylos.printNewick() <<endl;
    }

  //Extract the best phylogeny for further display 'n' such.
  Phylogeny bestphylos=runs[min].phylos;



  cout<< "mutation prob: "          << runs[min].mutation_probability
      << ", temerature sd: "        << runs[min].temperature_drift_sd
      << ", similarity threshold: " << runs[min].species_sim_thresh
      << ", number of species: "    << runs[min].nspecies
      <<endl;

  std::ofstream out_persistgraph(argv[1]);
  std::ofstream out_phylogeny   (argv[2]);
  bestphylos.persistGraph(out_persistgraph);
  out_phylogeny   <<bestphylos.printNewick() <<endl;

  return 0;
}
