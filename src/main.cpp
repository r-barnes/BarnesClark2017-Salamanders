#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <time.h>
using namespace std;


/**
  This function runs a simulation with the specified mutation probability and
  species similarity threshold. It returns the resulting phylogeny. The function
  is thread-safe.

  @param[in]  mutation_probability  Mutation probability in TODO
  @param[in]  temperature_drift_sd  Std.Dev. of change in temperature tolerance
                                    that occurs between parents and children.
  @param[in]  species_sim_thresh    Threshold for two individuals to be considered
                                    members of the same species. Has value [0,1].
  @param[out] avg_otempdegC         Average optimal temperature of all surviving
                                    species at the end of the run


  @returns   A phylogeny object
*/
Phylogeny RunSimulation(double mutation_probability, double temperature_drift_sd, double species_sim_thresh, double &avg_otempdegC){
  vector<MtBin> mts;

  //65Mya the Appalachian Mountains were 2.8km tall. We decide, arbitrarily to
  //use 1000 bins to represent the mountain.
  mts.reserve(1000);
  for(int m=0;m<1000;m++)
    mts.push_back(MtBin(m*2.8/1000.0));


  ////////////////////////////////////
  //INITIALIZE
  ////////////////////////////////////

  //Kill everything since new salamanders are, by default, constructed alive.
  for(auto &m: mts)
    m.killAll();

  //Eve is the first salamander species from which all others will emerge
  Salamander Eve;

  //Eve is her own ancestor so that her children have the correct parent
  //relationship
  Eve.parent = 0; 

  //This is the mean summer diurnal temperature at sea level 65 million years
  //ago in Greensboro, NC. Today, the optimal temperature for salamanders is
  //12.804279 degC. We assume that Eve is well-adapted for her time by setting
  //her optimal temperature to be this global average sea level temperature.
  Eve.otempdegC = 33.5618604122814; //degC

  //We set Eve initially to have a genome in which all of the bits are off.
  //Since the genomes are used solely to determine speciation and speciation is
  //determined by the number of bits which are different between two genomes,
  //any starting value could be used with equal validity.
  Eve.genes = (Salamander::genetype)0;

  Eve.mutation_probability = mutation_probability;
  Eve.temperature_drift_sd = temperature_drift_sd;

  //We populate the first (lowest) mountain bin with some Eve-clones. We
  //populate only the lowest mountain bin because that mountain bin will have a
  //temperature close to the global average which is optimal for Eve (see above).
  for(unsigned int s=0;s<10;++s)
    mts[0].addSalamander(Eve);

  //Begin a new phylogeny with Eve as the root
  Phylogeny phylos(Eve, 0);

  ////////////////////////////////////
  //MAIN LOOP
  ////////////////////////////////////

  //Loop over years, starting at t=0, which corresponds to 65 million years ago.
  //tMyrs is in units of millions of years
  for(double tMyrs=0;tMyrs<65.001;tMyrs+=0.5){
    //Increment up the mountain
    for(unsigned int m=0;m<mts.size();++m){
      mts[m].mortaliate(tMyrs);                          //Kill individuals in the bin
      mts[m].breed(tMyrs, species_sim_thresh);           //Breed individuals in the bin
      if(m>0)            mts[m].diffuse(tMyrs,mts[m-1]); //Diffuse salamanders "down" the mountain
      if(m<mts.size()-1) mts[m].diffuse(tMyrs,mts[m+1]); //Diffuse salamanders "up" the mountain
    }

    //Updates the phylogeny based on the current time, living salamanders, and
    //species similarity threshold
    phylos.UpdatePhylogeny(tMyrs, mts, species_sim_thresh);
  }

  int alive_count=0;
  for(const auto &m: mts)
  for(const auto &s: m.bin){
    avg_otempdegC+=s.otempdegC;
    alive_count++;
  }

  avg_otempdegC/=alive_count;
  
  return phylos;
}

//TODO What is the "argc" argumnet?
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
    double avg_otempdegC;
    Phylogeny phylos=RunSimulation(0.001, 0.01, 0.96,avg_otempdegC);
    phylos.persistGraph(out_persistgraph);
    out_phylogeny   <<phylos.printNewick() <<endl;
    return 0;
  }

  //Create a run struct. This struct will store the mutation and species
  //similarity thresholds used to run each of the simulations, along with the
  //outputs of those simulations. This allows us to easily multi-thread the
  //program.
  struct run {
    //A value [0,1] indicating the probability TODO
    double mutation_probability;
    //A value [0, Inf] dictating the standard deviation of the random, normal
    //change in temperature tolerance that occurs between parents and offspring.
    double temperature_drift_sd;
    //A value [0,1] indicating how similar to salamanders must be to be the same species
    double sim_thresh;
    //Number of living species at the end of the simulation
    int nalive;
    //Empirical cumulative distribution (ECDF) of average branch lengths between extant taxa
    double ecdf;
    //Average optimal temperature of all the salamanders alive at the end of the run
    double avg_otempdegC;
    //Phylogeny resulting from running the simulation
    Phylogeny phylos;
  };

  //Create a vector to hold the runs
  std::vector<struct run> runs;

  //Set up the runs
  for(double mutation_probability=1e-4; mutation_probability<1e-3; mutation_probability+=1e-4)
  for(double temperature_drift_sd=1e-3; temperature_drift_sd<1e-1; temperature_drift_sd+=1e-2)
  for(double sim_thresh=0.95; sim_thresh<1; sim_thresh+=0.01){
    struct run temp;
    temp.mutation_probability = mutation_probability;
    temp.temperature_drift_sd = temperature_drift_sd;
    temp.sim_thresh = sim_thresh;
    runs.push_back(temp);
  }

  //Run the runs and stash the results.
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    Phylogeny phylos = RunSimulation(runs[i].mutation_probability, runs[i].temperature_drift_sd, runs[i].sim_thresh, runs[i].avg_otempdegC);
    runs[i].nalive   = phylos.numAlive(65);    //Record number alive at present day
    runs[i].ecdf     = phylos.compareECDF(65); //Record ECDF of those species alive at present day
    runs[i].phylos   = phylos;
  }


  int min=0;
  for(unsigned int i=1;i<runs.size();++i)
    //Mark as best match if there are 80-120 species alive at the end of the
    //simulation, and if the ECDF is more similar to the Kozak and Wiens data
    //based on ECDF of mean branch distances than previous results.
    //otemp limits are from Gifford (&Martin) 1970: Mean +/- 2SD for salamander
    //feeding optimum
    if(runs[i].ecdf<runs[min].ecdf && (80<runs[i].nalive && runs[i].nalive<120) &&
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
      << ", similarity threshold: " << runs[min].sim_thresh
      << ", number of species: "    << runs[min].nalive
      <<endl;

  std::ofstream out_persistgraph(argv[1]);
  std::ofstream out_phylogeny   (argv[2]);
  bestphylos.persistGraph(out_persistgraph);
  out_phylogeny   <<bestphylos.printNewick() <<endl;

  return 0;
}
