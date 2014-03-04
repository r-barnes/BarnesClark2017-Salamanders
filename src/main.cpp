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

  @param[in] mutation_probability  Mutation probability in TODO
  @param[in] species_sim_thresh    Threshold for two individuals to be considered
                                   members of the same species. Has value [0,1].

  @returns   A phylogeny object
*/
Phylogeny RunSimulation(double mutation_probability, double species_sim_thresh){
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

  //This is the global average temperature at sea level 65 million years ago.
  //Today, the optimal temperature for salamanders is 12.804279 degC.
  //We assume that Eve is well-adapted for her time by setting her optimal
  //temperature to be this global average sea level temperature.
  Eve.otemp = 33.5618604122814; //degC

  //We set Eve initially to have a genome in which all of the bits are off.
  //Since the genomes are used solely to determine speciation and speciation is
  //determined by the number of bits which are different between two genomes,
  //any starting value could be used with equal validity.
  Eve.genes = (Salamander::genetype)0;

  Eve.mutation_probability = mutation_probability;

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
  
  return phylos;
}

int main(){
  //srand (time(NULL)); //TODO: Uncomment this line before production

  Phylogeny phylos=RunSimulation(0.001, 0.96);
  phylos.print("");
  cout<<phylos.printNewick()<<endl;

  return 0;

  //Create a run struct. This struct will store the mutation and species
  //similarity thresholds used to run each of the simulations, along with the
  //outputs of those simulations. This allows us to easily multi-thread the
  //program.
  struct run {
    //A value [0,1] indicating the probability TODO
    double mutation_probability;
    //A value [0,1] indicating how similar to salamanders must be to be the same species
    double sim_thresh;
    //Number of living species at the end of the simulation
    int nalive;
    //Empirical cumulative distribution (ECDF) of average branch lengths between extant taxa
    double ecdf;
    //Phylogeny resulting from running the simulation
    Phylogeny phylos;
  };

  //Create a vector to hold the runs
  std::vector<struct run> runs;

  //Set up the runs
  for(double mutation_probability=1e-4; mutation_probability<1e-3; mutation_probability+=1e-4)
  for(double sim_thresh=0.95; sim_thresh<1; sim_thresh+=0.01){
    struct run temp;
    temp.mutation_probability = mutation_probability;
    temp.sim_thresh = sim_thresh;
    runs.push_back(temp);
  }

  //Run the runs and stash the results.
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    Phylogeny phylos = RunSimulation(runs[i].mutation_probability, runs[i].sim_thresh);
    runs[i].nalive   = phylos.numAlive(65);    //Record number alive at present day
    runs[i].ecdf     = phylos.compareECDF(65); //Record ECDF of those species alive at present day
    runs[i].phylos   = phylos;
  }

  int min=0;
  for(unsigned int i=1;i<runs.size();++i)
    //Mark as best match if there are 80-120 species alive at the end of the
    //simulation, and if the ECDF is more similar to the Kozak and Wiens data
    //based on ECDF of mean branch distances than previous results.
    if(runs[i].ecdf<runs[min].ecdf && (80<runs[i].nalive && runs[i].nalive<120))
      min=i;

  //Extract the best phylogeny for further display 'n' such.
  Phylogeny bestphylos=runs[min].phylos;

  bestphylos.print("");
  cout<< "mutation prob: "          << runs[min].mutation_probability
      << ", similarity threshold: " << runs[min].sim_thresh
      << ", number of species: "    << runs[min].nalive
      <<endl;

  cout<<bestphylos.printNewick()<<endl;
}
