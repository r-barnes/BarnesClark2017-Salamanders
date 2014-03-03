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

//Profile of the mountain
//We allocate this in the global namespace to prevent a stack overflow

Phylogeny RunSimulation(double mutation_probability, double species_sim_thresh){
  vector<MtBin> mts;
  mts.resize(1000);

  ////////////////////////////////////
  //INITIALIZE
  ////////////////////////////////////

  //Start off by killing everything
  for(auto &m: mts)
    m.killAll();

  Salamander Eve;
  //Eve is her own ancestor so that her children have the correct parent
  //relationship
  Eve.parent=0; 
  Eve.otemp =33.5618604122814;//degC
  //This is the global average temperature at sea level 65 million years ago.
  //Today, the optimal temperature for salamanders is 12.804279 degC
  Eve.genes =(Salamander::genetype)169113934842922489;
  Eve.dead  =false;
  Eve.mutation_probability=mutation_probability;

  //populate the first bin (m=0) with s populations of Eve.
  for(unsigned int m=0;m<1;++m)
  for(unsigned int s=0;s<10;++s){
    mts[m].addSalamander(Eve);
  }
  Phylogeny phylos(Eve, 0);

  ////////////////////////////////////
  //MAIN LOOP
  ////////////////////////////////////


  //Find cell corresponding to optimal happy temperature
  //Make a salamander in that cell alive
  //Go forth

  //Salamanders are at their very happiest at 12.804279 degC
  //At sea level 65MYA the temp was 33.5618604122814 degC
  //Temperature drops at 9.8 degC/1000m (adiabatic lapse rate of dry air)
  //So (33.5618604122814-12.804279)/9.8=2.1181*1000m=2.11km

  //Loop over years, starting at t=0, with t in units of millions of years
  for(double t=0;t<65.001;t+=0.5){
    //cerr<<"#"<<t<<endl;
    unsigned int population_size=0;                  //Sum up population size at time = t.
    for(unsigned int m=0;m<mts.size();++m){          //Increment across all mountain bins
      mts[m].mortaliate(t);                          //Kill individuals in the bin
      mts[m].breed(t, species_sim_thresh);           //Breed individuals in the bin
      if(m>0)            mts[m].diffuse(t,mts[m-1]); //Diffuse salamanders "up" the mountain
      if(m<mts.size()-1) mts[m].diffuse(t,mts[m+1]); //Diffuse salamanders "down" the mountain
      population_size+=mts[m].startofdead;           //Indicate beginning of vector to show where "living" salamanders end
    }
    phylos.UpdatePhylogeny(t, mts, species_sim_thresh);

    //cerr<<"#Species count="<<phylos.nodes.size()<<endl;
    //cerr<<"#Population size="<<population_size<<endl;
  }
  
  return phylos;
}

int main(){
  //srand (time(NULL)); //TODO: Uncomment this line before production

  //Phylogeny phylos=RunSimulation(1e-4, 0.98);
  //cout<<phylos.printACL2(65)<<endl;


  //Create a run object
  //Vary mutation rate and the speciation threshold
  //Record information on number of living taxa and their relationship
  struct run {
    //mutation rate
    double mutation_probability;
    //threshold for genetic differences before declaring something a new species
    double sim_thresh;
    //number of living species (count during simulations)
    int nalive;
    //empirical cumulative distribution of average branch lengths between extant taxa
    double ecdf;
    Phylogeny phylos;
  };

  //Create a vector to hold the runs
  std::vector<struct run> runs;

 //Set up some runs
  for(double mutation_probability=1e-4; mutation_probability<1e-3; mutation_probability+=1e-4)
  for(double sim_thresh=0.95; sim_thresh<1; sim_thresh+=0.01){
    struct run temp;
    temp.mutation_probability=mutation_probability;
    temp.sim_thresh=sim_thresh;
    runs.push_back(temp);
  }

  //Run the runs, stash the results
  #pragma omp parallel for
  for(unsigned int i=0;i<runs.size();++i){
    // Make sure that processors don't compete when writing position in simulation
    #pragma omp critical
      cout<<"Run #"<<i<<endl;
    Phylogeny phylos=RunSimulation(runs[i].mutation_probability, runs[i].sim_thresh);
    runs[i].nalive=phylos.numAlive(65); // Record output after 65my
    runs[i].ecdf=phylos.compareECDF(65);
    runs[i].phylos=phylos;
  }

  int min=0;
  for(unsigned int i=1;i<runs.size();++i)
    // Mark as best match if there are 80-120 species alive at the end of the
    //simulation, and if the ecdf is more similar to the Kozak and Wiens data
    //based on ecdf of mean branch distances than previous results.
    if(runs[i].ecdf<runs[min].ecdf && (80<runs[i].nalive && runs[i].nalive<120))
      min=i;

  Phylogeny bestphylos=runs[min].phylos;

  bestphylos.print("");
  cout<< "mutation prob: " << runs[min].mutation_probability
      << ", similarity threshold: " << runs[min].sim_thresh
      << ", number of species: " << runs[min].nalive
      <<endl;

}
