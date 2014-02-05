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
  Eve.otemp =33.5618604122814;//degC //This is the global average temperature at sea level 65 million years ago.
                                     //Today, the optimal temperature for salamanders is 12.804279 degC
  Eve.genes =(Salamander::genetype)169113934842922489;
  Eve.dead  =false;
  Eve.mutation_probability=mutation_probability;

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

  for(double t=0;t<65.001;t+=0.5){
          cerr<<"#"<<t<<endl;
    unsigned int population_size=0;
    for(unsigned int m=0;m<mts.size();++m){
      mts[m].mortaliate(t);
      mts[m].breed(t, species_sim_thresh);
      if(m>0)            mts[m].diffuse(t,mts[m-1]);
      if(m<mts.size()-1) mts[m].diffuse(t,mts[m+1]);
      population_size+=mts[m].startofdead;
    }
    phylos.UpdatePhylogeny(t, mts, species_sim_thresh);

    cerr<<"#Species count="<<phylos.nodes.size()<<endl;
    cerr<<"#Population size="<<population_size<<endl;
  }
  
  return phylos;
}

int main(){
  //srand (time(NULL)); //TODO: Uncomment this line before production

  //RunSimulation(1e-4, 0.98);

  struct run {
    double mutation_probability;
    double sim_thresh;
    int nalive;
    double ecdf;
    Phylogeny phylo;
  };

  std::vector<struct run> runs;

  for(double mutation_probability=1e-4; mutation_probability<1e-3; mutation_probability+=1e-5)
<<<<<<< HEAD
  for(double sim_thresh=0.95; sim_thresh<1; sim_thresh+=0.01) {
    Phylogeny phylos=RunSimulation(mutation_probability, sim_thresh);
         
    int n_alive=phylos.numAlive(65);
    if( !(50<=n_alive && n_alive<=150) ) continue;
    
    double ecdf=phylos.compareECDF(65);
    if(ecdf<run_result.ecdf){    
      run_result.mutation_probability=mutation_probability;
      run_result.sim_tresh=sim_thresh;
      run_result.nalive=n_alive;
      run_result.ecdf=phylos.compareECDF(65);
      run_result.phylo=phylos;
    }
=======
  for(double sim_thresh=0.5; sim_thresh<1; sim_thresh+=0.01){
    struct run temp;
    temp.mutation_probability=mutation_probability;
    temp.sim_thresh=sim_thresh;
    runs.push_back(temp);
  }

  for(unsigned int i=0;i<runs.size();++i){
    Phylogeny phylos=RunSimulation(runs[i].mutation_probability, runs[i].sim_thresh);
    runs[i].nalive=phylos.numAlive(65);
    runs[i].ecdf=phylos.compareECDF(65);
    runs[i].phylo=phylos;
>>>>>>> 86d8c681c5e191e42b9e227fc0df1d58850a0a4b
  }

//run_result.phylo.print("");
//cout<< "mutation prob: " << run_result.mutation_probability << ", similarity threshold: " << run_result.sim_tresh << ", number of species: " << run_result.nalive<<endl;
}
