#include "simulation.hpp"
#include <iostream>

Simulation::Simulation(double mutation_probability0, double temperature_drift_sd0, double species_sim_thresh0){
  mutation_probability = mutation_probability0;
  temperature_drift_sd = temperature_drift_sd0;
  species_sim_thresh   = species_sim_thresh0;
  avg_otempdegC        = 0;
  salive               = 0;
  nspecies             = 0;
}

void Simulation::runSimulation(){
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
  phylos=Phylogeny(Eve, 0);

  ////////////////////////////////////
  //MAIN LOOP
  ////////////////////////////////////

  //Loop over years, starting at t=0, which corresponds to 65 million years ago.
  //tMyrs is in units of millions of years
  for(double tMyrs=0;tMyrs<65.001;tMyrs+=0.5){
    std::cerr<<"Alive at begin "<<tMyrs<<": "<<alive()<<std::endl;

    //Increment up the mountain
    for(auto &m: mts)
      m.mortaliate(tMyrs);                               //Kill individuals in the bin

    std::cerr<<"Alive after mortaliate "<<tMyrs<<": "<<alive()<<std::endl;

    for(auto &m: mts)
      m.breed(tMyrs, species_sim_thresh);                //Breed individuals in the bin

    for(unsigned int m=0;m<mts.size();++m){
      if(m>0)            mts[m].diffuse(tMyrs,mts[m-1]); //Diffuse salamanders "down" the mountain
      if(m<mts.size()-1) mts[m].diffuse(tMyrs,mts[m+1]); //Diffuse salamanders "up" the mountain
    }

    //Updates the phylogeny based on the current time, living salamanders, and
    //species similarity threshold
    phylos.UpdatePhylogeny(tMyrs, mts, species_sim_thresh);
  }

  avg_otempdegC = AvgOtempdegC();
  salive        = alive();
  nspecies      = phylos.numAlive(65);    //Record number of species alive at present day
  ecdf          = phylos.compareECDF(65); //Record ECDF of those species alive at present day

  mts=std::vector<MtBin>();
}

int Simulation::alive() const {
  int sum=0;
  for(const auto &m: mts)
    sum+=m.alive();
  return sum;
}

double Simulation::AvgOtempdegC() const {
  int alive=0;
  double avg=0;
  for(const auto &m: mts)
  for(const auto &s: m.bin){
    avg+=s.otempdegC;
    alive++;
  }

  return avg/alive;
}