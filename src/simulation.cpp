#include "simulation.hpp"
#include "params.hpp"
#include <iostream>
#include <limits>

void Simulation::runSimulation(){
  //65Mya the Appalachian Mountains were 2.8km tall. Initialize each bin to
  //point to its given elevation band.
  mts.reserve(numbins);
  for(int m=0;m<numbins;m++)
    mts.push_back(MtBin(m*2.8/numbins, TheParams::get().pVaryHeight()));

  ////////////////////////////////////
  //INITIALIZE
  ////////////////////////////////////

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

  if(TheParams::get().initialAltitude()<0 || (int)mts.size()<=TheParams::get().initialAltitude()){
    std::cerr<<"Initial bin was outside of range. Should be in [0,"<<(mts.size()-1)<<"]."<<std::endl;
    throw std::runtime_error("Initial bin outside of range.");
  }

  //We populate the first (lowest) mountain bin with some Eve-clones. We
  //populate only the lowest mountain bin because that mountain bin will have a
  //temperature close to the global average which is optimal for Eve (see above).
  for(unsigned int s=0;s<10;++s)
    mts[TheParams::get().initialAltitude()].addSalamander(Eve);

  //Begin a new phylogeny with Eve as the root
  phylos=Phylogeny(Eve, 0);

  ////////////////////////////////////
  //MAIN LOOP
  ////////////////////////////////////

  //Loop over years, starting at t=0, which corresponds to 65 million years ago.
  //tMyrs is in units of millions of years
  double tMyrs=0;
  for(tMyrs=0;tMyrs<65.000;tMyrs+=TheParams::get().timestep()){
    //This requires a linear walk of all the bins on the mountain. Hence, it's a
    //little expensive. But it prevents many walks below if all the salamanders
    //go extinct early on. Therefore, in a parameter space where many
    //populations won't make it, this is a worthwhile thing to do.
    if(alive()==0) break;

    //Visit death upon each bin
    for(auto &m: mts)
      m.mortaliate(tMyrs);

    //Let the salamanders in each bin be fruitful, and multiply
    for(auto &m: mts)
      m.breed(tMyrs);

    //For each bin, offer some salamanders therein the opportunity to migrate up
    //or down the mountain. NOTE: Another way of doing this would be to visit
    //bins in a random order. Since carrying capacities determine migration
    //success, it is easy to move up the mountain than down. However, this would
    //be true anyway becaues the bottom of the mountain typically has larger
    //populations.
    if(TheParams::get().dispersalType()==DISPERSAL_BETTER){
      for(unsigned int m=0;m<mts.size();++m){
        if(m==0)                 mts[m].diffuseToBetter(tMyrs, TheParams::get().dispersalProb(), nullptr,   &mts[m+1]);
        else if(m==mts.size()-1) mts[m].diffuseToBetter(tMyrs, TheParams::get().dispersalProb(), &mts[m-1], nullptr  );
        else                     mts[m].diffuseToBetter(tMyrs, TheParams::get().dispersalProb(), &mts[m-1], &mts[m+1]);
      }
    } else if(TheParams::get().dispersalType()==DISPERSAL_MAYBE_WORSE) {
      for(unsigned int m=0;m<mts.size();++m){
        if(m==0)                 mts[m].diffuseLocal(tMyrs, TheParams::get().dispersalProb(), nullptr,   &mts[m+1]);
        else if(m==mts.size()-1) mts[m].diffuseLocal(tMyrs, TheParams::get().dispersalProb(), &mts[m-1], nullptr  );
        else                     mts[m].diffuseLocal(tMyrs, TheParams::get().dispersalProb(), &mts[m-1], &mts[m+1]);
      }
    } else if(TheParams::get().dispersalType()==DISPERSAL_GLOBAL) {
      for(auto &m: mts)
        m.diffuseGlobal(tMyrs, TheParams::get().dispersalProb(), mts);
    }

    //Updates the phylogeny based on the current time, living salamanders, and
    //species similarity threshold
    phylos.UpdatePhylogeny(tMyrs, TheParams::get().timestep(), mts, TheParams::get().speciesSimthresh());
  }

  //Records the time at which the simulation ended
  if(tMyrs>=65) tMyrs-=TheParams::get().timestep(); //Since the last step goes past the end of time
  endtime = tMyrs;

  //Records the average optimal temperature of the salamanders alive at present
  //day
  avg_otempdegC = AvgOtempdegC();

  //Records number of salamanders alive at present day
  salive        = alive();

  //Record number of species alive at present day
  nspecies      = phylos.livingSpecies(endtime);

  //Record mean branch distance ECDF of those species alive at present day
  ecdf          = phylos.compareECDF(endtime);

  //Record the average elevation at which salamanders are found at the end of
  //the simulation
  avg_elevation = AvgElevation();

  //Destroy all of the salamanders and mountain bins so that the simulation is
  //not using excessive memory when it is not being run. We don't need this
  //information anyway because we capture it in the summary statistics above.
  mts.clear();
  mts.shrink_to_fit();
}


//This calculates the total number of living salamanders
int Simulation::alive() const {
  int sum = 0;
  for(const auto &m: mts)
    sum += m.alive();
  return sum;
}


//This finds the average optimal temperature of all of the salamanders
double Simulation::AvgOtempdegC() const {
  int alive  = 0;
  double avg = 0;
  for(const auto &m: mts)
  for(const auto &s: m.bin){
    avg += s.otempdegC;
    alive++;
  }

  return avg/alive;
}


//Replaces the stored phylogeny (which may take up quite a bit of memory!) with
//a new, empty phylogeny object. This causes the old phylogeny to get recycled,
//thereby freeing the memory.
void Simulation::dumpPhylogeny() {
  phylos = Phylogeny();
}


//Returns the average elevation of all the salamanders on the mountain
double Simulation::AvgElevation() const {
  double salamander_count   = 0;
  double weighted_elevation = 0;
  for(const auto &m: mts){
    salamander_count   += m.alive();
    weighted_elevation += m.alive()*m.height();
  }
  return weighted_elevation/salamander_count;
}