#include "simulation.hpp"
#include "params.hpp"
#include "random.hpp"
#include "temp.hpp"
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>

void Simulation::runSimulation(){
  //65Mya the Appalachian Mountains were 2.8km tall. Initialize each bin to
  //point to its given elevation band.
  mts.reserve(TheParams.numBins());
  for(int m=0;m<TheParams.numBins();m++)
    mts.push_back(MtBin(m*2.8/TheParams.numBins()));

  //This vector is shuffled before each dispersion event to ensure that there is
  //no bias towards upwards or downards movement on the mountain
  std::vector<unsigned int> mtbin_order;
  for(unsigned int i=0;i<mts.size();i++)
    mtbin_order.push_back(i);

  if(TheParams.initialAltitude()<0 || (int)mts.size()<=TheParams.initialAltitude()){
    std::cerr<<"Initial bin was outside of range. Should be in [0,"<<(mts.size()-1)<<"]."<<std::endl;
    throw std::runtime_error("Initial bin outside of range.");
  }

  ////////////////////////////////////
  //INITIALIZE
  ////////////////////////////////////

  //Eve is the first salamander species from which all others will emerge
  Salamander Eve;

  //Eve is her own ancestor so that her children have the correct parent
  //relationship
  Eve.parent = 0;

  //The simulation assumes that it is given a time series of mean summer diurnal
  //temperatures corresponding to the  temperatures at the base of the
  //mountains. We ensure that Eve is well-adapted for her time by setting her
  //optimal temperature to be equal to the temperature of the bin she starts in.
  Eve.otempdegC = Temperature.getTemp(0)-9.8*mts[TheParams.initialAltitude()].heightkm();

  //We set Eve initially to have a genome in which all of the bits are off.
  //Since the genomes are used solely to determine speciation and speciation is
  //determined by the number of bits which are different between two genomes,
  //any starting value could be used with equal validity.
  Eve.genes = (Salamander::genetype)0;

  //We populate the first (lowest) mountain bin with some Eve-clones. We
  //populate only the lowest mountain bin because that mountain bin will have a
  //temperature close to the global average which is optimal for Eve (see above).
  for(int s=0;s<TheParams.initialPopSize();++s)
    mts[TheParams.initialAltitude()].addSalamander(Eve);

  //Begin a new phylogeny with Eve as the root
  phylos = Phylogeny(Eve, 0);

  ////////////////////////////////////
  //MAIN LOOP
  ////////////////////////////////////

  //Loop over years, starting at t=0, which corresponds to 65 million years ago.
  //tMyrs is in units of millions of years
  double tMyrs=0;
  for(tMyrs=0;tMyrs<65.001;tMyrs+=TheParams.timestep()){
    //This requires a linear walk of all the bins on the mountain. Hence, it's a
    //little expensive. But it prevents many walks below if all the salamanders
    //go extinct early on. Therefore, in a parameter space where many
    //populations won't make it, this is a worthwhile thing to do.
    if(alive()+surrounding_lowlands.alive()==0) break;

    if(TheParams.debug())
      printMt(tMyrs);

    //Visit death upon each bin
    for(auto &m: mts)
      m.mortaliate(tMyrs);

    //Ensure that there are no Sky Salamanders in the simulation. Mountains
    //erode over time, the bins which are above the mountains' actual heights
    //must be emptied of their inhabitants.
    for(auto &m: mts)
      if(m.heightkm()>=MtBin::heightMaxKm(tMyrs))
        m.killAll();

    //Let the salamanders in each bin be fruitful, and multiply
    for(auto &m: mts)
      m.breed(tMyrs);

    surrounding_lowlands.breed(tMyrs);

    //Randomize the order in which we visit bins so there is no upwards or
    //downwards bias to movement. Such a bias could arise, say, by always
    //considering bins from bottom to top. In this case, a salamander at the
    //bottom would have an opportunity to move up several bins whereas no
    //salamander would be able to move downwards more than one bin.
    //std::random_shuffle() is not thread safe, so we use the Fisher-Yates-
    //Durstenfeld-Knuth algorithm.
    assert(mtbin_order.size()>2);
    for(unsigned int i=0;i<mtbin_order.size()-2;i++){
      unsigned int j = uniform_rand_int(i,mtbin_order.size()-1);
      std::swap(mtbin_order[i],mtbin_order[j]);
    }

    //For each bin, offer some salamanders therein the opportunity to migrate up
    //or down the mountain.
    if(TheParams.dispersalType()==DISPERSAL_BETTER){
      for(unsigned int mo=0;mo<mtbin_order.size();++mo){
        unsigned int m = mtbin_order[mo];
        if(m==0)                 mts[m].diffuseToBetter(tMyrs, nullptr,   &mts[m+1]);
        else if(m==mts.size()-1) mts[m].diffuseToBetter(tMyrs, &mts[m-1], nullptr  );
        else                     mts[m].diffuseToBetter(tMyrs, &mts[m-1], &mts[m+1]);
      }
    } else if(TheParams.dispersalType()==DISPERSAL_MAYBE_WORSE) {
      for(unsigned int mo=0;mo<mtbin_order.size();++mo){
        unsigned int m = mtbin_order[mo];
        if(m==0)                 mts[m].diffuseLocal(tMyrs, nullptr,   &mts[m+1]);
        else if(m==mts.size()-1) mts[m].diffuseLocal(tMyrs, &mts[m-1], nullptr  );
        else                     mts[m].diffuseLocal(tMyrs, &mts[m-1], &mts[m+1]);
      }
    } else if(TheParams.dispersalType()==DISPERSAL_GLOBAL) {
      //We don't need to randomize the order for global dispersion since it
      //contains no bias.
      for(auto &m: mts)
        m.diffuseGlobal(tMyrs, mts);
    }

    //Randomize order of execution to smooth biases
    if(uniform_rand_real(0,1)>=0.5){
      mts[0].diffuseToLowlands(surrounding_lowlands);
      surrounding_lowlands.diffuseToLowlands(mts[0]);
    } else {
      surrounding_lowlands.diffuseToLowlands(mts[0]);
      mts[0].diffuseToLowlands(surrounding_lowlands);
    }

    //Updates the phylogeny based on the current time, living salamanders, and
    //species similarity threshold
    phylos.UpdatePhylogeny(tMyrs, TheParams.timestep(), mts);
  }

  //Records the time at which the simulation ended
  if(tMyrs>=65.001) tMyrs-=TheParams.timestep(); //Since the last step goes past the end of time
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
    weighted_elevation += m.alive()*m.heightkm();
  }
  return weighted_elevation/salamander_count;
}

void Simulation::printMt(double tMyrs) const {
  unsigned int nalive   = alive();
  unsigned int maxalive = 0;
  for(const auto &m: mts)
    maxalive = std::max(maxalive,m.alive());
  std::cout<<"##########################\n";
  std::cout<<"t:           "<<tMyrs<<"\n";
  std::cout<<"Total Pop:   "<<alive()+surrounding_lowlands.alive()<<"\n";
  std::cout<<"Lowland Pop: "<<surrounding_lowlands.alive()<<"\n";
  std::cout<<std::setw(2)<<"i"<<" "<<std::setw(5)<<"elev"<<"  "<<std::setw(20)<<"Dist"<<"  "<<"       %  Count\n";
  std::cout<<std::setw(2)<<"-"<<" "<<std::setw(5)<<"LOW" <<" |"<<std::setw(20)<<" "<<"| "<<"         "<<std::setw(6)<<surrounding_lowlands.alive()<<"\n";
  for(unsigned int i=0;i<mts.size();i++){
    std::cout<<std::setw(2)<<i<<" "<<std::setw(5);
    if(mts[i].heightkm()<MtBin::heightMaxKm(tMyrs))
      std::cout<<mts[i].heightkm()<<" |"<<std::setw(20)<<std::string(20*mts[i].alive()/maxalive,'#')<<"| "<<std::setw(7)<<std::setprecision(4)<<(mts[i].alive()/(double)nalive*100.0)<<"% "<<std::setw(6)<<mts[i].alive();
    else
      std::cout<<"XXXXX"<<" |"<<std::setw(20)<<std::string(20,'-')<<"|";
    std::cout<<"\n";
  }
  std::cout<<"\n\n";
}