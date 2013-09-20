#include "mtbin.hpp"
#include "data.hpp"
#include "utility.hpp"
#include <algorithm>
#include <cassert>
#include <random>
#include <cmath>

std::default_random_engine rgen;

MtBin::MtBin(){
  startofdead=0;
}

MtBin::MtBin(double binx) : MtBin(binx) {}

void MtBin::killSalamander(int i) {
  assert(startofdead>0);
  bin[i].dead=true;
  std::swap(bin[i],bin[startofdead-1]);
  --startofdead;
}

void MtBin::mortaliate(double t) {
  ///If there are no living salamanders, then don't do anything
  if(startofdead==0) return;

  double mytemp=temp(t); //Current temperature of bin
  unsigned int maxalive=kkap(t);   //Current carrying capacity of the bin

  //For each salamander, check to see if it dies
  for(unsigned int s=0;s<startofdead && s<maxalive;++s){
    if(bin[s].pDie(mytemp))
      killSalamander(s);
  }

  //Since the carrying capacity of the bin may have been reduced since we last
  //mortaliated, kill at random individuals until we are within carrying capacity

  //Make a random number generator which points uniformly at all living individuals
  std::uniform_int_distribution<int> rdist(0, alive()-1);
  while(alive()>maxalive){
    int i=rdist(rgen);          //Get an individual
    if(bin[i].dead) continue;   //If the individual is already dead, ignore it
    killSalamander(i);
  }
}

double MtBin::temp(double t) const {
  assert(t>=0);
  assert(t<=65000);

  double h=area(t); //Get height of the mountain at this time
  //Find out the temperature adjustment for that height
  //assuming a dry air adiabatic lapse rate of 9.8 degC per vertical kilometer
  double altitude_temp_adjust=-9.8*h;

  //Interpolate temperature for this time
  int t0=(int)t;
  double ta=temps[t0];
  double tb=temps[t0+1];
  return ta+(tb-ta)*(t-t0)+altitude_temp_adjust;
}

double MtBin::area(double t) const {
  return MountainArea(height, t);
}

unsigned int MtBin::kkap(double t) const {
  //Maximum elevation of the mountain range over time
  double maxelevation = 2.8-1.846154e-05*t;
  double minarea = MountainArea(maxelevation, t);
  
  return std::min(area(t)/minarea, (double) binmax); //Returns a number [1, binmax], with 1 being the size of the smallest bin at time t
}

void MtBin::killAll() {
  for(auto &s: bin)
    s.dead=true;
  startofdead=0;
}

void MtBin::addSalamander(const Salamander &s) {
  if(startofdead==bin.size()) return;
  bin[startofdead]=s;
  ++startofdead;
}

unsigned int MtBin::alive() const {
  return bin.size()-startofdead;
}

void MtBin::breed(double t){
  unsigned int maxalive=kkap(t);   //Current carrying capacity of the bin

  if(alive()>=kkap(t)) return;     //The bin is too full for us to breed

  ///Make a random number generator that considers only the parents
  std::uniform_int_distribution<int> rdist(0, alive()-1);
  for(int i=0;i<10;++i){
    int parenta=rdist(rgen);
    int parentb=rdist(rgen);
    if(parenta.similar
  }
}
