#include "mtbin.hpp"
#include "data.hpp"
#include <algorithm>
#include <cassert>

MtBin::MtBin(){
  startofdead=0;
}

MtBin::MtBin(double binx) : MtBin(binx) {}

void MtBin::mortaliate(double t) {
  ///If there are no living salamanders, then don't do anything
  if(startofdead==0) return;

  double mytemp=temp(t);

  for(unsigned int s=0;s<startofdead;++s){
    if(bin[s].pDie(mytemp)){
      bin[s].dead=true;
      std::swap(bin[s],bin[startofdead-1]);
      --startofdead;
    }
  }
}

double MtBin::temp(double t) const {
  assert(t>=0);
  assert(t<=65000);

  double h=height(t); //Get height of the mountain at this time
  //Find out the temperature adjustment for that height
  //assuming a dry air adiabatic lapse rate of 9.8 degC per vertical kilometer
  double altitude_temp_adjust=-9.8*h;

  //Interpolate temperature for this time
  int t0=(int)t;
  double ta=temps[t0]  +altitude_temp_adjust;
  double tb=temps[t0+1]+altitude_temp_adjust;
  return ta+(tb-ta)*(t-t0);
}

double MtBin::area(double t) const {
  return MountainArea(height, t);
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
