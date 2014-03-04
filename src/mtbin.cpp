#include "mtbin.hpp"
#include "data.hpp"
#include "utility.hpp"
#include <algorithm>
#include <cassert>
#include <random>
#include <cmath>


///Give elevation as altitude in kilometers, and t in kiloyears ago (i.e. 0.001MYA)
double MountainArea(double elevation, double timeMyrs) {
  ///Constants defining a normal distribution that describes area available at
  ///different height bands in the Appalachian mountains. Paramteres are fit to
  ///contemporary height distributions presented in Kozak and Weins 2010.
  ///deltasd describes change in sd per year, which  
  const double elek     = 12.029753;
  const double elesigma = 0.211410;
  const double elemu    = 0.245547;
  const double pi       = 3.141593;
  const double deltasd  = 2.881572e-09;
  //Input "t" is in millions of years - transform this into thousands of years
  double timeKyrs = timeMyrs*1000;

  double area = elek
                * ( 1/sqrt( 2*pi*pow(elesigma+timeKyrs*deltasd*1000, 2) ))
                * exp(
                        ( -pow(elevation-elemu, 2) )
                      / (2 * pow(elesigma+timeKyrs*deltasd*1000, 2))
                  );
  area *= 100000; //Convert area to units of km2
  return area;
}



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

  double       mytemp  =temp(t);  //Current temperature of bin
  unsigned int maxalive=kkap(t);  //Current carrying capacity of the bin

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

  //Find out the temperature adjustment for that height
  //assuming a dry air adiabatic lapse rate of 9.8 degC per vertical kilometer
  double altitude_temp_adjust=-9.8*height;

  //Interpolate temperature for this time
  int t0    = (int)t;
  double ta = temps[t0];
  double tb = temps[t0+1];
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
  return startofdead;
}

void MtBin::breed(double t, double species_sim_thresh){
  if(startofdead==0) return;       //No one is alive here. There can be no breeding.

  unsigned int maxalive=kkap(t);   //Current carrying capacity of the bin
  if(alive()>=maxalive) return;    //The bin is too full for us to breed

  //Maximum number of tries to find a pair to mate
  int maxtries=4*(maxalive-alive());

  int max_babies=10;

  ///Make a random number generator that considers only the parents
  std::uniform_int_distribution<int> rdist(0, alive()-1);
  while(alive()<maxalive && max_babies>=0 && maxtries>=0){ //TODO: Should be maxtries--
    Salamander &parenta=bin[rdist(rgen)];
    Salamander &parentb=bin[rdist(rgen)];
    if(parenta.pSimilar(parentb, species_sim_thresh)){
      addSalamander(parenta.breed(parentb));
      max_babies--;
    }
  }
}

void MtBin::diffuse(double t, MtBin &b) {
  std::uniform_int_distribution<int> myguys   (0,   alive()-1);
  std::uniform_int_distribution<int> otherguys(0, b.alive()-1);
  unsigned int aswapn=  alive()/10;
  unsigned int bswapn=b.alive()/10;

  int swapc=std::min(aswapn,bswapn);
  for(int i=0;i<swapc;++i)
    std::swap(bin[myguys(rgen)], b.bin[otherguys(rgen)]);

  unsigned int ka=  kkap(t);
  unsigned int kb=b.kkap(t);

  if(aswapn>bswapn){
    aswapn-=swapc;
    aswapn=std::min(aswapn,alive()-ka);
    for(unsigned int i=0;i<aswapn;++i){
      int temp=myguys(rgen);
      b.addSalamander(bin[temp]);
      killSalamander(temp);
    }
  } else if(aswapn<bswapn) {
    bswapn-=swapc;
    bswapn=std::min(bswapn,b.alive()-kb);
    for(unsigned int i=0;i<bswapn;++i){
      int temp=otherguys(rgen);
      addSalamander(b.bin[temp]);
      b.killSalamander(temp);
    }
  }
}
