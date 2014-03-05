#include "mtbin.hpp"
#include "data.hpp"
#include "random.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>


///Generic Gaussian distribution function
double Gaussian(double x, double mean, double sigma){
  const double PI=3.14159265359;
  return 1/sigma/std::sqrt(2*PI)*std::exp(-std::pow(x-mean,2)/2/std::pow(sigma,2));
}


MtBin::MtBin(double heightkm0){
  heightkm    = heightkm0;
  bin.reserve(1000); //Maximum salamander populations per mountain bin
}


double MtBin::height() const {
  return heightkm;
}


void MtBin::killSalamander(int s) {
  assert(bin.size()>0);
  std::swap(bin.at(s),bin.back());
  bin.pop_back();
}


void MtBin::killSalamander(MtBin::container::iterator s) {
  assert(bin.size()>0);
  std::swap(*s,bin.back());
  bin.pop_back();
}


void MtBin::mortaliate(double tMyrs) {
  ///If there are no living salamanders, then don't do anything
  if(bin.empty()) return;

  double       mytemp   = temp(tMyrs);  //Current temperature of bin
  unsigned int maxalive = kkap(tMyrs);  //Current carrying capacity of the bin

  //For each salamander, check to see if it dies
  for(unsigned int s=0;s<bin.size() && s<maxalive;++s)
    if(bin.at(s).pDie(mytemp)) {
      killSalamander(s);
      //If we kill a salamander, we swap the last living salamander in the list
      //with the salamander we just killed. Therefore, we need to make sure
      //that we still run the mortaliate function for the living salamander that
      //now inhabits the spot that we just filled.
      s--;
    }

  //Since the carrying capacity of the bin may have been reduced since we last
  //mortaliated, kill individuals at random until we are within carrying capacity

  //Make a random number generator which points uniformly at all living individuals
  while(alive()>maxalive)
    killSalamander(randomSalamander());
}


double MtBin::temp(double tMyrs) const {
  //Input "t" is in millions of years - transform this into thousands of years
  double timeKyrs = tMyrs*1000;

  assert(timeKyrs >= 0);
  assert(timeKyrs <= 65000);

  //Find out the temperature adjustment for that height assuming a dry air
  //adiabatic lapse rate of 9.8 degC per vertical kilometer
  double altitude_temp_adjust = -9.8*heightkm;

  //Interpolate temperature for this time. Temperature is binned into 1kyr bins.
  int    t0    = (int)timeKyrs;  //Time of the start of the 1kyr bin
  double ta    = temps[t0];      //Temperature of this 1kyr bin
  double tb    = temps[t0+1];    //Temperature of the next 1kyr bin
  double tdiff = tb-ta;          //Temperature difference between the two bins
  double sea_level_temp = ta + tdiff*(timeKyrs-t0);
  return sea_level_temp + altitude_temp_adjust;
}


unsigned int MtBin::kkap(double tMyrs) const {
  //Input "t" is in millions of years - transform this into thousands of years
  double timeKyrs = tMyrs*1000;

  //TODO: The area isn't integrated over the height of the band. Is this bad?

  //Maximum elevation of the mountain range over time
  //Based on linear shrinking of mountain hight from 2.8km at 65Mya (according
  //to the USGS website on "Geologic Provinces of the Untied States: Appalachian
  //Highlands Province") to current elevation (1.6km) from Kozak and Wiens 2010.
  double const height_65mya = 2.8; //km
  double const height_0mya  = 1.6; //km
  double const erosion_rate = (height_65mya-height_0mya)/65000; //Erosion rate per 1kyr
  double maxelevation       = 2.8-erosion_rate*timeKyrs; //km
  double minarea_today      = area(maxelevation, tMyrs);
  
  //Returns a number [1, binmax], with 1 being the size of the smallest bin
  //at time tMyrs. This ensures that the smallest area (at the top of the
  //mountain) will always have a carrying capacity of at least 1 salamander.
  //TODO: Think about this more later.
  return std::min( area(heightkm, tMyrs)/minarea_today, (double) binmax);
}


void MtBin::killAll() {
  bin.clear();
}


void MtBin::addSalamander(const Salamander &s) {
  bin.push_back(s);
}


unsigned int MtBin::alive() const {
  return bin.size();
}


void MtBin::breed(double t, double species_sim_thresh){
  if(bin.empty()) return;          //No one is alive here. There can be no breeding.

  unsigned int maxalive=kkap(t);   //Current carrying capacity of the bin
  if(alive()>=maxalive) return;    //The bin is too full for us to breed

  //Maximum number of tries to find a pair to mate - prevents infinite loop
  int maxtries=40*(maxalive-alive());

  //Maximum number of new offspring per bin per unit time
  int max_babies=10;

  ///Make a random number generator that considers only the parents
  std::uniform_int_distribution<int> rdist(0, alive()-1);

  //As long as there's room in the bin, and we still have to make babies, and we
  //are not caught in an infinite loop, then try to make more babies.
  while(alive()<maxalive && max_babies>=0 && maxtries-->0){
    Salamander &parenta=bin.at(rdist(rand_engine()));
    Salamander &parentb=bin.at(rdist(rand_engine()));
    //If parents are genetically similar enough to be classed as the same species
    //based on species_sim_thresh, then they can breed.
    if(parenta.pSimilar(parentb, species_sim_thresh)){
      addSalamander(parenta.breed(parentb));
      max_babies--;
    }
  }
}


void MtBin::diffuse(double t, MtBin &b) {
  //Setting up random number generator that can draw living individuals from the
  //bin with uniform probability.
  std::uniform_int_distribution<int> myguys   (0,   alive()-1);
  std::uniform_int_distribution<int> otherguys(0, b.alive()-1);

  //We are now swapping 1/10 of the population in each bin.
  unsigned int aswapn=  alive()/10;
  unsigned int bswapn=b.alive()/10;

  //Find minimum of the two swap sizes.
  int swapc=std::min(aswapn,bswapn);
  for(int i=0;i<swapc;++i)
    //Up to the shared number of swaps (swapc), swap individuals between bins
    std::swap(bin.at(myguys(rand_engine())), b.bin.at(otherguys(rand_engine())));

  //Find carrying capacity in each bin
  unsigned int ka = kkap(t);
  unsigned int kb = b.kkap(t);

  //If more individuals are leaving A than are leaving B, take the extra individuals
  //from A and move them to B.
  //Else, if more individuals are leaving B than are leaving A, take the extra individuals
  //from B and move them to A.
  if(aswapn>bswapn){
    aswapn-=swapc; //These are the excess salamanders left in A after swapping
    //Make sure that we don't exceed the carrying capacity in B.
    aswapn=std::min(aswapn,kb-b.alive());
    for(unsigned int i=0;i<aswapn;++i){
      int temp=myguys(rand_engine()); //Choose a random salamander to push into B
      b.addSalamander(bin.at(temp));
      killSalamander(temp);
    }
  } else if(aswapn<bswapn) {
    bswapn-=swapc; //These are the excess salamanders left in B after swapping
    //Make sure that we don't exceed the carrying capacity in A.
    bswapn=std::min(bswapn,ka-alive());
    for(unsigned int i=0;i<bswapn;++i){
      int temp=otherguys(rand_engine()); //Choose a random salamander to push into A
      addSalamander(b.bin.at(temp));
      b.killSalamander(temp);
    }
  }
}



///Given a time tMyrs in millions of years ago returns area at that elevation
///IN SQUARE KILOMETERS
double MtBin::area(double elevationkm, double tMyrs) const {
  ///Constants defining a normal distribution that describes area available at
  ///different height bands in the Appalachian mountains. Paramteres are fit to
  ///contemporary height distributions presented in Kozak and Wiens 2010.
  ///deltasd describes change in sd per year, which  

  //Scaling parameter transforming standard normal dist to total area of the
  //Appalachians at present (0Mya) IN SQUARE KILOMETERS, as described by the
  //calculate_area_parameters.R script
  const double elek     = 1202975;

  //Standard deviation of an analogous normal distribution to above, but 65Mya, IN KILOMETERS  
  const double elesigma = 0.3858345;

  //Mean of the above normal distribution IN KILOMETERS
  const double elemu    = 0.1455467;

  //Change per thousand years in standard deviation of between today's standard
  //deviation and the standard deviation of 65Mya
  const double deltasd  = 2.683452e-9 * 1000;

  //Input "t" is in millions of years - transform this into thousands of years
  double timeKyrs = tMyrs*1000;

  //Area of of mountain at elevation in SQUARE KILOMETERS
  double area = elek * Gaussian(elevationkm, elemu, elesigma-timeKyrs*deltasd);

  return area;
}

MtBin::container::iterator MtBin::randomSalamander(){
  assert(!bin.empty());
  std::vector<Salamander>::iterator temp=bin.begin();
  int pos=uniform_rand_int(0, bin.size()-1);
  std::advance(temp,pos);
  return temp;
}

/*
container::const_iterator MtBin::randomSalamander(){
  assert(!bin.empty());
  std::vector<Salamander>::const_iterator temp=bin.begin();
  int pos=uniform_rand_int(0, bin.size()-1);
  std::advance(temp,pos);
  return temp;
}*/