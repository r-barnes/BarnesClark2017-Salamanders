#include "mtbin.hpp"
#include "temp.hpp"
#include "random.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <stdexcept>

///Generic Gaussian distribution function
double Gaussian(double x, double mean, double sigma){
  const double PI=3.14159265359;
  return 1/sigma/std::sqrt(2*PI)*std::exp(-std::pow(x-mean,2)/2/std::pow(sigma,2));
}

MtBin::MtBin(){}

MtBin::MtBin(double heightkm_val){
  this->heightkm_val = heightkm_val;
  //Reserve enough space to hold the maximum population. This keeps things
  //running fast by reducing the need to dynamically reallocate memory.
  bin.reserve(2000);
}


//Returns the nominal height of the bin.
double MtBin::heightkm() const {
  return heightkm_val;
}


///Returns the maximum height of the mountain range at the given time
double MtBin::heightMaxKm(double tMyrs) {
  if(!TheParams.pVaryHeight()) tMyrs=65;

  //Maximum elevation of the mountain range over time
  //Based on linear shrinking of mountain height from 2.8km at 65Mya (according
  //to the USGS website on "Geologic Provinces of the Untied States: Appalachian
  //Highlands Province") to current elevation (1.6km) from Kozak and Wiens 2010.
  const double height_65mya = 2.8; //km
  const double height_0mya  = 1.6; //km
  const double erosion_rate = (height_65mya-height_0mya)/65.0; //Erosion per 1Myr
  return 2.8-erosion_rate*tMyrs;   //km
}


//Removes a salamander from this bin. The removed salamander ceases to exist,
//which, for our purposes, is the same as killing it.
void MtBin::killSalamander(MtBin::container::iterator s) {
  assert(!bin.empty());

  //We overwrite the indicated salamander, which is now dead, with the
  //salamander at the back of the bin, which is still alive. If this method is
  //called by an iterator the iterator must decrement and then advance so that
  //the swapped salamander is still considered
  *s = bin.back();

  //Pop the salamander at the back of the container, because it is now a
  //duplicate and should never be accessed.
  bin.pop_back();
}


void MtBin::mortaliate(double tMyrs, int max_species, int species_sim_thresh) {
  ///If there are no living salamanders, then don't do anything
  if(bin.empty()) return;

  double mytemp = temp(tMyrs);  //Current temperature of bin

  if(alive()>30000){
    std::cerr<<"30ksals found in a bin. Killing the simulation."<<std::endl;
    throw std::runtime_error("30ksals found in a bin. Killing the simulation.");
  }

  //If individuals have the same parent species they are part of the same
  //species. Cache this here to maintain O(N) operation
  std::vector<int> species_abundance(max_species,0);
  for(const auto &s: bin)
    species_abundance.at(s.species)++;

  //For each salamander, check to see if it dies
  for(auto s=bin.begin();s!=bin.end();s++){
    //These are both initially used to count individuals. Then area is divided
    //to produce abundance.
    double conspecific_abundance    = species_abundance[s->species]-1;
    double heterospecific_abundance = bin.size()-species_abundance[s->species];

    //Turn counts into abundances, as promised
    conspecific_abundance    /= area(heightkm(), tMyrs);
    heterospecific_abundance /= area(heightkm(), tMyrs);

    //Now that we've calculated CA and HA, see if the salamander is affected by
    //it.
    if(s->pDie(mytemp, conspecific_abundance, heterospecific_abundance)){
      killSalamander(s);
      //If we kill a salamander, we swap the last living salamander in the list
      //with the salamander we just killed. Therefore, we need to make sure that
      //we still run the mortaliate function for the living salamander that now
      //inhabits the spot that we just filled.
      s--;
    }
  }
}


//Get the temperature of this bin at tMyrs, taking into account its elevation.
double MtBin::temp(double tMyrs) const {
  //Find out the temperature adjustment for that height assuming a dry air
  //adiabatic lapse rate of 9.8 degC per vertical kilometer
  double altitude_temp_adjust = -9.8*heightkm();

  return Temperature.getTemp(tMyrs) + altitude_temp_adjust;
}


//Add the indicated salamander to the bin
void MtBin::addSalamander(const Salamander &s) {
  bin.push_back(s);
}


//Return the number of living salamanders in this bin
unsigned int MtBin::alive() const {
  //Since the bin vector only ever contains living salamanders, its size is
  //equal to that number
  return bin.size();
}


//Give salamanders in this bin the opportunity to breed
void MtBin::breed(double tMyrs, int species_sim_thresh){
  if(bin.empty()) return;          //No one is alive here; there can be no breeding.

  //Maximum number of tries to find a pair to mate; prevents infinite loops.
  int maxtries = TheParams.maxTriesToBreed();

  //Maximum number of new offspring per bin per unit time
  int max_babies = TheParams.maxOffspringPerBinPerDt();

  //randomSalamaner() chooses a salamander randomly in the range [0,maxsal].
  //Baby salamanders will be added at maxsal+1, maxsal+2, ... So, by noting
  //maxsal now, we prevent baby salamanders from breeding in the timestep in
  //which they are born. This means that the only mating criteria is that two
  //salamanders be part of the same species at the beginning of the timestep.
  const int maxsal = bin.size()-1;

  //As long as there's room in the bin, and we still have to make babies, and we
  //are not caught in an infinite loop, then try to make more babies.
  while(max_babies>0 && maxtries-->0){
    auto parenta = randomSalamander(maxsal);
    auto parentb = randomSalamander(maxsal);
    //If parents are genetically similar enough to be classed as the same
    //species based on species_sim_thresh, then they can breed.
    if(parenta->species == parentb->species){
      addSalamander(parenta->breed(*parentb));
      max_babies--;
    }
  }
}


//Move salamanders from this bin to a different bin
void MtBin::moveSalamanderTo(const MtBin::container::iterator &s, MtBin &b){
  b.addSalamander(*s); //Add salamander to the indicated bin
  killSalamander(s);   //Remove the salamander from this bin
}


//Give salamanders in this bin the opportunity to move to neighbouring bins if
//advantageous
void MtBin::diffuseToBetter(double tMyrs, MtBin *lower, MtBin *upper) {
  if(bin.empty()) return;

  for(container::iterator s=bin.begin();s!=bin.end();s++){
    //Do I want to migrate?
    if(uniform_rand_real(0,1)>=TheParams.dispersalProb())
      continue;

    //Higher bins are cooler. If the salamander's optimal temperature is cooler
    //than the current bin and closer to the upper neighbour than the current
    //bin, the salamander tries to migrate up the mountain.
    if( upper
        && s->otempdegC<temp(tMyrs)
        && std::abs( s->otempdegC - upper->temp(tMyrs) ) 
                              < std::abs( s->otempdegC - temp(tMyrs) )
        && upper->heightkm()<heightMaxKm(tMyrs)
    ){
      moveSalamanderTo(s,*upper);
      --s;
    //Lower bins are warmer. If the salamander's optimal temperature is warmer
    //than the current bin and closer to the lower neighbour than the current
    //bin, the salamander tries to migrate down the mountain.
    } else if(
        lower
        && s->otempdegC>temp(tMyrs)
        && std::abs( s->otempdegC - lower->temp(tMyrs) )
                              < std::abs( s->otempdegC - temp(tMyrs) )
        && lower->heightkm()<heightMaxKm(tMyrs)
    ){
      moveSalamanderTo(s,*lower);
      --s;
    }
  }
}

//Give salamanders in this bin the opportunity to move to neighbouring bins.
void MtBin::diffuseLocal(double tMyrs, MtBin *lower, MtBin *upper) {
  if(bin.empty()) return;

  for(container::iterator s=bin.begin();s!=bin.end();s++){
    //Does the salamander want to migrate?
    if(uniform_rand_real(0,1)>=TheParams.dispersalProb())
      continue; //No

    //Am I moving up or down? Be sure not to move off the bottom or top
    if(uniform_rand_real(0,1)>0.5){
      if(upper && upper->heightkm()<heightMaxKm(tMyrs)){
        moveSalamanderTo(s,*upper);
        --s;
      }
    } else {
      if(lower && lower->heightkm()<heightMaxKm(tMyrs)){
        moveSalamanderTo(s,*lower);
        --s;
      }
    }
  }
}


//Method for moving salamanders into a special separate bin representing the
//surrounding lowlands.
void MtBin::diffuseToLowlands(MtBin &lowlands){
  for(container::iterator s=bin.begin();s!=bin.end();s++){
    //Do I want to migrate?
    if(uniform_rand_real(0,1)>=TheParams.toLowlandsProb())
      continue;

    moveSalamanderTo(s,lowlands);
    --s;
  }
}


//Method to be used by the surrounding lowlands to move salamanders back into
//the active simulation.
void MtBin::diffuseFromLowlands(MtBin &frontrange){
  for(container::iterator s=bin.begin();s!=bin.end();s++){
    //Do I want to migrate?
    if(uniform_rand_real(0,1)>=TheParams.fromLowlandsProb())
      continue;

    moveSalamanderTo(s,frontrange);
    --s;
  }
}



//Give salamanders in this bin the opportunity to move all over
void MtBin::diffuseGlobal(double tMyrs, std::vector<MtBin> &mts) {
  if(bin.empty()) return;

  for(container::iterator s=bin.begin();s!=bin.end();s++){
    //Do I want to migrate?
    if(uniform_rand_real(0,1)>=TheParams.dispersalProb())
      continue;

    int to_bin = -1;
    //Choose a bin to migrate to. Loop until the chosen bin is valid, in the
    //sense of not being above the top of the mountain.
    while(to_bin==-1 || mts[to_bin].heightkm() >= heightMaxKm(tMyrs))
      to_bin = uniform_rand_int(0,mts.size()-1);

    moveSalamanderTo(s,mts[to_bin]);
    --s;
  }
}



///Given a time tMyrs in millions of years ago returns area at that elevation
///IN SQUARE KILOMETERS
double MtBin::area(double elevationkm, double tMyrs) const {
  if(!TheParams.pVaryHeight()) tMyrs=65;

  ///Constants defining a normal distribution that describes area available at
  ///different height bands in the Appalachian mountains. Parameters are fit to
  ///contemporary height distributions presented in Kozak and Wiens 2010.
  ///deltasd describes change in sd per year, which

  //Scaling parameter transforming standard normal dist to total area of the
  //Appalachians at present (0Mya) IN SQUARE KILOMETERS, as described by the
  //calculate_area_parameters.R script
  const double elek     = 1202975;

  //Standard deviation of an analogous normal distribution to above, but 65Mya,
  //IN KILOMETERS
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

void MtBin::killAll() {
  for(auto s=bin.begin();s!=bin.end();s++){
    killSalamander(s);
    //If we kill a salamander, we swap the last living salamander in the list
    //with the salamander we just killed. Therefore, we need to make sure that
    //we still run the mortaliate function for the living salamander that now
    //inhabits the spot that we just filled.
    s--;
  }
}


//Choose a random salamander [0,maxsal] from the bin
MtBin::container::iterator MtBin::randomSalamander(int maxsal){
  //Cannot run this on an empty bin
  assert(!bin.empty());
  //Generate an iterator to the beginning of the bin
  std::vector<Salamander>::iterator temp = bin.begin();
  //Choose a random member of the bin
  int pos=uniform_rand_int(0, maxsal);
  //Advance the iterator so that it points at this member
  std::advance(temp,pos);
  //Return the iterator
  return temp;
}