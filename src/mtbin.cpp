#include "mtbin.hpp"
#include "temp.hpp"
#include "random.hpp"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iterator>
#include <iostream>
#include <iomanip>


///Generic Gaussian distribution function
double Gaussian(double x, double mean, double sigma){
  const double PI=3.14159265359;
  return 1/sigma/std::sqrt(2*PI)*std::exp(-std::pow(x-mean,2)/2/std::pow(sigma,2));
}


MtBin::MtBin(double heightkm0, bool vary_height0){
  heightkm    = heightkm0;
  vary_height = vary_height0;
  bin.reserve(binmax); //Maximum salamander populations per mountain bin
}


double MtBin::height() const {
  return heightkm;
}


void MtBin::killSalamander(MtBin::container::iterator s) {
  assert(!bin.empty());

  //We swap the indicated salamander, which is now dead, with the salamander at
  //the back of the bin, which is still alive. If this method is called by an
  //iterator the iterator must decrement and then advance so that the swapped
  //salamander is still considered
  std::swap(*s,bin.back());

  //Pop the dead salamander out of the container, so that it can never be
  //accessed.
  bin.pop_back();
}


void MtBin::mortaliate(double tMyrs) {
  ///If there are no living salamanders, then don't do anything
  if(bin.empty()) return;

  double       mytemp   = temp(tMyrs);  //Current temperature of bin
  unsigned int maxalive = kkap(tMyrs);  //Current carrying capacity of the bin

  //Since the carrying capacity of the bin may have been reduced or exceeded
  //since we last mortaliated, kill individuals at random until we are within
  //carrying capacity.
  while(alive()>maxalive)
    killSalamander(randomSalamander());

  //For each salamander, check to see if it dies
  for(auto s=bin.begin();s!=bin.end();s++)
    if(s->pDie(mytemp)){
      killSalamander(s);
      //If we kill a salamander, we swap the last living salamander in the list
      //with the salamander we just killed. Therefore, we need to make sure that
      //we still run the mortaliate function for the living salamander that now
      //inhabits the spot that we just filled.
      s--;
    }
}


double MtBin::temp(double tMyrs) const {
  //Find out the temperature adjustment for that height assuming a dry air
  //adiabatic lapse rate of 9.8 degC per vertical kilometer
  double altitude_temp_adjust = -9.8*heightkm;

  return Temperature::getInstance().getTemp(tMyrs) + altitude_temp_adjust;
}


unsigned int MtBin::kkap(double tMyrs) const {
  if(!vary_height) tMyrs=65;

  //Input "t" is in millions of years - transform this into thousands of years
  double timeKyrs = tMyrs*1000;

  //NOTE: Area isn't integrated over the height of the band. This is the
  //simplest way of handling area, but could be refined. However, we do not
  //expect such a refinement would alter our conclusions.

  //Maximum elevation of the mountain range over time
  //Based on linear shrinking of mountain hight from 2.8km at 65Mya (according
  //to the USGS website on "Geologic Provinces of the Untied States: Appalachian
  //Highlands Province") to current elevation (1.6km) from Kozak and Wiens 2010.
  double const height_65mya = 2.8; //km
  double const height_0mya  = 1.6; //km
  double const erosion_rate = (height_65mya-height_0mya)/65000; //Erosion rate per 1kyr
  double maxelevation       = 2.8-erosion_rate*timeKyrs; //km
  double minarea_today      = area(maxelevation, tMyrs);

  //Returns a number [0, binmax].The smallest area (at the top of the mountain)
  //will always have a carrying capacity of at least 1 salamander and all other
  //bins are scaled to this bin's size.
  return std::min( area(heightkm, tMyrs)/minarea_today, (double) binmax);
}


void MtBin::addSalamander(const Salamander &s) {
  bin.push_back(s);
}


unsigned int MtBin::alive() const {
  return bin.size();
}


void MtBin::breed(double t, double species_sim_thresh){
  if(bin.empty()) return;          //No one is alive here; there can be no breeding.

  unsigned int maxalive=kkap(t);   //Current carrying capacity of the bin
  if(alive()>=maxalive) return;    //The bin is too full for us to breed

  //Maximum number of tries to find a pair to mate; prevents infinite loops.
  int maxtries=40*(maxalive-alive());

  //Maximum number of new offspring per bin per unit time
  int max_babies=10;

  //As long as there's room in the bin, and we still have to make babies, and we
  //are not caught in an infinite loop, then try to make more babies.
  while(alive()<maxalive && max_babies>0 && maxtries-->0){
    auto parenta=randomSalamander();
    auto parentb=randomSalamander();
    //If parents are genetically similar enough to be classed as the same
    //species based on species_sim_thresh, then they can breed.
    if(parenta->pSimilar(*parentb, species_sim_thresh)){
      addSalamander(parenta->breed(*parentb));
      max_babies--;
    }
  }
}

void MtBin::moveSalamanderTo(const MtBin::container::iterator &s, MtBin &b){
  b.addSalamander(*s);
  killSalamander(s);
}

void MtBin::diffuse(double t, MtBin *lower, MtBin *upper) {
  if(bin.empty()) return;

  //The following code allows salamanders to move into neighbouring bins in a
  //way which may exceed their carrying capacity; however, in the Mortaliate()
  //step, random salamanders are killed until the bin is back at the carrying
  //capacity. Salamanders move without respecting carrying capacities and nature
  //does not choose who wins and loses in this process.

  for(container::iterator s=bin.begin();s!=bin.end();s++){
    //10% chance of wanting to migrate
    if(uniform_rand_real(0,1)>=0.1)
      continue;

    //Higher bins are cooler. If the salamander's optimal temperature is cooler
    //than the current bin and closer to the upper neighbour than the current
    //bin, the salamander tries to migrate up the mountain.
    if( upper && 
        s->otempdegC<temp(t) &&
        std::abs(s->otempdegC-upper->temp(t)) < std::abs(s->otempdegC-temp(t))
    ){
      moveSalamanderTo(s,*upper);
      --s;
    //Lower bins are warmer. If the salamander's optimal temperature is warmer
    //than the current bin and closer to the lower neighbour than the current
    //bin, the salamander tries to migrate down the mountain.
    } else if(lower && 
              s->otempdegC>temp(t) && 
              std::abs(s->otempdegC-lower->temp(t)) < std::abs(s->otempdegC-temp(t))
    ){
      moveSalamanderTo(s,*lower);
      --s;
    }
  }
}



///Given a time tMyrs in millions of years ago returns area at that elevation
///IN SQUARE KILOMETERS
double MtBin::area(double elevationkm, double tMyrs) const {
  if(!vary_height) tMyrs=65;

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


//This runs a series of tests designed to ensure that the simulation is behaving
//as we expect.
void MtBinUnitTest::run() const {
  unsigned int num_to_test=100;

  {
    std::cerr<<"MtBin: Simple addition and removal test."<<std::endl;
    std::cerr<<std::endl;
    MtBin bin(0, true);
    for(unsigned int i=0;i<num_to_test;i++)
      bin.addSalamander(Salamander());
    for(unsigned int i=0;i<num_to_test;i++)
      bin.killSalamander(bin.bin.begin());
    assert(bin.alive()==0);
    std::cerr<<"Passed"<<std::endl;
  }

  {
    std::cerr<<"MtBin: Simple addition and random removal test."<<std::endl;
    std::cerr<<std::endl;
    MtBin bin(0,true);
    for(unsigned int i=0;i<num_to_test;i++)
      bin.addSalamander(Salamander());
    for(unsigned int i=0;i<num_to_test;i++)
      bin.killSalamander(bin.randomSalamander());
    assert(bin.alive()==0);
    std::cerr<<"Passed"<<std::endl;
  }

  {
    std::cerr<<std::endl;
    std::cerr<<"MtBin: Test whether diffusion UP works."<<std::endl;
    MtBin bin1(0,true), bin2(2.8/10,true); //A bin at sea-level and a bin just above it
    std::cerr<<"Bin1 temp="<<bin1.temp(0)<<". Bin2 temp="<<bin2.temp(0)<<std::endl;
    std::cerr<<"Remaining salamanders: ";
    //Fill bin1 with a number of salamanders optimally adapted to bin2
    for(unsigned int i=0;i<num_to_test;i++){
      Salamander temp;
      temp.otempdegC=bin2.temp(0);
      bin1.addSalamander(temp);
    }
    //Test 100 times to show effect of diffusion
    for(unsigned int i=0;i<100;i++){
      std::cerr<<" "<<bin1.alive();
      bin1.diffuse(0,nullptr,&bin2);
    }
    std::cerr<<std::endl;
  }

  {
    std::cerr<<std::endl;
    std::cerr<<"MtBin: Test whether diffusion DOWN works."<<std::endl;
    MtBin bin1(2.8/10, true), bin2(0, true); //A bin at sea-level and a bin just above it
    std::cerr<<"Bin1 temp="<<bin1.temp(0)<<". Bin2 temp="<<bin2.temp(0)<<std::endl;
    std::cerr<<"Remaining salamanders: ";
    //Fill bin1 with a number of salamanders optimally adapted to bin2
    for(unsigned int i=0;i<num_to_test;i++){
      Salamander temp;
      temp.otempdegC=bin2.temp(0);
      bin1.addSalamander(temp);
    }
    //Test 100 times to show effect of diffusion
    for(unsigned int i=0;i<100;i++){
      std::cerr<<" "<<bin1.alive();
      bin1.diffuse(0,&bin2,nullptr);
    }
    std::cerr<<std::endl;
  }

  {
    std::cerr<<std::endl;
    std::cerr<<"MtBin: Test whether diffusion both ways works."<<std::endl;
    //A bin at sea-level and a bin just above it
    MtBin bin1(2.8/10, true), bin2(0, true), bin3(2.8/10*2, true); 
    std::cerr<<"Bin1 temp="<<bin1.temp(0)<<". Bin2 temp=";
    std::cerr<<bin2.temp(0)<<". Bin3 temp="<<bin3.temp(0)<<"."<<std::endl;
    std::cerr<<"Output of type: bin1(bin2,bin3)"<<std::endl;
    std::cerr<<"Remaining salamanders: ";
    //Fill bin1 with a number of salamanders optimally adapted to bin2
    for(unsigned int i=0;i<num_to_test;i++){
      Salamander temp;
      temp.otempdegC=bin2.temp(0);
      bin1.addSalamander(temp);
    }
    //Fill bin1 with a number of salamanders optimally adapted to bin3
    for(unsigned int i=0;i<num_to_test;i++){
      Salamander temp;
      temp.otempdegC=bin3.temp(0);
      bin1.addSalamander(temp);
    }
    //Test 100 times to show effect of diffusion
    for(unsigned int i=0;i<100;i++){
      std::cerr<<" "<<bin1.alive()<<"("<<bin2.alive()<<","<<bin3.alive()<<") ";
      bin1.diffuse(0,&bin2,&bin3);
    }
    std::cerr<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Test salamander freezing death response"<<std::endl;
    Temperature::getInstance().testOn(-10);
    for(unsigned int i=0;i<num_to_test;i++)
      bin.addSalamander(Salamander());
    bin.mortaliate(0);
    assert(bin.alive()==0);
    std::cerr<<"Passed"<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Test salamander heat death response"<<std::endl;
    Temperature::getInstance().testOn(100);
    for(unsigned int i=0;i<num_to_test;i++)
      bin.addSalamander(Salamander());
    bin.mortaliate(0);
    assert(bin.alive()==0);
    std::cerr<<"Passed"<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Test salamander heat death response. Global temp=33C. Salamander=25C.\n";
    Temperature::getInstance().testOn(33);
    for(unsigned int i=0;i<num_to_test;i++){
      Salamander temp;
      temp.otempdegC=25;
      bin.addSalamander(temp);
    }

    bin.mortaliate(0);
    std::cerr<<"Salamanders left alive: "<<bin.alive()<<"/"<<num_to_test<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Test salamander heat death response. Global temp=25C. Salamander=33C.\n";
    Temperature::getInstance().testOn(25);
    for(unsigned int i=0;i<num_to_test;i++){
      Salamander temp;
      temp.otempdegC=33;
      bin.addSalamander(temp);
    }

    bin.mortaliate(0);
    std::cerr<<"Salamanders left alive: "<<bin.alive()<<"/"<<num_to_test<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Breeding a sea-level bin at t=0 with "<<100<<" salamanders leaves ";
    for(unsigned int i=0;i<100;i++)
      bin.addSalamander(Salamander());
    bin.breed(0,0);
    std::cerr<<bin.alive()<<" salamanders afterwards."<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Breeding a sea-level bin at t=0 with "<<200<<" salamanders leaves ";
    for(unsigned int i=0;i<200;i++)
      bin.addSalamander(Salamander());
    bin.breed(0,0);
    std::cerr<<bin.alive()<<" salamanders afterwards."<<std::endl;
  }

  {
    MtBin bin(0, true);
    std::cerr<<std::endl;
    std::cerr<<"Breeding a sea-level bin at t=0 with "<<5<<" salamanders leaves ";
    for(unsigned int i=0;i<5;i++)
      bin.addSalamander(Salamander());
    bin.breed(0,0);
    std::cerr<<bin.alive()<<" salamanders afterwards."<<std::endl;
  }

  {
    //Construct bins at 0km, 1.6km, and 2.8km
    MtBin bin0(0, true), bin1(1.6, true), bin2(2.8, true);
    Temperature::getInstance().testOff();
    std::cerr<<std::endl;
    std::cerr<<"Temperature at 0Myr at: "<<std::endl;
    std::cerr<<"  0.0km="<<bin0.temp(0)<<std::endl;
    std::cerr<<"  1.6km="<<bin1.temp(0)<<std::endl;
    std::cerr<<"  2.8km="<<bin2.temp(0)<<std::endl;
    std::cerr<<"Temperature at 64.9Myr at: "<<std::endl;
    std::cerr<<"  0.0km="<<bin0.temp(64.9)<<std::endl;
    std::cerr<<"  1.6km="<<bin1.temp(64.9)<<std::endl;
    std::cerr<<"  2.8km="<<bin2.temp(64.9)<<std::endl;
  }

  {
    MtBin bin_novary(2.8,false), bin_vary(2.8,true);
    std::cerr<<std::endl;
    std::cerr<<"Checking kkap and area for mountain with non-varying height."<<std::endl;
    std::cerr<<std::setw(10)<<"kkap "<<std::setw(10)<<"area"<<std::endl;
    for(double tMyrs=0;tMyrs<65.000;tMyrs+=1)
      std::cerr<<std::setw(10)<<bin_novary.kkap(tMyrs)<<" "<<std::setw(10)
               <<bin_novary.area(2.8, tMyrs)<<std::endl;

    std::cerr<<std::endl;
    std::cerr<<"Checking area and kkap for mountain with varying height."<<std::endl;
    std::cerr<<std::setw(10)<<"kkap "<<std::setw(10)<<"area"<<std::endl;
    for(double tMyrs=0;tMyrs<65.000;tMyrs+=1)
      std::cerr<<std::setw(10)<<bin_vary.kkap(tMyrs)<<" "<<std::setw(10)
               <<bin_vary.area(2.8, tMyrs)<<std::endl;
  }
}
