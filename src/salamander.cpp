#include "salamander.hpp"
#include "utility.hpp"
#include <cstdlib>
#include <random>
#include <functional>
#include "data.hpp"
#include <iostream>
using namespace std;

std::default_random_engine generator;
std::normal_distribution<double> normdistribution(0,0.001);
std::uniform_real_distribution<double> unifdistribution(0.0,1.0);
auto normaldice=std::bind(normdistribution,generator);
auto unifdice=std::bind(unifdistribution,generator);

Salamander::Salamander(){
  genes=0;
  otemp=0;
  dead=false;
  parent=-1;
  mutation_probability=1e-4;
}

void Salamander::printGenome() const {
  Salamander::genetype selector=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    std::cerr<<((genes&selector)?'1':'0');
    selector=selector<<1;
  }
  std::cerr<<std::endl;
}

Salamander Salamander::breed(const Salamander &b) const {
  Salamander child;
  child.genes=genes & b.genes;
  //Child optimum temperature is the average of its parents, plus a mutation,
  //drawn from a standard normal distribution with mean = 0 and sd = 0.001.
  child.otemp= (otemp+b.otemp)/2+normaldice();

  //Combine genomes of parents. This results in a child with bitfield that matches
  //its parents wherever their bitfields match, and is chosen randomly from one
  //of its parents where they do not match.
  Salamander::genetype selector=1;
  Salamander::genetype shared_genes=genes ^ b.genes;
  //Taking an XOR of the parent genomes will return a bitfield where all the
  //active bits represent places where one or the other of the parent genomes
  //was active, but not both. We push the selector along one bit at a time and
  //if that bit is one of the places were the genomes differed, we choose with
  //50% probability whether to turn it on or off.
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(shared_genes&selector && rand()%2==0)
      child.genes|=selector;
    selector=selector<<1;
  }

  //Mutate child genome
  child.mutate();
/*
  printGenome();
  cerr<<"Bits on="<<countbits(genes)<<endl;
  b.printGenome();
  temp.printGenome();
*/
  child.parent=parent;

  return child;
}

void Salamander::mutate(){
  Salamander::genetype mutator=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(unifdice()<=mutation_probability){
//      std::cerr<<"GAAAH! MUTATION"<<std::endl;
      genes^=mutator;
    }
    mutator=mutator<<1;
  }
}

void Salamander::randomizeGenome(){
  Salamander::genetype selector=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(rand()%2==0)
      genes|=selector;
    selector=selector<<1;
  }
}

bool Salamander::pSimilarGenome(const Salamander::genetype &b, double species_sim_thresh) const {
  Salamander::genetype combined=(genes & b) | (~genes & ~b);
  unsigned int shared_genes=countbits(combined);
  return shared_genes > species_sim_thresh*8*sizeof(Salamander::genetype);
}

bool Salamander::pSimilar(const Salamander &b, double species_sim_thresh) const {
  return pSimilarGenome(b.genes, species_sim_thresh);
}

//Calculate probability of death given square distance between temperature at time
//t, and topt for this salamander
bool Salamander::pDie(double temp) const {
  //start with a living salamader
  bool dead=false;
  //parameters for a logit curve, that kills a salamander with ~50% probability
  //if it is more than 8 degrees C from its optimum temperature, and with ~90%
  //probability if it is more than 12 degrees from its optimum temperature.
  const double logitslope =  0.03051701;
  const double logitint   = -2.197225;
  //cerr<<"temp = "<<temp<<endl;
  if(!(0<=temp && temp<=50)) { //Temperature bounds
    //for temperatures outside of these limits (<0C and >50C), salamaders always die
    dead=true;
        //cerr<<"we have killed with 100-percent prob.!"<<endl;
  } else {                     //Find probability of death if bounds are not exceeded
    double dtemp  = pow(otemp-temp, 2); //Squared distance between temp and otemp

    //Logit function, centered at f(dtemp=0)=0.1; f(dtemp=12**2)=0.9
    double pdeath = 1/(1+exp(-(dtemp*logitslope+logitint)));
    //cerr<<"pdeath = "<<pdeath<<endl;
    //cerr<<"otemp = "<<otemp<<endl;
    //cerr<<"temp = "<<temp<<endl;
    if(unifdice()<pdeath)      //Kill individual with probability pdeath
      dead=true;
  }
  
  return dead;
}

