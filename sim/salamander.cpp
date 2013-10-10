#include "salamander.hpp"
#include "utility.hpp"
#include <cstdlib>
#include <random>
#include <functional>
#include "data.hpp"

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
}

Salamander Salamander::breed(const Salamander &b) const {
  Salamander temp;
  temp.genes=genes & b.genes;
  temp.otemp= (otemp+b.otemp)/2+normaldice();

  Salamander::genetype selector=1;
  Salamander::genetype shared_genes=genes ^ b.genes;
  //Taking an XOR of the parent genomes will return a bitfield where all the
  //active bits represent places where one or the other of the parent genomes
  //was active, but not both. We push the selector along one bit at a time and
  //if that bit is one of the places were the genomes differed, we choose with
  //50% probability whether to turn it on or off.
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(shared_genes&selector && rand()%2==0)
      temp.genes|=selector;
    selector=selector<<1;
  } 

  return temp;
}

void Salamander::mutate(){
  Salamander::genetype mutator=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(rand()%10==0)
      genes^=mutator;
    mutator=mutator<<1;
  }
}

void Salamander::randomizeGeneome(){
  Salamander::genetype selector=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(rand()%2==0)
      genes|=selector;
    selector=selector<<1;
  }
}

bool Salamander::pSimilarGenome(const Salamander::genetype &b) const {
  Salamander::genetype combined=genes & b;
  return countbits(combined) > 98*8*sizeof(Salamander::genetype)/100;
}

bool Salamander::pSimilar(const Salamander &b) const {
  return pSimilarGenome(b.genes);
}

bool Salamander::pDie(double temp) const {
  bool dead=false;
  double pdeath; // Calculate probability of death give square distance of t, topt
  double dtemp;
  const double logitslope = 0.03051701;
  const double logitint = -2.197225;
  
  if(!(0<=temp && temp<=50)) { //Maximual temperature bounds
    dead=true;
  } else { //Calculate probability of death if bounds are not exceeded
    dtemp = pow(otemp-temp, 2);
    pdeath = 1/(1+exp(-(dtemp*logitslope+logitint))); //Logit function, centered at f(dtemp=0)=0.1; f(dtemp=12**2)=0.9
    
    if(unifdice()<pdeath) //Kill individual with probability pdeath
      dead=true;
  }
  
  return dead;
}

