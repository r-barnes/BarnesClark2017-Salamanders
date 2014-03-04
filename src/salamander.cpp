#include "salamander.hpp"
#include "random.hpp"
#include "data.hpp"
#include <cstdlib>
#include <functional>
#include <iostream>
using namespace std;

//Counts the number of bits that are "on"
//This might use a method developed By Peter Wegner (TODO: This could be verified, but who wants to?)
//See: Communications of the ACM, Vol. 3 No. 5, Page 322
//doi: 10.1145/367236.367286
//URL: http://cacm.acm.org/magazines/1960/5/14709-a-technique-for-counting-ones-in-a-binary-computer/abstract
template<class T>
T countbits(T a){
  unsigned int c; // c accumulates the total bits set in combined
  for (c = 0; a; ++c)
    a &= a - 1; // clear the least significant bit set
  return c;
}



Salamander::Salamander(){
  genes                = 0;
  otemp                = 0;
  dead                 = false;
  parent               = -1;
  mutation_probability = 1e-4;
}

void Salamander::printGenome() const {
  Salamander::genetype selector=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    std::cerr<<((genes&selector)?'1':'0');
    selector=selector<<1;
  }
  std::cerr<<std::endl;
}

//TODO: Rethink what the normal dice are which I should be using here.
//TODO: Think about this one more time.
Salamander Salamander::breed(const Salamander &b) const {
  Salamander child;

  child.parent = parent;

  //Child optimum temperature is the average of its parents, plus a mutation,
  //drawn from a standard normal distribution with mean = 0 and sd = 0.001.
  child.otemp = (otemp+b.otemp)/2+normal_rand(0,0.001);

  //Find those genes the parents do not have in common.
  Salamander::genetype not_common_genes = ~((genes & b.genes) | (~genes & ~b.genes));

  //Child gets genes from one parent, but we'll choose some genes from the other
  //parent below.
  child.genes = genes;

  //We push the selector along one bit at a time and if that bit is one of the
  //places were the genomes differed, we choose with 50% probability whether to
  //turn it on or off.
  Salamander::genetype selector = 1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(not_common_genes&selector && rand()%2==0)
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

  return child;
}

void Salamander::mutate(){
  Salamander::genetype mutator=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(uniform_rand_real(0,1)<=mutation_probability)
      genes^=mutator;
    mutator=mutator<<1;
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
//t, and topt for this salamander. Return TRUE if salamander dies.
bool Salamander::pDie(double tempdegC) const {
  //Parameters for a logit curve, that kills a salamander with ~50% probability
  //if it is more than 8 degrees C from its optimum temperature, and with ~90%
  //probability if it is more than 12 degrees from its optimum temperature.
  const double logitslope =  0.03051701;
  const double logitint   = -2.197225;

  //For temperatures outside of these limits, the salamander always dies
  if(!(0<=tempdegC && tempdegC<=50)) {
    return true; //Dies
  } else {      
  //Find probability of death if bounds are not exceeded
    //Squared difference between salamander's optimal temp and input temp
    double dtemp = pow(otemp-tempdegC, 2);

    //Logit function, centered at f(dtemp=0)=0.1; f(dtemp=12**2)=0.9
    double pdeath = 1/(1+exp(-(dtemp*logitslope+logitint)));

    //Kill individual with probability pdeath
    if(uniform_rand_real(0,1)<pdeath)
      return true; //Dies
  }
  
  return false; //Lives
}

