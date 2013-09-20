#include "salamander.hpp"
#include "utility.hpp"
#include <cstdlib>

Salamander::Salamander(){
  genes=1;
  otemp=1;
}

Salamander Salamander::breed(const Salamander &b) const {
  Salamander temp;
  temp.genes=genes & b.genes;
  temp.otemp= (otemp+b.otemp)/2+gaussrand();

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

unsigned int Salamander::similarity(const Salamander &b) const {
  Salamander::genetype combined=genes & b.genes;
  return countbits(combined);
}
