#include "salamander.hpp"
#include "utility.hpp"

/**
  @brief Generate a random mask of genes to mutate the parent genomes
*/
Salamander::genetype RandomGeneMask(int prob){
//  for(int i=0;i<sizeof(Salamander::genetype)*8;++i)
    
}

Salamander::Salamander(){
  genes=1;
  otemp=1;
}

Salamander Salamander::breed(const Salamander &b) const {
  Salamander temp;
  temp.genes=genes ^ b.genes;
  temp.otemp= (otemp+b.otemp)/2+gaussrand();

  return temp;
}

unsigned int Salamander::similarity(const Salamander &b) const {
  Salamander::genetype combined=genes & b.genes;
  return countbits(combined);
}
