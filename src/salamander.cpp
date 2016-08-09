#include "salamander.hpp"
#include "random.hpp"
#include "params.hpp"
#include <cstdlib>
#include <functional>
#include <iostream>
using namespace std;

Salamander::Salamander(){
  genes                = 0;
  otempdegC            = 0;
  parent               = -1;
  mutation_probability = 1e-4;
  temperature_drift_sd = 1e-3;
}


Salamander Salamander::breed(const Salamander &b) const {
  //The child starts out as a copy of one of its parents. We modify that copy to
  //build up the child.
  Salamander child=*this;

  //Set the child's parent species to be the parent species of one of its two
  //parents. Since both of the parents must belong to the same species, just
  //choose one.
  child.parent = parent;

  //Child optimum temperature is the average of its parents, plus a mutation,
  //drawn from a standard normal distribution with mean = 0 and sd = 0.001.
  child.otempdegC = (otempdegC+b.otempdegC)/2+normal_rand(0,temperature_drift_sd);

  //Find those genes the parents do not have in common.
  Salamander::genetype not_common_genes = ~((genes & b.genes) | (~genes & ~b.genes));

  //Child gets genes from one parent, but we'll choose some genes from the other
  //parent below.
  child.genes = genes;

  //We push the selector along one bit at a time and if that bit is one of the
  //places were the genomes differed, we choose with 50% probability whether to
  //turn it on or off (that is, which parent it should be inherited from).
  Salamander::genetype selector = 1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(not_common_genes&selector && rand()%2==0)
      child.genes|=selector;
    selector=selector<<1;
  }

  //Mutate child genome
  child.mutate();

  return child;
}


//Walk through the salamander's genome and with a probability
//`mutation_probability` flip the gene.
void Salamander::mutate(){
  Salamander::genetype mutator=1;
  for(unsigned int i=0;i<sizeof(Salamander::genetype)*8;++i){
    if(uniform_rand_real(0,1)<=TheParams::get().mutationProb())
      genes^=mutator;
    mutator=mutator<<1;
  }
}


//Determine whether the genomes of two salamanders are more similar than the
//given threshold.
bool Salamander::pSimilarGenome(const Salamander::genetype &b, double species_sim_thresh) const {
  Salamander::genetype combined=(genes & b) | (~genes & ~b); //TODO: What does this do?
  unsigned int shared_genes=__builtin_popcountll(combined);  //Counts the number of 1 bits NOTE: Ensure that the popcount type matches the genetype
  return shared_genes > species_sim_thresh*8*sizeof(Salamander::genetype);
}


//Determine whether two salamanders are similar. This only takes genomes into
//account right now, but could, presumably, include other properties of the
//salamanders.
bool Salamander::pSimilar(const Salamander &b, double species_sim_thresh) const {
  return pSimilarGenome(b.genes, species_sim_thresh);
}


//Calculate probability of death given square distance between temperature at
//time t, and topt for this salamander. Return TRUE if salamander dies.
bool Salamander::pDie(
  const double tempdegC,
  const double conspecific_abundance,
  const double heterospecific_abundance
) const {
  //Parameters for a logit curve, that kills a salamander with ~50% probability
  //if it is more than 8 degrees C from its optimum temperature, and with ~90%
  //probability if it is more than 12 degrees from its optimum temperature.
  //TODO: Or it did!

  //For temperatures outside of these limits, the salamander always dies
  if(!(0<=tempdegC && tempdegC<=50))
    return true; //Dies

  //Find probability of death if bounds are not exceeded
  //Squared difference between salamander's optimal temp and input temp
  const double dtemp = pow(otempdegC-tempdegC, 2);

  //Logit function, centered at f(dtemp=0)=0.1; f(dtemp=12**2)=0.9
  const double pdeath = 1/(1+exp(-
    (
      TheParams::get().logitOffset()+dtemp*TheParams::get().logitTempWeight()
      +conspecific_abundance   *TheParams::get().logitCAweight()
      +heterospecific_abundance*TheParams::get().logitHAweight()
    )
  ));

  //Kill individual with probability pdeath
  return uniform_rand_real(0,1)<pdeath; //If true, salamander dies
}