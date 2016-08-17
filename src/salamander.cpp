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
  species              = -1;
}


Salamander Salamander::breed(const Salamander &b) const {
  //The child starts out as a copy of one of its parents. We modify that copy to
  //build up the child.
  Salamander child=*this;

  //Set the child's parent species to be the parent species of one of its two
  //parents. Since both of the parents must belong to the same species, just
  //choose one.
  child.species = species;

  //Child optimum temperature is the average of its parents, plus a mutation,
  //drawn from a standard normal distribution with mean = 0 and sd = 0.001.
  child.otempdegC = (otempdegC+b.otempdegC)/2+normal_rand(0,TheParams.tempDrift());

  //Find those genes the parents do not have in common.
  Salamander::genetype not_common_genes = (genes ^ b.genes);

  //Child gets genes from one parent, but we'll choose some genes from the other
  //parent below.
  child.genes = genes;

  //We push the selector along one bit at a time and if that bit is one of the
  //places were the genomes differed, we choose with 50% probability whether to
  //turn it on or off (that is, which parent it should be inherited from).

  //For 64-bit types
  Salamander::genetype selector = uniform_bits<uint64_t>();

  //For 128-bit types, generate two 64-bit fields;
  //Salamander::genetype left_selector  = uniform_bits<uint64_t>();
  //Salamander::genetype right_selector = uniform_bits<uint64_t>();
  //Generate a 128-bit selector where each bit now has 50% probability
  //Salamander::genetype selector = (left_selector<<64) | right_selector;
  
  //Select the uncommon genes, each with 50% probability
  Salamander::genetype selected_uncommon = not_common_genes & selector;
  //We now have a bitfield, selected_uncommon, where each uncommon bit from the
  //parents is represented by a 1 with 50% probability. Now, flip those bits in
  //the child.
  child.genes ^= selected_uncommon;

  //Mutate child genome
  child.mutate();

  return child;
}


//Walk through the salamander's genome and with a probability
//`mutation_probability` flip the gene.
void Salamander::mutate(){
  std::bernoulli_distribution d(TheParams.mutationProb());
  for(unsigned int i=0;i<genes.size();i++)
    if(d(rand_engine()))
      genes[i].flip();
}


//Determine whether the genomes of two salamanders are more similar than the
//given threshold.
bool Salamander::pSimilarGenome(const Salamander::genetype &b, int species_sim_thresh) const {
  //Create genetype where only the bits which match in genes and b are on
  Salamander::genetype combined=~(genes ^ b);

  //Were enough bits shared? 
  return combined.count() >= (unsigned int)species_sim_thresh;
}


//Determine whether two salamanders are similar. This only takes genomes into
//account right now, but could, presumably, include other properties of the
//salamanders.
bool Salamander::pSimilar(const Salamander &b, int species_sim_thresh) const {
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
  //TODO: Improve documentation now that we are using params

  //For temperatures outside of these limits, the salamander always dies
  if(!(0<=tempdegC && tempdegC<=50))
    return true; //Dies

  //Find probability of death if bounds are not exceeded
  //Squared difference between salamander's optimal temp and input temp
  const double dtemp = pow(otempdegC-tempdegC, 2);

  //Logit function, centered at f(dtemp=0)=0.1; f(dtemp=12**2)=0.9
  const double pdeath = 1/(1+exp(-
    (
      TheParams.logitOffset()+dtemp*TheParams.logitTempWeight()
      +conspecific_abundance   *TheParams.logitCAweight()
      +heterospecific_abundance*TheParams.logitHAweight()
    )
  ));

  //Kill individual with probability pdeath
  return uniform_rand_real(0,1)<pdeath; //If true, salamander dies
}