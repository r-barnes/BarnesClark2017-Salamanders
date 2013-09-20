#ifndef _phylo
#define _phylo

#include "salamander.hpp"

class Phylo {
  public:
    Phylo(Salamander::genetype &genes, double t, int parent);

    ///Genome of the strain
    Salamander::genetype genes;
    ///When this strain emerged
    double emergence;
    ///When this strain disappeared
    double extinction;
    ///Which strain this one emerged from
    int parent;
};

#endif
