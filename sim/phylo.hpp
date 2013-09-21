#ifndef _phylo
#define _phylo

#include "salamander.hpp"

class Phylo {
  public:
    ///Innitialize new salamandar
    Phylo(const Salamander &s, double t);

    ///Genome of the strain
    Salamander::genetype genes;
    ///When this strain emerged
    double emergence;
    ///Last time this strain had a child
    double lastchild;
    ///Which strain this one emerged from
    int parent;
    ///Which otemp did the first parent inherit?
    double otemp;
};

#endif
