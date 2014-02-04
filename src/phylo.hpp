#ifndef _phylo
#define _phylo

#include "salamander.hpp"
#include "mtbin.hpp"
#include <vector>
#include <array>

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

typedef std::vector<Phylo> phylolist;

void UpdatePhylogeny(double t, std::vector<MtBin> &mts, phylolist &plist);

#endif
