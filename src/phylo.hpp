#ifndef _phylo
#define _phylo

#include "salamander.hpp"
#include "mtbin.hpp"
#include <vector>
#include <array>

class PhyloNode {
  public:
    ///Innitialize new salamandar
    PhyloNode(const Salamander &s, double t);

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

class Phylogeny {
 public:
  std::vector<PhyloNode> nodes;
  void UpdatePhylogeny(double t, std::vector<MtBin> &mts);
};

#endif
