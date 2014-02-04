#ifndef _phylo
#define _phylo

#include "salamander.hpp"
#include "mtbin.hpp"
#include <vector>

class PhyloNode {
  public:
    ///Innitialize new salamandar
    PhyloNode(const Salamander &s, double t);
    void addChild(int childNode);

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
    ///Which nodes store my children?
    std::vector<int> children;
};

class Phylogeny {
 private:
  ///Adds a new node to the phylogeny
  void addNode(const Salamander &s, double t);
 public:
  ///Initialize using a single salamander as the parent
  Phylogeny(const Salamander &s, double t);
  ///The collection of phylogenic nodes compromising the tree
  std::vector<PhyloNode> nodes;
  ///Updates the phylogeny based on the current state of the salamanders
  void UpdatePhylogeny(double t, std::vector<MtBin> &mts);
  ///Counts the number of species which are alive at a given point in time
  int numAlive(double t) const;
  ///Calculate the mean branch distance for the phylogeny
  typedef std::vector< std::pair<double, int> > mbdStruct;
  mbdStruct meanBranchDistance(double t) const;
  ///Calculate empirical cumulative distribution function of mean branch distances
  void getECDF(const Phylogeny &p, double t) const;
};

#endif
