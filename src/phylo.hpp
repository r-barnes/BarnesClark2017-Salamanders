#ifndef _phylo
#define _phylo

#include "salamander.hpp"
#include "mtbin.hpp"
#include <vector>
#include <string>

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
    ///Children of this node
    std::vector<int> children;
    ///Adds a child to this node
    void addChild(int n);
};

class Phylogeny {
 private:
  ///Adds a new node to the phylogeny
  void addNode(const Salamander &s, double t);
 public:
  //Empty constructor. Avoid using this whenever possible.
  Phylogeny();
  ///Initialize using a single salamander as the parent
  Phylogeny(const Salamander &s, double t);
  ///The collection of phylogenic nodes compromising the tree
  std::vector<PhyloNode> nodes;
  ///Updates the phylogeny based on the current state of the salamanders
  void UpdatePhylogeny(double t, std::vector<MtBin> &mts, double species_sim_thresh);
  ///Counts the number of species which are alive at a given point in time
  int numAlive(double t) const;
  ///Calculate the mean branch distance for the phylogeny
  typedef std::vector< std::pair<double, int> > mbdStruct;
  mbdStruct meanBranchDistance(double t) const;
  ///Calculate empirical cumulative distribution function of mean branch distances
  double compareECDF(double t) const;
  ///Print graphs of relatedness
  void print(std::string prefix) const;
  ///Prints a Newick representation of the tree's living members at a given time
  ///See: https://en.wikipedia.org/wiki/Newick_format
  std::string printNewick(int n=0, int depth=0) const;
};

#endif
