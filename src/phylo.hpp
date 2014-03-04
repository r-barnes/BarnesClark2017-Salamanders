#ifndef _phylo
#define _phylo

#include "salamander.hpp"
#include "mtbin.hpp"
#include <vector>
#include <string>

///PhyloNode is used to store information about distinct species, the time that
///the node came into being, the time that it does extict, its genetic attributes
///the parent species from which the node arose, and any child species
///that arise from the node. This is used to figure out which species
///and lineage each salamander belongs to. 
class PhyloNode {
  public:
    ///Innitialize new salamandar in the phylogeny.
    PhyloNode(const Salamander &s, double t);
    ///Genome of the strain. Inherited from salamandar::Salamander().
    Salamander::genetype genes;
    ///When this strain emerged. emergence describes the time that
    ///the species arises, and is taken as the current t.
    double emergence;
    ///Last time this strain had a child. lastchild is used to identify species
    ///that have gone extinct, and is updated with the current time t every
    ///time that a living salamader is identified a member of the
    ///species recorded by PhyloNode.
    double lastchild;
    ///Which strain this one emerged from? Inherited from salamandar::Salamander().
    int parent;
    ///Which otemp did the first parent inherit? Inherited from salamandar::Salamander().
    double otemp;
    ///Children of this node. Lists all species that have branched off of the
    ///species described by PhyloNode.
    std::vector<int> children;
    ///Adds a child to this node. Child is represented by an integer, which lists
    ///all species as sequential integers starting at 0, based on their time of
    ///emergence.
    void addChild(int n);
};

///Phylogeny stores all of the nodes corresponding to living and extinct species.
///It also includes functions for making comparing generated phylogenies to the
///results in Kozak and Wiens 2010, as well as make .tre files in Newick format.
class Phylogeny {
 private:
  ///Adds a new node to the phylogeny
  void addNode(const Salamander &s, double t);
 public:
  ///Empty constructor -- creates a phylogeny without any attributes.
  ///Avoid using this whenever possible!!
  Phylogeny();
  ///Initialize using a single salamander as the parent
  Phylogeny(const Salamander &s, double t);
  ///The collection of phylogenic nodes compromising the tree
  std::vector<PhyloNode> nodes;
  ///Updates the phylogeny based on the current state of the salamanders
  void UpdatePhylogeny(double t, std::vector<MtBin> &mts, double species_sim_thresh);
  ///Counts the number of species which are alive at a given point in time
  int numAlive(double t) const;
  ///Calculate the mean branch distance for the phylogeny. Finds the number of
  ///years separating them from their last common ancestor (e.g., if 2MY
  ///separates each from the ancestor, distance = 4MY), and then takes the
  ///average of this number across all unique species pairs. Is used as a summary
  ///statistic for comparing to the Kozak and Wiens phylogeny.
  typedef std::vector< std::pair<double, int> > mbdStruct;
  mbdStruct meanBranchDistance(double t) const;
  ///Calculate empirical cumulative distribution function of branch distances.
  ///Creates evenly-spaced bins along the range of branch lenghts that result
  ///from the simulation, and finds the cumulative number of species pair for
  ///which branch length falls at or below the value represented by the bin.
  ///Compares this distribution to the observed distribution from the
  ///phylogeny in Kozak and Wiens 2010 to assess phylogeny similarity.
  double compareECDF(double t) const;
  ///Print graphs of relatedness
  void print(std::string prefix) const;
  ///Prints a Newick representation of the tree's living members at a given time
  ///See: https://en.wikipedia.org/wiki/Newick_format
  std::string printNewick(int n=0, int depth=0) const;
};

#endif
