#ifndef _phylo
#define _phylo

#include "salamander.hpp"
#include "mtbin.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <limits>
#include <algorithm>

//SpeciesStats is used to hold summary statistics above the distribution of
//salamander properties at each time step of a species' existence
class SpeciesStats {
 public:
  double t;             //Time at which these statistics were collected
  int num_alive;        //Number of species alive at that time
  double elev_min;      //Minimum elevation of the species at that time
  double elev_max;      //Maximum elevation of the species at that time
  double elev_avg;      //Average elevation of the species at that time
  double opt_temp_min;  //Minimum of the optimum temperature trait at that time
  double opt_temp_max;  //Maximum of the optimum temperature trait at that time
  double opt_temp_avg;  //Average of the optimum temperature trait at that time
  SpeciesStats(double t0){
    //In the following we set everything up so it is ready for a "reduction"
    //pattern
    t            = t0;
    num_alive    = 0;
    elev_min     = std::numeric_limits<double>::max();
    elev_max     = std::numeric_limits<double>::min();
    elev_avg     = 0;
    opt_temp_min = std::numeric_limits<double>::max();
    opt_temp_max = std::numeric_limits<double>::min();
    opt_temp_avg = 0;
  }

  void update(double elevation, double opt_temp) {
    num_alive++;
    elev_min     = std::min(elev_min,elevation);
    elev_max     = std::max(elev_max,elevation);
    elev_avg    += elevation;
    opt_temp_min = std::min(opt_temp_min,opt_temp);
    opt_temp_max = std::max(opt_temp_max,opt_temp);
    opt_temp_avg += opt_temp;
  }
};


///PhyloNode is used to store information about distinct species, the time that
///the node came into being, the time that it went extinct, the genetic
///attributes of the parent species from which the node arose, and any child
///species that arose from the node. This is used to figure out which species
///and lineage each salamander belongs to.
class PhyloNode {
 public:
  ///Copy relevant properties from the indicated salamander into this
  ///phylogenetic record.
  PhyloNode(const Salamander &s, double t);

  ///Genome of the strain.
  Salamander::genetype genes;

  ///When this strain emerged. Emergence describes the time that the species
  ///arose. Copied from Salamander on initialisation.
  double emergence;

  ///Last time this strain had a child. lastchild is used to identify species
  ///that have gone extinct, and is updated with the current time whenever a
  ///living salamander is identified as a member of the species.
  double lastchild;

  ///Which strain this one emerged from. Copied from Salamander on
  ///initialisation.
  int parent;

  ///Which optimal temperature did the first member of this strain have?
  ///Copied from Salamander on initialisation.
  double otempdegC;

  ///Children of this node. Lists all species that have branched off of this
  ///species.
  std::vector<int> children;

  std::vector<SpeciesStats> stats;

  ///Adds a child to this node. Child is represented by an integer, which
  ///corresponds to the child's placement in the Phylogeny class's list of
  ///Phylonodes.
  void addChild(int n);

  ///Returns true if the strain is alive at the indicated time, based on the
  ///emergence and lastchild data.
  bool aliveAt(double t) const;

  ///Sets the lastchild time and updates the species' statistics
  void updateWithSal(const MtBin &mt, const Salamander &s, double t);
};


///Phylogeny stores all of the nodes corresponding to living and extinct
///species. It also includes functions for comparing generated phylogenies to
///the results in Kozak and Wiens 2010, as well as to produce .tre files in
///Newick format.
class Phylogeny {
 private:
  ///Adds a new node to the phylogeny
  void addNode(const Salamander &s, double t);

  ///Calculate the mean branch distance for the phylogeny. Finds the
  //distance between each species and the last common ancestor of that species
  //and all other species in the phylogeny. (e.g., if 2MY
  ///separates each from the ancestor, distance = 4MY), and then takes the
  ///average of this number across all unique species pairs. Is used as a
  ///summary statistic for comparing to the Kozak and Wiens phylogeny.
  ///For each species i. Examine every other species j. Find distance between
  ///species i and last common ancestor of i and j. Average these distances.
  ///mbdStruct is a <Mean Branch Length, Species ID> pair
  typedef std::vector< std::pair<double, int> > mbdStruct;
  mbdStruct meanBranchDistance(double t) const;

 public:
  ///Empty constructor -- creates a phylogeny without any attributes.
  ///Avoid using this whenever possible!!
  Phylogeny();

  ///Initialize using a single salamander as the parent
  Phylogeny(const Salamander &s, double t);

  ///The collection of phylogenetic nodes compromising the tree
  std::vector<PhyloNode> nodes;

  ///Updates the phylogeny based on the current state of the salamanders
  void UpdatePhylogeny(double t, double dt, std::vector<MtBin> &mts);

  ///Counts the number of species which are alive at a given point in time
  int livingSpecies(double t) const;

  ///Calculate empirical cumulative distribution function of branch distances.
  ///Creates evenly-spaced bins along the range of branch lengths that result
  ///from the simulation, and finds the cumulative number of species pair for
  ///which branch length falls at or below the value represented by the bin.
  ///Compares this distribution to the observed distribution from the phylogeny
  ///in Kozak and Wiens 2010 to assess phylogeny similarity.
  double compareECDF(double t) const;

  ///Print species labels and their persistence to the specified output stream
  void persistGraph(std::ofstream &out) const;

  ///Returns a Newick representation of the tree's living members at a given
  ///time. See: https://en.wikipedia.org/wiki/Newick_format
  std::string printNewick(int n=0, int depth=0) const;

  ///Prints each species' SpeciesStats vector to the specified file
  void speciesSummaries(std::ofstream &out) const;
};

#endif