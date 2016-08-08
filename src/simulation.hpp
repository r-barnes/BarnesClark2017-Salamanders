#ifndef _simulation_hpp_
#define _simulation_hpp_

#include <vector>
#include <fstream>
#include "mtbin.hpp"
#include "phylo.hpp"
#include "params.hpp"

//This class will hold the parameters used to control a simulation. Running the
//simulation will result in the creation of a phylogeny and the setting of
//summary variables. Since the simulations are stored as instantiations of this
//class, it is simple to parallelize the program.
class Simulation {
 private:
   std::vector<MtBin> mts;

 public:
  Simulation(const Params& params);
  //All the simulation's parameters
  Params    params;
  //We use this many elevation bins to represent the mountain
  static const int numbins = 10;
  //Runs the simulations described by the following properties
  void      runSimulation();
  //Number of living salamanders
  int       alive() const;
  //Returns the average optimal temperature of the living salamanders
  double    AvgOtempdegC() const;
  //Returns the average elevation at which the living salamanders are found on
  //the mountain
  double    AvgElevation() const;
  //Average optimal temperature of all the salamanders alive at the end of the run
  double    avg_otempdegC = 0;
  //Number of living species at the end of the simulation
  int       nspecies = 0;
  //Number of salamanders alive at the end of the simulation
  int       salive = 0;
  //Phylogeny resulting from running the simulation
  Phylogeny phylos;
  //Dumps the phylogeny object to save space
  void dumpPhylogeny();
  //Empirical cumulative distribution (ECDF) of average branch lengths between extant taxa
  double    ecdf;
  //Time at which the simulation ended
  double    endtime = 0;
  //Average elevation of the salamanders at the end of the simulation
  double    avg_elevation = 0;
};

#endif
