#ifndef _simulation_hpp_
#define _simulation_hpp_

#include <vector>
#include <fstream>
#include "mtbin.hpp"
#include "phylo.hpp"

//This class will hold the parameters used to control a simulation. Running the
//simulation will result in the creation of a phylogeny and the setting of
//summary variables. Since the simulations are stored as instantiations of this
//class, it is simple to parallelize the program.
class Simulation {
 private:
   std::vector<MtBin> mts;

 public:
  Simulation(
    double mutation_probability0, 
    double temperature_drift_sd0, 
    double species_sim_thresh0, 
    double tempdeathfactor0, 
    double timestep0,
    bool vary_height0
  );

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
  //A value [0,1] indicating the probability of mutation (see salamander.hpp)
  double    mutation_probability;
  //A value [0, Inf] dictating the standard deviation of the random, normal
  //change in temperature tolerance that occurs between parents and offspring.
  double    temperature_drift_sd;
  //A value [0,1] indicating how similar to salamanders must be to be the same species
  double    species_sim_thresh;
  //Average optimal temperature of all the salamanders alive at the end of the run
  double    avg_otempdegC;
  //Number of living species at the end of the simulation
  int       nspecies;
  //Number of salamanders alive at the end of the simulation
  int       salive;
  //Phylogeny resulting from running the simulation
  Phylogeny phylos;
  //Dumps the phylogeny object to save space
  void dumpPhylogeny();
  //Empirical cumulative distribution (ECDF) of average branch lengths between extant taxa
  double    ecdf;
  //Length of a timestep in the simulation
  double    timestep;
  //Time at which the simulation ended
  double    endtime;
  //Average elevation of the salamanders at the end of the simulation
  double    avg_elevation;
  //Adjusts how salamander mortality relates to temperature. See
  //Salamander::pDie() for details.
  double    tempdeathfactor;
  //Determines whether the mountains erode over time
  bool      vary_height;
};

#endif
