#ifndef _simulation_hpp_
#define _simulation_hpp_

#include <vector>
#include <fstream>
#include "mtbin.hpp"
#include "phylo.hpp"

//This class will hold the parameters used to control a simulation. Running the
//simulation will cause certain summary variables to begin values. This allows
//us to easily multi-thread the program.
class Simulation {
 private:
   std::vector<MtBin> mts;

 public:
  Simulation(double mutation_probability0, double temperature_drift_sd0, double species_sim_thresh0, double tempdeathfactor0, double timestep0, bool vary_height0);

  //We use this many elevation bins to represent the mountain
  static const int numbins = 10;

  
  //This function runs a simulation with the specified mutation probability and
  //species similarity threshold. It returns the resulting phylogeny. The function
  //is thread-safe.
  //
  // @param[in]  mutation_probability  Mutation probability in per timestep
  // @param[in]  temperature_drift_sd  Std.Dev. of change in temperature 
  //                                   tolerance that occurs between parents and
  //                                   children.                        
  // @param[in]  species_sim_thresh    Threshold for two individuals to be
  //                                   considered members of the same species. 
  //                                   Has value [0,1].
  // @param[out] avg_otempdegC         Average optimal temperature of all
  //                                   surviving species at the end of the run
  //
  //
  // @returns   A phylogeny object
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
  bool vary_height;
};

#endif
