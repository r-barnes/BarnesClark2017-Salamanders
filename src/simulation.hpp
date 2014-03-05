#ifndef _simulation_hpp_
#define _simulation_hpp_

#include <vector>
#include "mtbin.hpp"
#include "phylo.hpp"

//Create a run struct. This struct will store the mutation and species
//similarity thresholds used to run each of the simulations, along with the
//outputs of those simulations. This allows us to easily multi-thread the
//program.
class Simulation {
 private:
   std::vector<MtBin> mts;
 public:
  Simulation(double mutation_probability0, double temperature_drift_sd0, double species_sim_thresh0);

  /**
    This function runs a simulation with the specified mutation probability and
    species similarity threshold. It returns the resulting phylogeny. The function
    is thread-safe.

    @param[in]  mutation_probability  Mutation probability in TODO
    @param[in]  temperature_drift_sd  Std.Dev. of change in temperature tolerance
                                      that occurs between parents and children.
    @param[in]  species_sim_thresh    Threshold for two individuals to be considered
                                      members of the same species. Has value [0,1].
    @param[out] avg_otempdegC         Average optimal temperature of all surviving
                                      species at the end of the run


    @returns   A phylogeny object
  */
  void      runSimulation();
  int       alive() const;
  double    AvgOtempdegC() const;
  //A value [0,1] indicating the probability TODO
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
  //Empirical cumulative distribution (ECDF) of average branch lengths between extant taxa
  double    ecdf;
};

#endif