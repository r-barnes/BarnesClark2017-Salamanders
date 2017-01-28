#ifndef _sal_params_
#define _sal_params_

#include <fstream>

const int DISPERSAL_BETTER      = 1;
const int DISPERSAL_MAYBE_WORSE = 2;
const int DISPERSAL_GLOBAL      = 3;

class Params {
 private:
  void        Input_CheckParamName(std::ifstream &fparam, const std::string &param_name) const;
  std::string Input_Filename      (std::ifstream &fparam, const std::string &param_name) const;
  bool        Input_YesNo         (std::ifstream &fparam, const std::string &param_name) const;
  double      Input_Double        (std::ifstream &fparam, const std::string &param_name) const;
  int         Input_Integer       (std::ifstream &fparam, const std::string &param_name) const;

  std::string out_summary;
  std::string out_persist;
  std::string out_phylogeny;
  std::string out_species_stats;

  ///If this is set to true, the mountains erode over time.
  bool vary_height; 

  ///If this is set to true, the temperature changes over time.
  bool vary_temp;

  //DEPRECATED
  bool run_once;    

  ///Shows real-time stats about simulation state. Forces the simulation to
  ///perform only one realization using only one thread.
  bool debug_val;   

  ///Mutation probability per timestep - used by salamander::mutate() to
  ///determine the probability of mutation in the genome. Note - this changes
  ///the genome that determines relatedness, speciation, and ability to breed.
  ///It does not directly alter optimum temperature. Should be [0,1].
  double mutation_probability;

  ///Drift rate for temperature optimum per timestep - used by
  ///salamander::breed() to determine the change in optimum temperature between
  ///children and parents. Note - this changes the temperature optimum, but does
  ///not directly influence relatedness, speciation, and ability to breed.
  ///Should be [0,Inf).
  double temperature_drift_sd;

  ///The number of bits which must be the same for two salamanders to be
  ///considered part of the same species. Can range from (-Inf,Inf). Values less
  ///than zero assure relatedness while values greater than the number of bits
  ///in Salamander::genetype assure unrelatedness.
  int species_sim_thresh;

  ///Affects the salamander's probability of death. See Salamander::pDie()
  double logit_offset;    
  ///This variable adjusts how the the salamander's probability of death is
  ///affected by temperature. For further details, please look at
  ///Salamander::pDie()
  double logit_temp_weight;

  ///Affects how conspecific abundance affects the salamander's probability of
  ///death. See Salamander::pDie()
  double logit_ca_weight;  

  ///Affects how heterospecific abundance affects the salamander's probability
  ///of death. See Salamander::pDie()
  double logit_ha_weight;  

  ///Duration of a timestep in the simulation
  double timestep_val;

  ///How many realizations of the simulation to run
  int maxiter_val;

  ///Random seed used to initialization the simulation. Using 0 causes the PRNGs
  ///to be seeded with entropy.
  int random_seed;

  ///Probability of a salamander deciding to move to another bin, if another bin
  ///exists (the top and bottom of the mountain are special cases).
  double dispersal_prob;

  ///What kind of movement decision metric a salamander uses. Options are
  ///DISPERSAL_BETTER (only move to bins with more favourable climates),
  ///MAYBE_WORSE (move to one of the neighbouring bins regardless of its
  ///temperature), or GLOBAL (move to any bin on the mountain).
  int dispersal_type;
  
  ///File to read the temperature time series from
  std::string temp_series_filename;

  ///The bin the salamanders are initially in
  int initial_altitude;

  ///Maximum number of new offspring per bin per unit time
  int max_offspring_per_bin_per_dt;

  ///Random pairings of salamanders in each bin are chosen for breeding without
  ///regard as to whether they can breed. The pairing is rejected if they are
  ///not of the same species, and two different salamanders are randomly chosen.
  ///This parameter prevents infinite loops by stopping the process if too many
  ///tries are made.
  int max_tries_to_breed;

  ///How many salamanders are placed in the initial bin initially
  int initial_pop_size;

  ///We use this many elevation bins to represent the mountain
  int numbins_val;

  ///Probability per salamander per timestep of a salamander migrating from the
  ///bottom bin out into the surrounding lowlands. Range is [0,1]. Negative
  ///values turn this effect off.
  double to_lowlands_prob;

  ///Probability per salamander per timestep of a salamander migrating from the
  ///bottom bin out into the surrounding lowlands. Range is [0,1]. Negative
  ///values turn this effect off.
  double from_lowlands_prob;

 public:
  Params();
  void load(std::string filename);

  //A large number of access methods which return the above variables
  std::string outSummaryFilename      () const;
  std::string outPersistFilename      () const;
  std::string outPhylogenyFilename    () const;
  std::string outSpeciesStatsFilename () const;
  bool        pVaryHeight             () const;
  bool        pVaryTemp               () const;
  bool        pRunOnce                () const;
  int         numBins                 () const;
  double      mutationProb            () const;
  double      tempDrift               () const;
  int         speciesSimthresh        () const;
  double      timestep                () const;
  int         maxiter                 () const;
  double      dispersalProb           () const;
  int         dispersalType           () const;
  std::string tempSeriesFilename      () const;
  int         initialAltitude         () const;
  int         initialPopSize          () const;
  int         randomSeed              () const;
  double      tempDeathFactor         () const;
  double      logitTempWeight         () const;
  double      logitOffset             () const;
  double      logitCAweight           () const;
  double      logitHAweight           () const;
  int         maxOffspringPerBinPerDt () const;
  int         maxTriesToBreed         () const;
  double      toLowlandsProb          () const;
  double      fromLowlandsProb        () const;
  bool        debug                   () const;
};

extern Params TheParams;

#endif