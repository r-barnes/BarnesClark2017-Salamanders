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

  bool vary_height; ///If this is set to true, the mountains erode over time.
  bool vary_temp;
  bool run_once;    //TODO: Cut
  bool debug_val;   ///Shows real-time stats about simulation state

  ///Mutation probability per timestep - used by salamander::mutate() to
  ///determine the probability of mutation in the genome. Note - this changes
  ///the genome that determines relatedness, speciation, and ability to breed.
  ///It does not directly alter optimum temperature. Should be [0,1].
  double    mutation_probability;
  ///Drift rate for temperature optimum per timestep - used by
  ///salamander::breed() to determine the change in optimum temperature between
  ///children and parents. Note - this changes the temperature optimum, but does
  ///not directly influence relatedness, speciation, and ability to breed.
  ///Should be [0,Inf).
  double    temperature_drift_sd;
  //A value [0,1] indicating how similar to salamanders must be to be the same species
  int       species_sim_thresh;


  ///This variable adjusts how the the salamander's probability of death is
  ///affected by temperature. For further details, please look at
  ///Salamander::pDie()
  double logit_temp_weight;
  double logit_offset;     
  double logit_ca_weight;  
  double logit_ha_weight;  

  //Length of a timestep in the simulation
  double timestep_val;
  int maxiter_val;
  int random_seed;

  double dispersal_prob;
  int    dispersal_type;
  
  std::string temp_series_filename;
  int initial_altitude;

  //Maximum number of new offspring per bin per unit time
  int max_offspring_per_bin_per_dt;

  int max_tries_to_breed;

  int initial_pop_size;

  //We use this many elevation bins to represent the mountain
  int numbins_val;

  double to_lowlands_prob;
  double from_lowlands_prob;

 public:
  Params();
  void load(std::string filename);

  // void setVaryHeight         (std::ifstream &fparam);
  // void setVaryTemp           (std::ifstream &fparam);
  // void setRunOnce            (std::ifstream &fparam);
  // void setMprob              (std::ifstream &fparam);
  // void setTdrift             (std::ifstream &fparam);
  // void setSimthresh          (std::ifstream &fparam);
  // void setTimestep           (std::ifstream &fparam);
  // void setMaxiter            (std::ifstream &fparam);
  // void setDispersalProb      (std::ifstream &fparam);
  // void setDispersalType      (std::ifstream &fparam);
  // void setTempSeriesFilename (std::ifstream &fparam);
  // void setInitialAltitude    (std::ifstream &fparam);

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