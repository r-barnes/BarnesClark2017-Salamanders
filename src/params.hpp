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

  //Determines whether the mountains erode over time
  bool vary_height;
  bool vary_temp;
  bool run_once;

  //A value [0,1] indicating the probability of mutation (see salamander.hpp)
  double    mutation_probability;
  //A value [0, Inf] dictating the standard deviation of the random, normal
  //change in temperature tolerance that occurs between parents and offspring.
  double    temperature_drift_sd;
  //A value [0,1] indicating how similar to salamanders must be to be the same species
  double    species_sim_thresh;


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
  double      mutationProb            () const;
  double      tempDrift               () const;
  double      speciesSimthresh        () const;
  double      timestep                () const;
  int         maxiter                 () const;
  double      dispersalProb           () const;
  int         dispersalType           () const;
  std::string tempSeriesFilename      () const;
  int         initialAltitude         () const;
  int         randomSeed              () const;
  double      tempDeathFactor         () const;
  double      logitTempWeight         () const;
  double      logitOffset             () const;
  double      logitCAweight           () const;
  double      logitHAweight           () const;
};


class TheParams {
 private:
  TheParams(){}
 public: 
  static Params& get(){
    static Params instance;     // This object is not created until the first
    return instance;                  // time getInstance() is called.
  }
};

#endif