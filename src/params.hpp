#ifndef _sal_params_
#define _sal_params_

#include <fstream>

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

  bool vary_height;
  bool vary_temp;
  bool run_once;

  double mortality_prob;
  double temp_drift;
  double simthresh;
  
  double timestep;
  int maxiter;
  int random_seed;

  double dispersal_prob;
  std::string dispersal_type;
  
  std::string temp_series_filename;
  int initial_altitude;

 public:
  Params(std::string filename);

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

  std::string getOutSummaryFilename      () const;
  std::string getOutPersistFilename      () const;
  std::string getOutPhylogenyFilename    () const;
  std::string getOutSpeciesStatsFilename () const;
  bool        pVaryHeight                () const;
  bool        pVaryTemp                  () const;
  bool        pRunOnce                   () const;
  double      getMutationProb            () const;
  double      getTempDrift               () const;
  double      getSimthresh               () const;
  double      getTimestep                () const;
  int         getMaxiter                 () const;
  double      getDispersalProb           () const;
  std::string getDispersalType           () const;
  std::string getTempSeriesFilename      () const;
  int         getInitialAltitude         () const;
  int         getRandomSeed              () const;
};

#endif