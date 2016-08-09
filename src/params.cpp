#include <fstream>
#include <iostream>
#include <stdexcept>
#include "params.hpp"

Params::Params(){}

void Params::load(std::string filename){
  std::ifstream fparam(filename);

  out_summary       = Input_Filename(fparam,"SummaryStatsFilename");
  out_persist       = Input_Filename(fparam,"PersistenceGraphFilename");
  out_phylogeny     = Input_Filename(fparam,"PhylogenyFilename");
  out_species_stats = Input_Filename(fparam,"SpeciesStatsFilename");

  vary_height = Input_YesNo(fparam,"VaryHeight");
  vary_temp   = Input_YesNo(fparam,"VaryTemp");
  run_once    = Input_YesNo(fparam,"RunOnce");

  numbins_val = Input_Integer(fparam,"NumBins");

  mutation_probability = Input_Double(fparam, "MutationProb");
  temperature_drift_sd = Input_Double(fparam, "TemperatureDrift");
  species_sim_thresh   = Input_Double(fparam, "SpeciesSimilarity");
  timestep_val         = Input_Double(fparam, "timestep");

  dispersal_prob = Input_Double(fparam, "DispersalProb");

  {
    std::string temp = Input_Filename(fparam,"DispersalType");
    if(temp=="Global")
      dispersal_type = DISPERSAL_GLOBAL;
    else if(temp=="MaybeWorse")
      dispersal_type = DISPERSAL_MAYBE_WORSE;
    else if(temp=="Better")
      dispersal_type = DISPERSAL_BETTER;
    else {
      std::cerr<<"Unrecognised dispersal type! Expected: Global, MaybeWorse, Better"<<std::endl;
      throw std::runtime_error("Unrecognised dispersal type! Expected: Global, MaybeWorse, Better");
    }
  }

  maxiter_val = Input_Integer(fparam,"maxiter");
  random_seed = Input_Integer(fparam,"PRNGseed");

  temp_series_filename = Input_Filename(fparam,"TempSeries");

  initial_altitude = Input_Integer(fparam,"InitialAltitude");

  logit_temp_weight = Input_Double(fparam,"LogitTempWeight");
  logit_offset      = Input_Double(fparam,"LogitOffset");
  logit_ca_weight   = Input_Double(fparam,"LogitCAweight");
  logit_ha_weight   = Input_Double(fparam,"LogitHAweight");  

  max_offspring_per_bin_per_dt = Input_Integer(fparam,"MaxOffspringPerBinPerDt");
  max_tries_to_breed           = Input_Integer(fparam,"MaxTriesToBreed");

  {
    //std::string test_if_file_is_empty;
    //if(fparam>>test_if_file_is_empty){
    //  std::cerr<<fparam<<std::endl;
    //  std::cerr<<"File contained unexpected parameters at the end!"<<std::endl;
    //  throw std::runtime_error("Unexpected parameters in parameter file!");
    //}
  }
}




void Params::Input_CheckParamName(std::ifstream &fparam, const std::string &param_name) const {
  std::string in_param_name;
  fparam>>in_param_name;
  if(in_param_name!=param_name){
    std::cerr<<"Unexpected parameter name! Found '"<<in_param_name<<"' expected '"<<param_name<<"'"<<std::endl;
    throw std::runtime_error("Unexpected parameter name!");
  }
}

std::string Params::Input_Filename(std::ifstream &fparam, const std::string &param_name) const {
  Input_CheckParamName(fparam, param_name);
  std::string temp;
  fparam>>temp;
  return temp;
}

bool Params::Input_YesNo(std::ifstream &fparam, const std::string &param_name) const {
  Input_CheckParamName(fparam, param_name);
  std::string temp;
  fparam>>temp;
  if(! (temp=="YES" || temp=="NO")){
    std::cerr<<"Paremeter '"<<param_name<<"' must be YES or NO!"<<std::endl;
    throw std::runtime_error("Bad parameter value!");
  }
  return (temp=="YES");
}

double Params::Input_Double(std::ifstream &fparam, const std::string &param_name) const {
  Input_CheckParamName(fparam, param_name);
  double temp;
  fparam>>temp;
  return temp;
}

int Params::Input_Integer(std::ifstream &fparam, const std::string &param_name) const {
  Input_CheckParamName(fparam, param_name);
  int temp;
  fparam>>temp;
  return temp;
}

std::string Params::outSummaryFilename      () const {return out_summary;                  }
std::string Params::outPersistFilename      () const {return out_persist;                  }
std::string Params::outPhylogenyFilename    () const {return out_phylogeny;                }
std::string Params::outSpeciesStatsFilename () const {return out_species_stats;            }
bool        Params::pVaryHeight             () const {return vary_height;                  }
bool        Params::pVaryTemp               () const {return vary_temp;                    }
bool        Params::pRunOnce                () const {return run_once;                     }
double      Params::mutationProb            () const {return mutation_probability;         }
double      Params::tempDrift               () const {return temperature_drift_sd;         }
double      Params::speciesSimthresh        () const {return species_sim_thresh;           }
double      Params::timestep                () const {return timestep_val;                 }
int         Params::maxiter                 () const {return maxiter_val;                  }
double      Params::dispersalProb           () const {return dispersal_prob;               }
int         Params::dispersalType           () const {return dispersal_type;               }
std::string Params::tempSeriesFilename      () const {return temp_series_filename;         }
int         Params::initialAltitude         () const {return initial_altitude;             }
int         Params::randomSeed              () const {return random_seed;                  }
double      Params::logitTempWeight         () const {return logit_temp_weight;            }
double      Params::logitOffset             () const {return logit_offset;                 }
double      Params::logitCAweight           () const {return logit_ca_weight;              }
double      Params::logitHAweight           () const {return logit_ha_weight;              }
int         Params::maxOffspringPerBinPerDt () const {return max_offspring_per_bin_per_dt; }
int         Params::maxTriesToBreed         () const {return max_tries_to_breed;           }
int         Params::numBins                 () const {return numbins_val;                  }