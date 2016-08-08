#include <fstream>
#include <iostream>
#include <stdexcept>
#include "params.hpp"

Params::Params(std::string filename){
  std::ifstream fparam(filename);

  out_summary       = Input_Filename(fparam,"SummaryStatsFilename");
  out_persist       = Input_Filename(fparam,"PersistenceGraphFilename");
  out_phylogeny     = Input_Filename(fparam,"PhylogenyFilename");
  out_species_stats = Input_Filename(fparam,"SpeciesStatsFilename");

  vary_height = Input_YesNo(fparam,"VaryHeight");
  vary_temp   = Input_YesNo(fparam,"VaryTemp");
  run_once    = Input_YesNo(fparam,"RunOnce");

  mortality_prob = Input_Double(fparam, "MutationProb");
  temp_drift     = Input_Double(fparam, "TemperatureDrift");
  simthresh      = Input_Double(fparam, "SpeciesSimilarity");
  timestep       = Input_Double(fparam, "timestep");

  dispersal_prob = Input_Double(fparam, "DispersalProb");
  dispersal_type = Input_Filename(fparam,"DispersalType");

  maxiter     = Input_Integer(fparam,"maxiter");
  random_seed = Input_Integer(fparam,"PRNGseed");

  std::string temp_series_filename = Input_Filename(fparam,"TempSeries");

  initial_altitude = Input_Integer(fparam,"InitialAltitude");

  {
    std::string test_if_file_is_empty;
    if(fparam>>test_if_file_is_empty){
      std::cerr<<"File contained unexpected parameters at the end!"<<std::endl;
      throw std::runtime_error("Unexpected parameters in parameter file!");
    }
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

std::string Params::getOutSummaryFilename      () const {return out_summary;          }
std::string Params::getOutPersistFilename      () const {return out_persist;          }
std::string Params::getOutPhylogenyFilename    () const {return out_phylogeny;        }
std::string Params::getOutSpeciesStatsFilename () const {return out_species_stats;    }
bool        Params::pVaryHeight                () const {return vary_height;          }
bool        Params::pVaryTemp                  () const {return vary_temp;            }
bool        Params::pRunOnce                   () const {return run_once;             }
double      Params::getMutationProb            () const {return mortality_prob;       }
double      Params::getTempDrift               () const {return temp_drift;           }
double      Params::getSimthresh               () const {return simthresh;            }
double      Params::getTimestep                () const {return timestep;             }
int         Params::getMaxiter                 () const {return maxiter;              }
double      Params::getDispersalProb           () const {return dispersal_prob;       }
std::string Params::getDispersalType           () const {return dispersal_type;       }
std::string Params::getTempSeriesFilename      () const {return temp_series_filename; }
int         Params::getInitialAltitude         () const {return initial_altitude;     }
int         Params::getRandomSeed              () const {return random_seed;          }