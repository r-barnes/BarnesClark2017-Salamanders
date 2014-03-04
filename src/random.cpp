#include "random.hpp"
#include <cassert>
#include <random>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num()  0
	#define omp_get_num_threads() 1
	#define omp_get_max_threads()	1
#endif

std::default_random_engine& rand_engine(){
  static std::default_random_engine e[PRNG_THREAD_MAX];
  return e[omp_get_thread_num()];
}

void seed_rand(){
  static std::random_device rd{};
  #pragma omp parallel
  {
    #pragma omp critical
    rand_engine().seed( rd() );
  }
}

int uniform_rand(int from, int thru){
  static std::uniform_int_distribution<> d[PRNG_THREAD_MAX];
  using parm_t = std::uniform_int_distribution<>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}

double uniform_rand(double from, double thru){
  static std::uniform_real_distribution<> d[PRNG_THREAD_MAX];
  using parm_t = std::uniform_real_distribution<>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}

double normal_rand(double mean, double stddev){
  static std::normal_distribution<double> d[PRNG_THREAD_MAX];
  using parm_t = std::normal_distribution<double>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{mean, stddev} );
}