#include "random.hpp"
#include "randutil.hpp"
#include <cassert>
#include <random>
#include <cstdint>

#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num()  0
	#define omp_get_num_threads() 1
	#define omp_get_max_threads()	1
#endif

our_random_engine& rand_engine(){
  static our_random_engine e[PRNG_THREAD_MAX];
  return e[omp_get_thread_num()];
}


//Be sure to read: http://www.pcg-random.org/posts/cpp-seeding-surprises.html
//and http://www.pcg-random.org/posts/cpps-random_device.html
void seed_rand(unsigned long seed){
  #pragma omp critical
  {
    if(seed==0)
      rand_engine().seed( randutils::auto_seed_128{}.base() );
    else
      rand_engine().seed( seed*omp_get_thread_num() );
  }
}


int uniform_rand_int(int from, int thru){
  static std::uniform_int_distribution<> d[PRNG_THREAD_MAX];
  using parm_t = std::uniform_int_distribution<>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}


double uniform_rand_real(double from, double thru){
  static std::uniform_real_distribution<> d[PRNG_THREAD_MAX];
  using parm_t = std::uniform_real_distribution<>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{from, thru} );
}


double normal_rand(double mean, double stddev){
  static std::normal_distribution<double> d[PRNG_THREAD_MAX];
  using parm_t = std::normal_distribution<double>::param_type;
  return d[omp_get_thread_num()]( rand_engine(), parm_t{mean, stddev} );
}