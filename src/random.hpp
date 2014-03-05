#ifndef _prng_header
#define _prng_header

///Maximum number of threads this class should deal with
#define PRNG_THREAD_MAX 32

#include <random>

//Returns a PRNG engine specific to the calling thread
std::default_random_engine& rand_engine();

//Seeds the PRNG engines using entropy from the computer's random device
void seed_rand();

//Returns an integer value on the closed interval [from,thru]
//Thread-safe
int uniform_rand_int(int from, int thru);

//Returns an floating-point value on the interval [from,thru)
//Thread-safe
double uniform_rand_real(double from, double thru);

//Returns a Gaussian-distributed value with specified mean and standard
//deviation. Thread-safe
double normal_rand(double mean, double stddev);

#endif
