#ifndef _prng_header
#define _prng_header

///Maximum number of threads this class should deal with
#define PRNG_THREAD_MAX 32

#include <random>

std::default_random_engine& rand_engine();

void seed_rand();

int uniform_rand(int from, int thru);

double uniform_rand(double from, double thru);

double normal_rand(double mean, double stddev);

#endif
