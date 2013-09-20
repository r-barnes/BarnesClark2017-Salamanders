#ifndef _utility
#define _utility

double ele_fxn(double elevation, double time); //Give elevation as altitude in 1000m, and t in 1000 years ago (i.e. 0.001MYA)

template<class T>
T countbits(T a){
  unsigned int c; // c accumulates the total bits set in combined
  for (c = 0; a; ++c)
    a &= a - 1; // clear the least significant bit set
  return c;
}

#endif
