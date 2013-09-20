#ifndef _utility
#define _utility

///Give elevation as altitude in kilometers, and t in kiloyears ago (i.e. 0.001MYA)
double MountainElevation(double elevation, double time);

template<class T>
T countbits(T a){
  unsigned int c; // c accumulates the total bits set in combined
  for (c = 0; a; ++c)
    a &= a - 1; // clear the least significant bit set
  return c;
}

#endif
