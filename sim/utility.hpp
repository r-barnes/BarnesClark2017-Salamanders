#ifndef _utility
#define _utility

template<class T>
T countbits(T a){
  unsigned int c; // c accumulates the total bits set in combined
  for (c = 0; a; ++c)
    a &= a - 1; // clear the least significant bit set
  return c;
}

#endif
