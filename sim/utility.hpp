#ifndef _utility
#define _utility

#include <iostream>

template<class T>
void printbits(T a){
  T selector=1;
  for(unsigned int i=0;i<sizeof(T)*8;++i){
    std::cerr<<((a&selector)?'1':'0');
    selector=selector<<1;
  }
  std::cerr<<std::endl;
}

template<class T>
T countbits(T a){
  unsigned int c; // c accumulates the total bits set in combined
  for (c = 0; a; ++c)
    a &= a - 1; // clear the least significant bit set
  return c;
}

#endif
