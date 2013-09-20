#include "salamander.hpp"
#include <array>
using namespace std;

int main(){
  array< array<Salamander, 1000>, 10000> mts;

  //Start off by killing everything
  for(unsigned int x=0;x<mts.size();++x)
    for(auto &s: mts[x])
      s.dead=true;

  //Find cell corresponding to optimal happy temperature
  //Make a salamander in that cell alive
  //Go forth
}
