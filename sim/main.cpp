#include "salamander.hpp"
#include "mtbin.hpp"
#include <array>
using namespace std;

int main(){
  //One elevation bin of the mountain
  typedef array<Salamander, 1000> mtbin;

  //Profile of the mountain
  array<mtbin, 10000> mts;

  //Start off by killing everything
  for(auto &m: mts)
    m.killAll();

  //Find cell corresponding to optimal happy temperature
  //Make a salamander in that cell alive
  //Go forth

  //Salamanders are at their very happiest at 12.804279 degC
  //At sea level 65MYA the temp was 33.5618604122814 degC
  //Temperature drops at 9.8 degC/1000m (adiabatic lapse rate of dry air)
  //So (33.5618604122814-12.804279)/9.8=2.1181*1000m=2.11km


  
}
