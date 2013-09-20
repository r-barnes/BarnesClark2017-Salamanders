#include "salamander.hpp"
#include "mtbin.hpp"
#include <array>
#include <deque>
using namespace std;

//Profile of the mountain
//We allocate this in the global namespace to prevent a stack overflow
static array< MtBin, 1000> mts;

typedef std::deque<Phylo> phylolist;

int main(){
  phylolist phylos;

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

  for(t=0;t<65001;t+=0.5){
    for(unsigned int m=0;m<mts.size();++m){
      mts[m].mortaliate();
      mts[m].breed();
      if(m>0) mts[m].diffuse(mts[m-1]);
      if(m<mts.size()-1) mts[m].diffuse(mts[m+1]);
    }
  }
}
