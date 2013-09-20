#include "mtbin.hpp"
#include "data.hpp"
#include <algorithm>

MtBin::MtBin(){
  startofdead=0;
}

MtBin::MtBin(double binx) : MtBin(binx) {}

void MtBin::mortaliate(double t) {
  ///If there are no living salamanders, then don't do anything
  if(startofdead==0) return;

  double mytemp=temp(t);

  for(unsigned int s=0;s<startofdead;++s){
    if(bin[s].pDie(mytemp)){
      bin[s].dead=true;
      std::swap(bin[s],bin[startofdead-1]);
      --startofdead;
    }
  }
}

double MtBin::temp(double t) const {
  
}

double MtBin::height(double t) const {

}

void MtBin::killAll() {
  for(auto &s: bin)
    s.dead=true;
  startofdead=0;
}

void MtBin::addSalamander(const Salamander &s) {
  if(startofdead==bin.size()) return;
  bin[startofdead]=s;
  ++startofdead;
}
