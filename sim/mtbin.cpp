#include "mtbin.hpp"
#include "data.hpp"

MtBin::MtBin(){
  startofdead=0;
}

MtBin::MtBin(double binx) : MtBin(binx) {}

void MtBin::mortaliate() {
  double mytemp=temp();

  for(auto &s: bin){
    if(s.pDie(mytemp)){
      s.dead=true;
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
