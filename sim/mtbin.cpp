#include "mtbin.hpp"

MtBin::MtBin(){}

MtBin::MtBin(double binx) : MtBin(binx) {}

void MtBin::Mortaliate() {
  double mytemp=temp();

  for(auto &s: bin){
    if(s.pDie(mytemp)){}
  }
}

double MtBin::temp() const {
  
}
