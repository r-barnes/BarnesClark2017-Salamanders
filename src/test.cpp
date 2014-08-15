#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include "simulation.hpp"
#include "temp.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
using namespace std;

int main(){
  Temperature::getInstance().init("../data/temp_series_degreesC_0_65MYA_by_0.001MY.csv");

  cerr<<"Verifying temperatures"<<endl;
  for(double tMyrs=0;tMyrs<65;tMyrs+=0.5)
    cerr<<"Temperature at "<<tMyrs<<"\t=\t"<<Temperature::getInstance().getTemp(tMyrs)<<endl;

  cerr<<endl<<"Carrying Capacity"<<endl;
  cerr<<"Elevation   KKap(0)  KKap(64.9)"<<endl;
  for(int i=0;i<Simulation::numbins;i++){
    MtBin bin(2.8/Simulation::numbins*i, true);
    cerr<<setw(9)<<bin.height()<<setw(10)<<bin.kkap(0)<<setw(12)<<bin.kkap(64.9)<<endl;
  }

  MtBinUnitTest mtbintest;
  mtbintest.run();
}