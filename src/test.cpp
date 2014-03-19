#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include "simulation.hpp"
#include "data.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
using namespace std;

int main(){
  Temperature::getInstance().init("../data/temp_series_degreesC_0_65MYA_by_0.001MY.csv");

  cerr<<"Verifying temperatures"<<endl;
  for(double tMyrs=0;tMyrs<65;tMyrs+=0.5)
    cerr<<"Temperature at "<<tMyrs<<"\t=\t"<<Temperature::getInstance().getTemp(tMyrs)<<endl;

  MtBinUnitTest mtbintest;
  mtbintest.run();
}