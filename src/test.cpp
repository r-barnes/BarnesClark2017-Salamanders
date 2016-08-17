#include "salamander.hpp"
#include "mtbin.hpp"
#include "phylo.hpp"
#include "simulation.hpp"
#include "temp.hpp"
#include "random.hpp"
#include <array>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <string>
#include <bitset>
using namespace std;

int main(int argc, char **argv){
  Temperature.init("../data/temp_series_degreesC_0_65MYA_by_0.001MY.csv");

  if(argc!=2){
    cerr<<"Syntax: "<<argv[0]<<" <Parameters File>"<<endl;
    return -1;
  }

  TheParams.load(argv[1]);

  seed_rand(0);

  cout<<"Verifying temperatures"<<endl;
  for(double tMyrs=0;tMyrs<65.0001;tMyrs+=0.5)
    cout<<"Temperature at "<<tMyrs<<"\t=\t"<<Temperature.getTemp(tMyrs)<<endl;

  //Test quality of random number generator using
  // ./test.exe  | grep -i Rand | sed 's/.*://' | tr " " "\n" | sort | 
  //                  uniq -c | awk '{print $0" "($1/100000*100)}'

  cout<<"100000 random integers in the range [0,9] from Mersenne: ";
  for(int i=0;i<100000;i++)
    cout<<uniform_rand_int(0,9)<<" ";
  cout<<endl;

  cout<<"Maximum height over time: \n";
  for(double tMyrs;tMyrs<65.001;tMyrs+=0.5)
    cout<<"Maximum elevation at "
        <<setw(4)<<tMyrs<<"\t=\t"<<MtBin::heightMaxKm(tMyrs)<<endl;

  cout<<"10000 random uint64 bit fields: ";
  for(int i=0;i<10000;i++)
    cout<<std::bitset<64>(uniform_bits<uint64_t>())<<" ";
  cout<<endl;

  MtBinUnitTest mtbintest;
  mtbintest.run();
}