//Define a singleton class for Temperatures

#include "temp.hpp"
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>

Temperature::Temperature(){
  test = false;
}

void Temperature::init(const char filename[]) {
  assert(temps.size()==0);

  //Read in temperatures
  std::ifstream fin(filename);
  double temp;
  while(fin>>temp)
    temps.push_back(temp);
  std::reverse(temps.begin(),temps.end());
}

double Temperature::getTemp(double tMyrs) const {
  assert(temps.size()!=0);

  if(test) return testTemp;
  
  double timeKyrs = tMyrs*1000;

  unsigned int t0=(int)timeKyrs; //Time of the start of the 1kyr bin 

  assert(t0 >= 0);
  assert(t0 < temps.size()-2);

  double ta    = temps[t0];      //Temperature of this 1kyr bin
  double tb    = temps[t0+1];    //Temperature of the next 1kyr bin
  double tdiff = tb-ta;          //Temperature difference between the two bins
  return ta + tdiff*(timeKyrs-t0);
}

void Temperature::testOn(double temp){
  test=true;
  testTemp=temp;
  std::cerr<<"Temperature test mode on. Temp set to "<<testTemp<<"degC"<<std::endl;
}

void Temperature::testOff(){
  test=false;
  std::cerr<<"Temperature test mode off."<<std::endl;
}