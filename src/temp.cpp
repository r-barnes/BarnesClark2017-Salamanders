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

  //Copy the last value of the array to ensure that we can interpolate right up
  //to the end of the time series
  temps.push_back(temps.back());
}

double Temperature::getTemp(double tMyrs) const {
  //Ensure that we've read temperature data
  assert(temps.size()!=0);

  //If we are in testing mode, return the temperature the user wanted to test
  //with
  if(test) return testTemp;
  
  double timeKyrs = tMyrs*1000;

  //Assumes that time is always positive
  unsigned int t0=(unsigned int)timeKyrs; //Time of the start of the 1kyr bin 

  //Ensure that we have points to interpolate with
  assert(t0 >= 0);

  //temps.size() is 1 larger than the position of the last element of the array.
  //temps.size()-1 is that last element. But we have copied the last element so
  //that there are two instances of it. This means that temps.size()-2 is the
  //last temperature period for which we can perform an interpolation.
  assert(t0 < temps.size()-2);

  double ta    = temps[t0];      //Temperature of this 1kyr bin
  double tb    = temps[t0+1];    //Temperature of the next 1kyr bin
  double tdiff = tb-ta;          //Temperature difference between the two bins
  return ta + tdiff*(timeKyrs-t0); //Perform the interpolation
}

void Temperature::testOn(double temp){
  test     = true;
  testTemp = temp;
  std::cerr<<"Temperature test mode on. Temp set to "<<testTemp<<"degC"<<std::endl;
}

void Temperature::testOff(){
  test = false;
  std::cerr<<"Temperature test mode off."<<std::endl;
}