//Define a singleton class for Temperatures

#include "temp.hpp"
#include <fstream>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <string>

Temperature::Temperature(){
  test = false;
}


//Use the RAII pattern. Load a file and read in the temperature data.
void Temperature::init(const std::string filename) {
  assert(temps.size()==0);

  //Read in temperatures
  std::ifstream fin(filename);
  if(!fin.good()){
    std::cerr<<"Could not open temperature file '"<<filename<<"'!"<<std::endl;
    assert(fin.good());
  }

  double temp;                              //Temporary variable for reading
  while(fin>>temp)                          //Read data, if available
    temps.push_back(temp);                  //Put in back of series

  //The data file we used needs to be reversed (if it is read in the above
  //manner) in order to be in chronological order.
  std::reverse(temps.begin(),temps.end());

  //Copy the last value of the array to ensure that we can interpolate right up
  //to the end of the time series
  temps.push_back(temps.back());
}


//Get the temperature at tMyrs, interpolating if necessary.
double Temperature::getTemp(double tMyrs) const {
  //Ensure that we've read temperature data
  assert(temps.size()!=0);

  //If we are in testing mode, return the temperature the user wanted to test
  //with.
  if(test) return testTemp;

  double timeKyrs = tMyrs*1000;

  //Assumes that time is always positive
  unsigned int t0 = (unsigned int)timeKyrs; //Time of the start of the 1kyr bin

  //Ensure that we have points to interpolate with
  assert(t0 >= 0);

  //temps.size() is 1 larger than the position of the last element of the array.
  //temps.size()-1 is that last element. But we have copied the last element so
  //that there are two instances of it. This means that temps.size()-2 is the
  //last temperature period for which we can perform an interpolation.
  assert(t0 < temps.size()-2);

  //The following performs the interpolation
  double ta    = temps[t0];        //Temperature of this 1kyr bin
  double tb    = temps[t0+1];      //Temperature of the next 1kyr bin
  double tdiff = tb-ta;            //Temperature difference between the two bins
  return ta + tdiff*(timeKyrs-t0); //Perform the interpolation
}


//Turn testing mode on
void Temperature::testOn(double temp){
  test     = true;
  testTemp = temp;
  std::cerr<<"Temperature test mode on. Temp set to "<<testTemp<<"degC"<<std::endl;
}


//Turn testing mode off
void Temperature::testOff(){
  test = false;
  std::cerr<<"Temperature test mode off."<<std::endl;
}