#ifndef _temp_hpp_
#define _temp_hpp_

#include <vector>
#include <string>

//Define a singleton class for Temperatures which loads them, stores them, and
//interpolates between them as necessary.
class TemperatureClass {
 private:
  //The temperature series
  std::vector<double> temps;

  ///If this is TRUE then any calls to getTemp() return testTemp instead. This
  ///is FALSE by default.s
  bool test;

  ///testTemp is the value to return from getTemp() if test is TRUE
  double testTemp;

  //The following ensure that this class is a singleton
  TemperatureClass(const TemperatureClass&);             ///Prevent copying
  TemperatureClass& operator=(const TemperatureClass&);  ///Prevent assignment
 public:
  TemperatureClass();

  ///Load temperature data from the filename
  void init(const std::string filename);

  ///Get temperature at tMyrs performing interpolation if necessary. If test
  ///mode is on, then the actual temperature data is ignored and testTemp is
  ///returned.
  double getTemp(double tMyrs) const;

  ///Turn test mode on. If test mode is on, then calls to getTemp() return temp.
  void testOn(double temp);
  
  ///Turn test mode off.
  void testOff();
};

extern TemperatureClass Temperature;

#endif