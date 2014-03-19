#ifndef _temp_hpp_
#define _temp_hpp_

#include <vector>

//Define a singleton class for Temperatures
class Temperature {
 private:
  std::vector<double> temps;
  ///If this is TRUE then any calls to getTemp() return testTemp instead. This
  ///is FALSE by default.s
  bool test;
  ///testTemp is the value to return from getTemp() if test is TRUE
  double testTemp;
  Temperature();                               ///Prevent construction
  Temperature(const Temperature&);             ///Prevent copying
  Temperature& operator=(const Temperature&);  ///Prevent assignment
 public:
  ///Return the singleton instance
  static Temperature& getInstance() {
    static Temperature instance;
    return instance;
  }
  ///Load temperature data from the filename
  void init(const char filename[]);

  ///Get temperature at tMyrs performing interpolation if necessary. If test
  ///mode is on, then the actual temperature data is ignored and testTemp is
  ///returned.
  double getTemp(double tMyrs) const;

  ///Turn test mode on. If test mode is on, then calls to getTemp() return temp.
  void testOn(double temp);
  
  ///Turn test mode off.
  void testOff();
};

#endif