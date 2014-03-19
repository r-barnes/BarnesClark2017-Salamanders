#ifndef _temp_hpp_
#define _temp_hpp_

#include <vector>

//Define a singleton class for Temperatures
class Temperature {
 private:
  std::vector<double> temps;
  bool test;
  double testTemp;
  Temperature();                               //Prevent construction
  Temperature(const Temperature&);             //Prevent copying
  Temperature& operator=(const Temperature&);  //Prevent assignment
 public:
  static Temperature& getInstance() {
    static Temperature instance;
    return instance;
  }
  void init(const char filename[]);
  double getTemp(double tMyrs) const;
  void testOn(double temp);
  void testOff();
};

#endif