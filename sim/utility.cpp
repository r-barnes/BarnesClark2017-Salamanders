#include <cstdlib>
#include <cmath>


double ele_fxn(double elevation, double time) { //Give elevation as altitude in 1000m, and t in 1000 years ago (i.e. 0.001MYA)
  double elek=12.029753, elesigma=0.211410, elemu=0.245547, pi=3.141593, deltasd=2.881572e-09;
  double area;
  
  
  area = elek*(1/sqrt(2*pi*pow(elesigma+time*deltasd*1000,2)))*exp((-pow(elevation-elemu,2))/(2*pow(elesigma+time*deltasd*1000,2)));
  //Outputs area in units of km2
  area *= 100000;
  return area;
}
