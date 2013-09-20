#include <cstdlib>
#include <cmath>

///Give elevation as altitude in kilometers, and t in kiloyears ago (i.e. 0.001MYA)
double MountainElevation(double elevation, double time) {
  double elek=12.029753;
  double elesigma=0.211410;
  double elemu=0.245547;
  double pi=3.141593;
  double deltasd=2.881572e-09;

  double area = elek
                * ( 1/sqrt( 2*pi*pow(elesigma+time*deltasd*1000, 2) ))
                * exp(
                        ( -pow(elevation-elemu, 2) )
                      / (2 * pow(elesigma+time*deltasd*1000, 2))
                  );
  area *= 100000; //Convert area to units of km2
  return area;
}
