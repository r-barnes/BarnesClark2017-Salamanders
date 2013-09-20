#include <cstdlib>
#include <cmath>

///Give elevation as altitude in kilometers, and t in kiloyears ago (i.e. 0.001MYA)
double MountainElevation(double elevation, double time) {
  const double elek=12.029753;
  const double elesigma=0.211410;
  const double elemu=0.245547;
  const double pi=3.141593;
  const double deltasd=2.881572e-09;

  double area = elek
                * ( 1/sqrt( 2*pi*pow(elesigma+time*deltasd*1000, 2) ))
                * exp(
                        ( -pow(elevation-elemu, 2) )
                      / (2 * pow(elesigma+time*deltasd*1000, 2))
                  );
  area *= 100000; //Convert area to units of km2
  return area;
}
