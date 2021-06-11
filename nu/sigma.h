#include <math.h>
/*
 *  xF differential cross-sections for proton Air collisions
 *
 *--------------------------------------------------------------------*/

double dsig_1(double E, double xF) { // scaling function
  if (xF>1.) {return 1e-12;}
  return pow(1.-xF,6.)/xF;
}

double dsig_2(double E, double xF) { // MRST-1998
  if (xF>1.) {return 1e-12;}
  return pow(1.-xF,6.)*(1.-2.8*sqrt(xF)+2.8*xF)/xF;
}
