#include <math.h>
/*
 *  xF differential cross-sections for proton Air collisions
 *
 *--------------------------------------------------------------------*/

double dsig_1(double E, double xF) { // scaling function
  if (xF>1.) {return 1e-16;}
  return pow(1.-xF,6.)/xF;
}

double dsig_2(double E, double xF) { // MRST-1998
  if (xF>1.) {return 1e-12;}
  return pow(1.-xF,6.)*(1.-2.8*sqrt(xF)+2.8*xF)/xF;
}

double dsig_3(double E, double xF) { // Martin 2003
  double A  = 140. + pow( 11.*log( E/100. ), 1.65);
  double n  = 7.6 + .025*log(E/10000.);
  double be = .05 - .016*log(E/10000.);
  return A*pow(xF,be-1.)*pow(1.-pow(xF,1.2),n);
}
