#include <math.h>
/*
 *  xF differential cross-sections for proton Air collisions
 *
 *--------------------------------------------------------------------*/

double dsig_1(double E, double xF) { // scaling function
  if (xF>1.) {return 1e-16;}
  return pow(1.-xF,3.)/xF;
}

double dsig_2(double E, double xF) { // MRST-1998
  if (xF>1.) {return 1e-12;}
  return pow(1.-xF,6.)*(1.-2.8*sqrt(xF)+2.8*xF)/xF;
}

double dsig_3(double E, double xF) { // Martin 2003 (A.1)
  double be = .05 - .016*log(E/10000.), n, A;
  if (xF>1.) {return 1e-12;}
  else if (E<1e8) {
    n  = 7.6 ;//+ .025*log(E/10000.);
    A  = 140. + pow( 11.*log( E/100. ), 1.65);
  }
  else if (E>1e8) {
    n  = 7.6 ;//+ .012*log(E/10000.);
    A  = 4100. + 245.*log( E/1e8 );
  }
  return A*pow(xF,0.-1.)*pow(1.-pow(xF,1.2),n);
}
