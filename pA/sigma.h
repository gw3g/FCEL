#include <math.h>
/*
 *  y differential cross-sections for pA collisions
 *
 *--------------------------------------------------------------------*/

double dsig(double y, void *params) { // scaling function
  double n       = ((double *)params)[0]; // exponent
  double kappa   = ((double *)params)[1]; // = 2.*mT/sqrt(s)/(1-xi)
  return pow(1.-kappa*cosh(y),n);
}

