#include <math.h>
/*
 *  y differential cross-sections for pA collisions
 *
 *--------------------------------------------------------------------*/

double dsig(double y, void *params) {     // scaling function
  double n       = ((double *)params)[0]; // exponent
  double kappa   = ((double *)params)[1]; // = 2.*mT/sqrt(s)/(1-xi) ?
  double eps_g   = ((double *)params)[2];

  double xF = kappa*cosh(y);
  double gam = -eps_g;

  if (xF>1.) {return 1e-5;}
  else {
    //return pow(1.-xF,n)*(1.+eps_g*sqrt(xF)+gam*xF);
    //return pow((1.-xF)*(1.-sqrt(xF)),n);
    return pow((1.-xF/2.),n);
  }
}

