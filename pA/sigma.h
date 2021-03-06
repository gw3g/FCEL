#include <math.h>
/*
 *  y differential cross-sections for pA collisions
 *
 *--------------------------------------------------------------------*/

double dsig(double y, void *params) { // scaling function
  double n       = ((double *)params)[0]; // exponent
  double kappa   = ((double *)params)[1]; // = 2.*mT/sqrt(s)/(1-xi) ?
  double res =  1.-kappa*cosh(fmin(y,acosh(1./kappa)));
  //printf("y=%g, kappa=%g, res= %g\n",y,kappa,res);
  //if (res<1e-3) return 1e-5;
  return pow(res,n);
}

