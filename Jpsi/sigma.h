#include <math.h>
/*
 *  y differential cross-sections for pA collisions
 *
 *--------------------------------------------------------------------*/

double dsig(double xF, void *params) {     // scaling function
  double n       = ((double *)params)[0];  // exponent
  double kappa   = ((double *)params)[1];  // M2/s

  double xp = sqrt(xF*xF+4.*kappa);
  //printf("TEST - 3e\n"); 

  if (xp>1.) {return 1e-16;}
  else {
    //return pow(1.-xF,n)*(1.+eps_g*sqrt(xF)+gam*xF);
    return pow((1.-xp),n);
    //return pow((1.-xF/2.),n);
  }
}

