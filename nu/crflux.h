#include <math.h>
/*
 * cosmic ray flux(es)
 *
 *--------------------------------------------------------------------*/

extern double g; // spectral index

double phiN(double Ep) { // units?
  return pow(Ep,-g);
}

double phi_knee(double Ep) { // Thunman 1996
  if (Ep<5e6) return 1.7/pow(Ep,2.7);
  else return 174./pow(Ep,3.);
}

double phi_H3a(double Ep) { // Goncalves 2017
  double I0  = 1.15, g1  = 1.65, g2  = 2.4, Ek  = 1.2e6, eps = 3.0;
  return I0*pow(Ep,-1.-g1)*pow( 1. + pow(Ep/Ek,eps), (g1-g2)/eps );
}
