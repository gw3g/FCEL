#include <math.h>
/*
 * cosmic ray flux(es)
 *
 *--------------------------------------------------------------------*/

extern double g; // spectral index

double phiN(double Ep) { // units?
  return pow(Ep,-g);
}

double phi_knee(double Ep) {
  if (Ep<5e6) return 1.7/pow(Ep,2.7);
  else return 174./pow(Ep,3);
}

