#include <math.h>
/*
 *  xF differential cross-sections for proton Air collisions
 *
 *--------------------------------------------------------------------*/

double dsig_1(double E, double xF) { // scaling function
  return pow(1.-xF,7.)/xF;
}

