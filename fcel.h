#include <math.h>

/*--------------------------------------------------------------------*/
// QCD
double      Nc = 3.;

#define Cf  (Nc*Nc-1)/(2.*Nc)
#define Tf  (1./2.)
#define df  Nc

#define Ca  Nc
#define Ta  Nc
#define da  (Nc*Nc-1.)

// keep it dimensionless, but sometimes need:
double mp     = .93827208816          ; // GeV
double hbarc  = .19732698041522198110 ; // GeV.fm

/*--------------------------------------------------------------------*/

double dilog(double x) { double res, z2 = M_PI*M_PI/6.;
  if ( x<2. && 3*x>1 ) { // intermediate x
    double a1 = M_LN2*(1.-x),
           a2=a1*a1, a4=a2*a2, a6=a2*a4, a8=a4*a4, a10=a4*a6, a12=a6*a6,
           a3=a1*a2, a5=a1*a4, a7=a1*a6, a9=a1*a8, a11=a1*a10;
    res  = - 1.;
    res += ( 1.2158542*a1  + .24439311*a2  + .08293384*a3
           + .03486953*a4  + .01649813*a5  + .00841503*a6
           + .00452197*a7  + .00252503*a8  + .00145210*a9
           + .00085478*a10 + .00051279*a11 + .00031249*a12 );
  } else {
    double b1 = -x;
    if ( x>2.) { // large-x approx.
      res = log(x); res *= res; res = 2. + .607927*res;
      b1  = 1./b1;
      z2 *= -1.;
    }
    if (3*x<1) { // small-x approx.
      res = 0.;
    }
    double b2=b1*b1, b4=b2*b2, b6=b2*b4, b8=b4*b4, b10=b4*b6, b12=b6*b6,
           b3=b1*b2, b5=b1*b4, b7=b1*b6, b9=b1*b8, b11=b1*b10;

    res += ( 1.2158542*b1 + .30396355*b2  + .13509491*b3  + .075990888*b4
           + .04863417*b5 + .03377373*b6  + .02481335*b7  + .018997722*b8
           + .01501055*b9 + .01215854*b10 + .01004838*b11 + .008443432*b12 );
  }
  return .5*z2*res; // = Li2(-x), to absolute err < 1e-7
}

double phat(double x, double A, double B, double C) {
  /*  A = L2A/M2xi
   *  B = L2p/M2xi
   *  C = Fc*as
   */
  double x2 = x*x, onent =  C*( dilog( A/x2 ) - dilog( B/x2 ) );
  return (2./x)*C*exp(onent)*log( (A+x2)/(B+x2) ) ;
}

/*--------------------------------------------------------------------*/
// nuclear params

#define L_p 1.5 // fm, sets 'units'
#define L_eff(A) (double) L_p*1.12*pow(A,1./3.) // HS approx.

double Qs2(double L, double q0, double x2) {
  // broadening
  double xtilde = hbarc/(2.*mp*L);
  return L*q0*pow( .01/fmin(xtilde,x2) , .3 );
}

/*--------------------------------------------------------------------*/
