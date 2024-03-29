#include <math.h>
#define SQR(x) ((x)*(x))
#define SGN(x) (double) ((x>0)-(x<0))

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
double ootp   = .15915494309189533577 ; // 1/(2.pi)

/*--------------------------------------------------------------------*/
typedef double (*rho)();

// colour-bilities: g, g -> g, g

#define C27  2.*(Nc+1.)

double gg_gg_27(double xi, double *Fc) { *Fc = C27;
  double xibar = 1.-xi,
         denom = 1.+SQR(xi)+SQR(xibar); denom*=(Nc+1.);
  return .5*(Nc+3.)/denom;
}

double gg_gg_1(double xi, double *Fc) { *Fc = 0.;
  double dud = 0;
  return gg_gg_27(xi,&dud)*4./( (Nc+3.)*(Nc-1.) );
}

double gg_gg_8(double xi, double *Fc) { *Fc = Ca;
  double dud = 0;
  return 1.-gg_gg_27(xi,&dud)*2.*(Nc+1.)/(Nc+3.);
}
// */
// colour-bilities: q, g -> q, g

#define C15  .5*(Nc+1.)*(3.*Nc-1.)/Nc
#define C6   .5*(Nc-1.)*(3.*Nc+1.)/Nc

double qg_qg_15(double xi, double *Fc) { *Fc = Cf + C15 - Ca;
  double xibar = 1.-xi,
         denom = Cf*SQR(xi)+Nc*xibar; denom*=4.*(Nc+1.);
  return Nc*(Nc+2.)/denom;
}

double qg_qg_6(double xi, double *Fc) { *Fc = Cf + C6 - Ca;
  double dud   = 0,
         denom = (Nc-1.)*(Nc+2.);
  return qg_qg_15(xi,&dud)*(Nc+1.)*(Nc-2.)/denom;
}

double qg_qg_3(double xi, double *Fc) { *Fc = 2.*Cf - Ca ;
  double dud   = 0,
         denom = (Nc+2.)*Nc;
  //return qg_qg_15(xi,&dud)*4.*(Nc+1.)*Cf*SQR(xi-.5*Nc/Cf)/denom; // ???
  return 1. - qg_qg_6(xi,&dud) - qg_qg_15(xi,&dud);
}
// */
// colour-bilities: g, q -> q, g

double gq_qg_15(double xi, double *Fc) { *Fc = Ca + C15 - Cf;
  double dud   = 0;
  return qg_qg_15(xi,&dud);
}

double gq_qg_6(double xi, double *Fc) { *Fc = Ca + C6 - Cf;
  double dud   = 0;
  return qg_qg_6(xi,&dud);
}

double gq_qg_3(double xi, double *Fc) { *Fc = Ca ;
  double dud   = 0;
  return qg_qg_3(xi,&dud);
}

// */
// colour-bilities: g, g -> q, \bar{q}

double gg_qq_1(double xi, double *Fc) { *Fc = 0;
  double xibar = 1.-xi,
         denom = SQR(Nc*xi) + SQR(Nc*xibar) - 1.;
  return 1./denom;
}

double gg_qq_8(double xi, double *Fc) { *Fc = Ca;
  double dud   = 0;
  return 1.-gg_qq_1(xi,&dud);
}
// */
//
// colour-bilities: q, q -> Q, \bar{Q}

double qq_qq_8(double xi, double *Fc) { *Fc = Ca;
  return 1.;
}

// */


//#define IRREPS 2
rho reps_gg_gg[3] = {&gg_gg_1,&gg_gg_8,&gg_gg_27};
rho reps_qg_qg[3] = {&qg_qg_3,&qg_qg_6,&qg_qg_15};
rho reps_gq_qg[3] = {&gq_qg_3,&gq_qg_6,&gq_qg_15};
rho reps_gg_qq[2] = {&gg_qq_1,&gg_qq_8};
rho reps_qq_qq[1] = {&qq_qq_8};
rho *reps;
char reaction[10];
int IRREPS;

void channel(int r) {
  switch(r) {
    case 1 : strcpy(reaction, "{gg;gg}"); reps = reps_gg_gg; IRREPS = 3; break;
    case 2 : strcpy(reaction, "{qg;qg}"); reps = reps_qg_qg; IRREPS = 3; break;
    case 3 : strcpy(reaction, "{gq;qg}"); reps = reps_gq_qg; IRREPS = 3; break;
    case 4 : strcpy(reaction, "{gg;qq}"); reps = reps_gg_qq; IRREPS = 2; break;
    case 5 : strcpy(reaction, "{qq;qq}"); reps = reps_qq_qq; IRREPS = 1; break;
    default: reps = NULL;
  }
  printf("Switched channel to %d\n",r);
}


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

double lam = .3;

double Qs2(double L, double q0, double x2) {
  // broadening [GeV^2]
  double xtilde = hbarc/(2.*mp*L);
  return L*q0*pow( .01/fmin(xtilde,x2) , lam );
}

/*--------------------------------------------------------------------*/
