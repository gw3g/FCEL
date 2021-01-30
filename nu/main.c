#include <stdio.h>
#include <math.h>
#include <string.h>
/*
 *    code for neutrino flux w/ FCEL
 *    @author  GJ
 *    @version 0.1
 *
 */


/*--------------------------------------------------------------------*/
#include "sigma.h"
#include "crflux.h"
double g; // spectral index
double x_max = 100.; // soft-gluon approx. ?
FILE *in;
FILE *out;

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

#define SQR(x) x*x
#define SGN(x) (double) ((x>0)-(x<0))

/*--------------------------------------------------------------------*/

#include <gsl/gsl_integration.h>
size_t calls=1e5; double tol=1e-7;

void integrator(a,b,func,params,res,err)
  double a, b;
  double (*func)(double,void *);
  void (*params);
  double *res, *err;
{
  /*
   *  generic code: integrate func(x) from x=a to b
   */
  size_t limit = calls;
  gsl_integration_workspace * WS = gsl_integration_workspace_alloc (calls);
  gsl_function f_aux; f_aux.function = func; f_aux.params = params;
  gsl_integration_qag (&f_aux, a, b, tol, tol, calls, 3, WS, res, err); 
  gsl_integration_workspace_free (WS);
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

void Zpc_FCEL(Fc,E,dsig,phi,params,res)
  double Fc;                      // Ca + CR - Cb
  double E;                       // neutrino energy
  double (*dsig)(double,double);  // dsig/dxF(E,xF)
  double (*phi)(double);          // flux_CR(E)
  void *params;                   // {x,y,z} = (see below)
  double *res;
{
  /*
   *  Z_pc moment w/ FCEL effect
   */
  if (fabs(Fc)<tol) {*res=1.; return;} // no mods.
  if (Fc<0) { printf("FCEG still to be implemented\n"); return; }

  double as      = ((double *)params)[0]; // coupling/(2*pi)
  double Qs_M    = ((double *)params)[1]; // Qs/M
  double Qs_Qsp  = ((double *)params)[2]; // Qs/Qsp

  double res_outer, err, denom;
  double Qsp_M = Qs_M/Qs_Qsp;

  double _inner(double xF, void *params_inner) {
    // xF ~ x1-x2 (?)
    double x  = ((double *)params_inner)[0];
    double Ep = E*(1.+x)/xF;
    return dsig(Ep,xF)*phi(Ep)/(xF*phi(E));
  };

  double _outer(double x, void *params_outer) {
    double A  = ((double *)params_outer)[0];
    double B  = ((double *)params_outer)[1];
    double C  = ((double *)params_outer)[2];

    double  res_inner,
            P = phat(x,A,B,C);

    double in[1] = {x}; // x = eps/E
    integrator(1e-8,1.,_inner,in,&res_inner,&err);
    return res_inner*P;
  };

  double out[3] = {SQR(Qs_M),SQR(Qsp_M),fabs(Fc)*as};
  integrator(0.,x_max,_outer,out,&res_outer,&err);
  out[0] = 0.; // to calc. Zpc w/o FCEL
  integrator(1e-8,1.,_inner,out,&denom,&err);
  *res = res_outer/denom;

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


double prepare_input(double Fc, void *params) {
  double alpha_s  = ((double *)params)[0]; // coupling
  double qhat     = ((double *)params)[1]; // transport param.
  double x2       = ((double *)params)[2]; // x2 (?)
  double xi       = ((double *)params)[3]; // prob. fraction
  double mT       = ((double *)params)[4]; // transverse mass
  double A        = ((double *)params)[5]; // 

  double as   = alpha_s/(2.*M_PI),
         M2   = mT*mT/( xi*(1.-xi) ),
         Q2p  = Qs2(L_p,qhat,x2),
         Q2A  = Qs2(L_eff(A),qhat,x2);

  double Qs_M   = sqrt(Q2A/M2); // Qs/M
  double Qs_Qsp = sqrt(Q2A/Q2p); // Qs/Qsp
}

/*--------------------------------------------------------------------*/
// 1 -> 2 kinematics
void ymT_to_x1x2(double y,double pT,double *in, double *out) {
  double s       = ((double *)in)[0];
  double xi      = ((double *)in)[1];
}
void x1x2_to_ymT(double x1,double x2,double *out) {
}

/*--------------------------------------------------------------------*/
// colour-bilities: g, g -> c, cbar

double rho_1(double xi, double *Fc) { *Fc = 0;
  double xibar = 1.-xi,
         denom = SQR(Nc)*( SQR(xi)+SQR(xibar) ) - 1.;
  return 1./denom;
}

double rho_8(double xi, double *Fc) { *Fc = Nc;
  double dud = 0;
  return 1.-rho_1(xi,&dud);
}

/*--------------------------------------------------------------------*/

int main() {

  double R1, R2, R3;
  double E = 1., xi=.5, prob;
  double res;
  double Fc;

  double alpha_s = .5, ootp = 1./(2.*M_PI);
  // starting params
  double params[3] = {alpha_s*ootp,0.001,sqrt(L_eff(14.5)/L_p)};
  printf ("Qs/Qsp = % .10f\n", params[2] );

  char *prefix=(char*)"FCEL_scaling2";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,".dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL (w/scaling), alpha=%g\n",alpha_s);
  fprintf(out,"# columns: Qs/M, F(kappa=.3), F(kappa=.4), F(kappa=.5)\n");

  double frac;
  g=.01;

  while (g<4.) {
    //frac = params[1]/.64;
    frac = g/4.;

    R1=0.; R2=0.; R3=0.;

    params[1]=.3;
    prob = rho_1(xi,&Fc);
    Zpc_FCEL(Fc,E,dsig_1,phiN,params,&res);
    R1+= prob*res;
    prob = rho_8(xi,&Fc);
    Zpc_FCEL(Fc,E,dsig_1,phiN,params,&res);
    R1 = res;

    //g=2.5;
    params[1]=.4;
    prob = rho_1(xi,&Fc);
    Zpc_FCEL(Fc,E,dsig_1,phiN,params,&res);
    R2+= prob*res;
    prob = rho_8(xi,&Fc);
    Zpc_FCEL(Fc,E,dsig_1,phiN,params,&res);
    R2 = res;

    //g=3.;
    params[1]=.5;
    prob = rho_1(xi,&Fc);
    Zpc_FCEL(Fc,E,dsig_1,phiN,params,&res);
    R3+= prob*res;
    prob = rho_8(xi,&Fc);
    Zpc_FCEL(Fc,E,dsig_1,phiN,params,&res);
    R3 = res;

    printf(" Qs/M = %.5e , [%2.2f%]\n", params[1] , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", g, R1, R2, R3 );
    g += .1;
  }

  //double A=14.5;
  //printf ("% .10f   % .10f\n", A, L_eff(A) );
  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

  return 0;
}
