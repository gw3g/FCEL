#include <stdio.h>
#include <math.h>
#include <string.h>
/*
 *    code for LHC predictions, RpA
 *    @author  GJ
 *    @version 0.1
 *
 */


/*--------------------------------------------------------------------*/
#include "../fcel.h"
#include "sigma.h"
double SQRTS = 8.16; // TeV
double A = 208.;   // Pb
FILE *in;
FILE *out;

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

void RpA_FCEL(Fc,pT,y,sig,params,res)
  double Fc;                                // Ca + CR - Cb
  double pT;                                // transverse mom
  double y;                                 // rapidity
  double (*sig)(double,void*);              // l2_perp(Lp,q0,x2)
  void *params;                             // {u,v,w,x,y,z} = (see below)
  double *res;
{
  /*
   *  nuclear mod factor due to FCEL
   */
  if (fabs(Fc)<tol) {*res=1.; return;} // no mods.

  double alpha_s  = ((double *)params)[0]; // coupling
  double qhat     = ((double *)params)[1]; // transport param.
  double xi       = ((double *)params)[2]; // prob. fraction
  double z        = ((double *)params)[3]; // fragmentation var
  double n        = ((double *)params)[4]; // pp exponent
  double m        = ((double *)params)[5]; // meson mass (?)

  double res_outer, err;

  double kT = pT/z; // partonic momenta rescaled by z
  double mT = sqrt( SQR(pT) + SQR(m) )/z;

  double rs     = SQRTS*1e3, // root s
         x2     = ( mT/((1-xi)*rs) )*exp(-y),
         y_max  = log(rs/mT),
         x_max  = fmin(1.,xi*exp(y_max-y)-1.),
         kappa  = 2.*mT*z/(xi*rs);

  double as   = alpha_s/(2.*M_PI),
         M2   = mT*mT/( xi*(1.-xi) ),
         Q2p  = Qs2(L_p,qhat,x2),
         Q2A  = Qs2(L_eff(A),qhat,x2);

  double _integrand(double x, void *params_out) {
    double A  = ((double *)params_out)[0];
    double B  = ((double *)params_out)[1];
    double C  = ((double *)params_out)[2];

    double  res_integrand,
            dy = SGN(Fc)*log(1.+x),
            P = phat(x,A,B,C);

    double in[2] = {n,kappa}; // ~ (1.-kappa*cosh(y))^n
    return P*sig(y+dy,in)/( sig(y,in)*pow(1.+x,SGN(Fc)) );
  };

  double out[3] = {Q2A/M2,Q2p/M2,fabs(Fc)*as};
  integrator(0.,x_max,_integrand,out,&res_outer,&err);
  *res = res_outer ;

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
// colour-bilities: g, g -> g, g

double rho_27(double xi, double *Fc) { *Fc = 2.*(Nc+1.);
  double xibar = 1.-xi,
         denom = 1.+SQR(xi)+SQR(xibar); denom*=(Nc+1);
  return .5*(Nc+3.)/denom;
}

double rho_1(double xi, double *Fc) { *Fc = 0;
  double dud = 0;
  return rho_27(xi,&dud)*4.*(Nc+1.)/( (Nc+3.)*(SQR(Nc)-1.) );
}

double rho_8(double xi, double *Fc) { *Fc = Nc;
  double dud = 0;
  return 1.-rho_27(xi,&dud)*2.*(Nc+1.)/(Nc+3.);
}

/*--------------------------------------------------------------------*/

int main() {

  double R1, R8, R27;
  double xi=.5, prob;
  double res;
  double Fc;

  double alpha_s = .5; 
  double z =.6;
  // starting params
  double params[6] = {alpha_s,0.07,xi,z,7.,.0};
  RpA_FCEL(3.,2.,-6.,dsig,params,&res);
  printf ("res = % .10f\n", res);

  char *prefix=(char*)"R_gg_gg_pT2GeV";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,".dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pPb, z=%g, alpha=%g\n",z,alpha_s);
  fprintf(out,"# columns: y, R(1), R(8), R(27), R_ave\n");

  double frac;
  double y=-6.;
  double pT=2.;

  while (y<6.) {
    //frac = params[1]/.64;
    frac = (y+6.)/12.;
    res = 0.;

    R1=0.; R8=0.; R27=0.;

    prob = rho_1(xi,&Fc);
    RpA_FCEL(Fc,pT,y,dsig,params,&R1);
    res+= prob*R1;

    prob = rho_8(xi,&Fc);
    RpA_FCEL(Fc,pT,y,dsig,params,&R8);
    res+= prob*R8;

    prob = rho_27(xi,&Fc);
    RpA_FCEL(Fc,pT,y,dsig,params,&R27);
    res+= prob*R27;

    printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e   %.8e\n", y, R1, R8, R27, res );
    y += .1;
  }

  //double A=14.5;
  //printf ("% .10f   % .10f\n", A, L_eff(A) );
  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

  return 0;
}
