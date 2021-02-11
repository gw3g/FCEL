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
#include "../fcel.h"
#include "sigma.h"
#include "crflux.h"
double g; // spectral index
double x_max = 100.; // soft-gluon approx. ?
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
