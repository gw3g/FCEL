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

void Z_scan_g();
void r_scan_E(double(*)(double,double),double(*)(double));

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

  double as      = ((double *)params)[0]; // coupling/(2*pi)
  double Qs_M    = ((double *)params)[1]; // Qs/M
  double Qs_Qsp  = ((double *)params)[2]; // Qs/Qsp

  double res_outer, err, denom;
  double Qsp_M = Qs_M/Qs_Qsp;

  double _inner(double xF, void *params_inner) {
    // xF ~ x1-x2 (?)
    double x  = ((double *)params_inner)[0];
    double Ep = fabs(E)*(1.+x)/xF;
    return dsig(Ep,xF)*phi(Ep)/(xF*phi(fabs(E)));
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

  double inn[1] = {0.}; // to calc. Zpc w/o FCEL
  integrator(1e-8,1.,_inner,inn,&denom,&err);

  if ((E<0.)||(fabs(Fc)<tol)) { *res = denom; printf("E<0: %.5e\n", denom); return;}
  //else if (fabs(Fc)<tol) {*res=1.; return;} // no mods.

  if (Fc<0) { printf("FCEG still to be implemented\n"); return; }

  double out[3] = {SQR(Qs_M),SQR(Qsp_M),fabs(Fc)*as};
  integrator(0.,x_max,_outer,out,&res_outer,&err);
  *res = res_outer; printf("E>0: %.5e\n", res_outer);
  
}

double Z_sum_( double E, 
               double (*dsig)(double,double), double (*phi)(double),
               void *params) {

  double Fc, temp, prob, res=.0,
         xi  = ((double *)params)[1];

  for (int i=0;i<IRREPS;i++) {
    prob = reps[i](xi,&Fc);
    printf("%g\n",Fc);
    Zpc_FCEL(Fc,E,dsig,phiN,params,&res);
    res+= prob*temp;
  };

  return res;
}

/*--------------------------------------------------------------------*/

int main() {

  Z_scan_g();
  //r_scan_E(dsig_1,phi_H3a);

  return 0;
}

/*--------------------------------------------------------------------*/

void r_scan_E( double (*dsig)(double,double),
               double (*phi)(double)           ) { // SOMETHING WRONG!!
  int N_E;
  double E, Emin, Emax, step,
         Z_orig, Z_fcel;

  double kT = 0., x2 = 1e-7, xi = .5, qhat = .07, m = 1.5;
  double alpha_s = .5, ootp = 1./(2.*M_PI);

  char *prefix=(char*)"r_FCEL_E";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{kT=%.2f,x2=%.1f}.dat",kT,x2);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL (w/scaling), alpha=%g\n",alpha_s);
  fprintf(out,"# columns: E/GeV r (1)\n");


  double res_outer, err;

  double as   = alpha_s*ootp,
         M2   = ( SQR(kT) + SQR(m) )/( xi*(1.-xi) ),
         //Q2p  = Qs2(L_p,qhat,x2);
         Q2A  = Qs2(L_eff(14.5),qhat,x2);

  // starting params
  //
  double params[3] = {as,sqrt(Q2A/M2),sqrt(L_eff(14.5)/L_p)};

  // Here are some parameters that can be changed:
  N_E=50; 
  Emin=1e5;
  Emax=1e8;
  // don't change anything after that.

  double frac;
  E = Emin;
  step=pow(Emax/Emin,1./(double)(N_E-1));

  for (int i=0;i<N_E;i++) { 
    frac = (double)i/(double)(N_E-1);

    printf(" E = %.5e , [%2.2f%]\n", E, 100.*frac); 

    Z_orig = Z_sum_(-E,dsig,phi,params);
    Z_fcel = Z_sum_(+E,dsig,phi,params);

    printf( "%.8e   %.8e   %.8e\n", E, Z_orig, Z_fcel);
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", E, Z_orig, Z_fcel, Z_fcel/Z_orig);

    E *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

void Z_scan_g() {
  double gmin, gmax, dg;

  double R1, R2, R3;
  double E = 1., xi=.5;
  double alpha_s = .5, ootp = 1./(2.*M_PI);
  // starting params
  //
  double params[3] = {alpha_s*ootp,0.01,sqrt(L_eff(14.5)/L_p)};


  char *prefix=(char*)"FCEL_scaling_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{E=%.2f,xi=%.1f}.dat",E,xi);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL (w/scaling), alpha=%g\n",alpha_s);
  fprintf(out,"# columns: Qs/M, F(kappa=.3), F(kappa=.4), F(kappa=.5)\n");

  // Here are some parameters that can be changed:
  gmin=.01;
  gmax=4.;
  dg=.1;
  // don't change anything after that.

  g=gmin;

  printf(" Settings: Qs/Qsp=%g, with g_min=%g, g_max=%g\n",params[2],gmin,gmax); 
  double frac;

  while (g<gmax) {
    frac = g/4.;

    params[1]=.3; R1 = Z_sum_(E,dsig_1,phiN,params);
    params[1]=.4; R2 = Z_sum_(E,dsig_1,phiN,params);
    params[1]=.5; R3 = Z_sum_(E,dsig_1,phiN,params);

    printf(" gamma = %.5e , [%2.2f%]\n", g , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", g, R1, R2, R3 );
    g += dg;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

