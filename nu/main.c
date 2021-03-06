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
double x_max = 1.; // soft-gluon approx. ?
FILE *in;
FILE *out;

void Z_scan_g();
void Z_scan_k();
void r_scan_E(double(*)(double,double),double(*)(double),char*);

/*--------------------------------------------------------------------*/

#include <gsl/gsl_integration.h>
size_t calls=1e5; double tol=1e-4;

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
    double zz  = ((double *)params_inner)[0];
    double Ep = fabs(E)/(zz*xF);
    return dsig(Ep,xF)*phi(Ep)/(xF*phi(fabs(E)));
  };

  double _outer(double x, void *params_outer) {
    double A  = ((double *)params_outer)[0];
    double B  = ((double *)params_outer)[1];
    double C  = ((double *)params_outer)[2];

    double  res_inner, P;

    if (Fc>0)      { P = phat( x, A, B, C); }
    else if (Fc<0) { P = phat( x, A, B, C); }
    //printf(" zz = %g, P = %g\n",zz,P);

    double zz = pow(1.+x,-SGN(Fc));
    double in[1] = { zz }; // x = eps/E, zz = 1/(x+1)
    integrator(1e-8,fmin(1.,1./zz),_inner,in,&res_inner,&err);
    //return P;
    return res_inner*P;
  };

  double inn[1] = {1.}; // to calc. Zpc w/o FCEL
  integrator(1e-8,1.,_inner,inn,&denom,&err);

  if ((E<0.)||(fabs(Fc)<tol)) { *res = denom; /*printf("E<0: %.5e\n", denom);*/ return;}
  //else if (fabs(Fc)<tol) {*res=1.; return;} // no mods.

  //if (Fc<0) { printf("FCEG still to be implemented\n"); return; }

  double out[3] = {SQR(Qs_M),SQR(Qsp_M),fabs(Fc)*as};
  if (Fc>0) { integrator(.0,10.,_outer,out,&res_outer,&err); }
  else if (Fc<0) { integrator(.0,10.,_outer,out,&res_outer,&err); }
  *res = res_outer; //printf("res: %.5e\n", res_outer);
}

double Z_sum_( double E, // neutrino energy (in GeV)
               double (*dsig)(double,double), double (*phi)(double),
               void *params) {

  double Fc, temp, prob, res=.0,
         xi  = ((double *)params)[1];

  //int i=0;
  for (int i=0;i<IRREPS;i++) {
    prob = reps[i](xi,&Fc);
    Zpc_FCEL(Fc,E,dsig,phi,params,&temp);
    res += prob*temp;
  };

  return res;
}

double rFCEL_scaling(void *params) {
  double dummy = 1.;
  return Z_sum_(dummy,dsig_1,phi_SP,params)/Z_sum_(-dummy,dsig_1,phi_SP,params);
}

/*--------------------------------------------------------------------*/

int main() {
  channel(2);
  double Fc;
  printf("IRREPS = %d\n", IRREPS );

  //printf("p = %g\n", phat(2.2,1,.5,.3) );

  //Z_scan_g();
  //Z_scan_k();
  //r_scan_E(dsig_1,phi_H3a,"out/r_nu_FCEL_H3a_");
  //r_scan_E(dsig_1,phi_knee,"out/r_nu_FCEL_knee_");
  //r_scan_E(dsig_1,phi_GSF,"out/r_nu_FCEL_GSF_");
  g = 2.7;
  r_scan_E(dsig_1,phi_SP,"out/r_nu_scaling1_");
  g = 3.;
  r_scan_E(dsig_1,phi_SP,"out/r_nu_scaling2_");
  g = 3.4;
  r_scan_E(dsig_1,phi_SP,"out/r_nu_scaling3_");//*/

  return 0;
}

/*--------------------------------------------------------------------*/

void r_scan_E( double (*dsig)(double,double),
               double (*phi)(double),
               char    *prefix                ) {
  int N_E;
  double E, Emin, Emax, step,
         Z_orig, Z_fcel;

  double kT = .1, x2 = 1e-5, xi = .5, qhat = .07, m = 0.0;
  double alpha_s = .5;

  //char *prefix=(char*)"r_FCEL2_E";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{kT=%.2f,x2=%g}.dat",kT,x2);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL, xi=%g, alpha=%g\n",xi,alpha_s);
  fprintf(out,"# columns: E/GeV,    Zpc,    Zpc w/ FCEL,   ratio\n");


  double res_outer, err;

  double as   = alpha_s*ootp,
         M2   = ( SQR(kT) + SQR(m) )/( xi*(1.-xi) ),
         Q2A  = Qs2(L_eff(14.5),qhat,x2);

  double params[3] = {as,sqrt(Q2A/M2),sqrt(L_eff(14.5)/L_p)};

  // Here are some parameters that can be changed:
  N_E=400; 
  Emin=1e3;
  Emax=1e10;
  // don't change anything after that.

  printf("Settings: kT = %2.2f, x2 = %g \n",kT,x2); 
  double frac;

  E = Emin;
  step=pow(Emax/Emin,1./(double)(N_E-1));

  for (int i=0;i<N_E;i++) { 
    frac = (double)i/(double)(N_E-1);

    printf(" E = %.5e , [%2.2f%]\n", E, 100.*frac); 

    Z_orig = Z_sum_(-E,dsig,phi,params);
    Z_fcel = Z_sum_(+E,dsig,phi,params);

    //printf( "%.8e   %.8e   %.8e\n", E, Z_orig, Z_fcel);
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


  char *prefix=(char*)"out/rFCEL_gamma_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{E=%.2f,xi=%.1f}.dat",E,xi);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL (w/scaling, kappa = Qs/M), alpha=%g\n",alpha_s);
  fprintf(out,"# columns: gamma, F(kappa=.3), F(kappa=.4), F(kappa=.5)\n");

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

    params[1]=.3; R1 = rFCEL_scaling(params);
    params[1]=.4; R2 = rFCEL_scaling(params);
    params[1]=.5; R3 = rFCEL_scaling(params);

    printf(" gamma = %.5e , [%2.2f%]\n", g , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", g, R1, R2, R3 );
    g += dg;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


void Z_scan_k() {
  double kappa, k_min, k_max, dk;

  double R1, R2, R3;
  double E = 1., xi=.5;
  double alpha_s = .5;
  // starting params
  //
  double params[3] = {alpha_s*ootp,0.,sqrt(L_eff(14.5)/L_p)};


  char *prefix=(char*)"out/rFCEL_kappa_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{E=%.2f,xi=%.1f}.dat",E,xi);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL (w/scaling, kappa = Qs/M), alpha=%g\n",alpha_s);
  fprintf(out,"# columns: kappa, F(gamma=2), F(gamma=2.5), F(gamma=3)\n");

  // Here are some parameters that can be changed:
  k_min = .001;
  k_max = .6 ;
  dk    = .01 ;
  // don't change anything after that.

  kappa = k_min;

  printf(" Settings: gamma=%g, with g_min=%g, g_max=%g\n",g,k_min,k_max); 
  double frac;

  while (kappa<k_max) {
    frac = kappa/(k_max-k_min);
    params[1] = kappa;

    g = 2.0; R1 = rFCEL_scaling(params);
    g = 2.5; R2 = rFCEL_scaling(params);
    g = 3.0; R3 = rFCEL_scaling(params);

    printf(" Qs/M = %.5e , [%2.2f%]\n", kappa , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", kappa, R1, R2, R3 );
    kappa += dk;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

