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
double x_max = 1000.; // soft-gluon approx. ?
FILE *in;
FILE *out;

void Z_scan_g();
void Z_scan_k();
void r_scan_E(double(*)(double,double),double(*)(double),char*);
void r_errors(double,char*);

/*--------------------------------------------------------------------*/

#include <gsl/gsl_integration.h>
size_t calls=1e7; double tol=1e-3;

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

  double _inner(double x, void *params_inner) {

    double zz = pow(1.+x,-SGN(Fc)); // rescaling variable:
                                    // x = eps/E, zz = 1/(x+1) OR zz = (x+1)

    double A  = ((double *)params_inner)[0];
    double B  = ((double *)params_inner)[1];
    double C  = ((double *)params_inner)[2];
    double xF = ((double *)params_inner)[3];

    double P = phat(x, A, B, C); // quenchin weight
    double Ep = fabs(E)/(xF);    // proton energy

    if ( Ep > 9.832e10 ) return 0.; // GZK cutoff (?)

    return dsig(Ep,xF/zz)*P;
  };

  double _outer(double xF, void *params_outer) {
    double A  = ((double *)params_outer)[0];
    double B  = ((double *)params_outer)[1];
    double C  = ((double *)params_outer)[2];

    double Ep = fabs(E)/(xF);

    double  res_inner;

    x_max = 1./xF - 1.; // xF < 1

    double in[4] = { A, B, C, xF }; 

    if ((E<0)||(fabs(Fc)<1e-1)) { res_inner = dsig(Ep,xF); }
    else if (E>0) {
      integrator(0.,x_max,_inner,in,&res_inner,&err);
    }

    return res_inner*phi(Ep)/(xF*phi(fabs(E)));
  };
  double out[3] = {SQR(Qs_M),SQR(Qsp_M),fabs(Fc)*as};

  integrator(0.,1.,_outer,out,&res_outer,&err);

  *res = res_outer; //printf("res: %.5e , err: %.5e\n", res_outer, err);
}

double Z_sum_( double E, // neutrino energy (in GeV)
               double (*dsig)(double,double), double (*phi)(double),
               void *params) {

  double Fc, temp, prob, res=.0,
         xi  = ((double *)params)[3];

  //int i=0;
  for (int i=0;i<IRREPS;i++) {
    prob = reps[i](xi,&Fc);
    Zpc_FCEL(Fc,E,dsig,phi,params,&temp);
    //printf(" i = %d, temp = %g\n", i, temp);
    res += prob*temp;
  };

  return res;
}

double rFCEL_scaling(void *params) {
  double dummy = 1.;
  return Z_sum_(dummy,dsig_1,phi_SP,params)/Z_sum_(-dummy,dsig_1,phi_SP,params);
}

/*--------------------------------------------------------------------*/
// Hessian method
//
#define PARAMS 4  //  q0, xi, Kt, m_Q
double dS[PARAMS] = {.02,.25,.5,0.2};
double  S[PARAMS] = {.07,.50,2.,1.3};

double R_nu(double E, void *params) {
  double Fc, temp, prob, res=.0,
         alpha_s = .5,
         qhat= ((double *)params)[0],
         xi  = ((double *)params)[1],
         kT  = ((double *)params)[2],
         m   = ((double *)params)[3];

  double as   = alpha_s*ootp,
         M2   = ( SQR(kT) + SQR(m) )/( xi*(1.-xi) ),
         Q2A  ;//= Qs2(L_eff(14.5),qhat,x2);

  double x2 = ( SQR(kT) + SQR(m) )/(mp*E); // proper def. (?)
  Q2A  = Qs2(L_eff(14.5),qhat,x2);

  double new_params[4];
  new_params[0] = as;
  new_params[1] = sqrt(Q2A/M2);
  new_params[2] = sqrt(L_eff(14.5)/L_p);
  new_params[3] = xi;

  res = rFCEL_scaling(new_params);

  return res;
}


void R_limits(double E, 
              double *R_ave, double *R_min, double *R_max) {

  double Sp[PARAMS], Sm[PARAMS], res1=0.,res2=0.,
         RSp,        RSm,        RS = R_nu(E,S);

  memcpy(Sp,S,sizeof(double)*PARAMS);
  memcpy(Sm,S,sizeof(double)*PARAMS);

  for (int k=0; k<PARAMS; k++) {

    // modify running arrays ...
    Sp[k]+=dS[k]; Sm[k]-=dS[k];

    RSp = R_nu(E,Sp); RSm = R_nu(E,Sm);

    res1+= SQR( fmax( fmax( RSp-RS, RSm-RS ), 0. ) );
    res2+= SQR( fmax( fmax( RS-RSp, RS-RSm ), 0. ) );

    // ... undo those mods.
    Sp[k]-=dS[k]; Sm[k]+=dS[k];
  };

  *R_ave = RS;
  *R_max = RS + sqrt(res1);
  *R_min = RS - sqrt(res2);
}


/*--------------------------------------------------------------------*/

int main() {
  //channel(1); // gg - gg
  //channel(2); // qg - qg
  //channel(4); // gg - qq
  channel(5); // qq - QQ
  double Fc;
  printf("IRREPS = %d\n", IRREPS );

  //printf("p = %g\n", phat(2.2,1,.5,.3) );

  //Z_scan_g();
  //Z_scan_k();
  //r_scan_E(dsig_3,phi_H3a,"out/Rnu_prompt_H3a_");
  //r_scan_E(dsig_3,phi_knee,"out/Rnu_prompt_knee_");
  //r_scan_E(dsig_3,phi_GSF,"out/Rnu_prompt_GSF_");
  lam = .3;
  g = 2.0;
  r_scan_E(dsig_1,phi_SP,"out/Rnu_prompt_scaling0_");
  g = 2.7;
  r_scan_E(dsig_1,phi_SP,"out/Rnu_prompt_scaling1_");
  g = 3.;
  r_scan_E(dsig_1,phi_SP,"out/Rnu_prompt_scaling2_");
  g = 3.4;
  r_scan_E(dsig_1,phi_SP,"out/Rnu_prompt_scaling3_");
  g = 3.6;
  r_scan_E(dsig_1,phi_SP,"out/Rnu_prompt_scaling4_");//*/
  //g = 3.0;
  //r_errors(g,"out/Rnu_uncertainty_scaling_");

  return 0;
}

/*--------------------------------------------------------------------*/

void r_scan_E( double (*dsig)(double,double),
               double (*phi)(double),
               char    *prefix                ) {
  int N_E;
  double E, Emin, Emax, step,
         Z_orig, Z_fcel;

  double kT = 2., x2 = 1e-5, xi = .5, qhat = .07, m = 1.3;
  //double kT = 2., x2 = 1e-5, xi = .5, qhat = .07, m = 0.0;
  double alpha_s = .5;

  //char *prefix=(char*)"r_FCEL2_E";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  //sprintf(suffix,"_{kT=%.2f,x2=%g}.dat",kT,x2);
  sprintf(suffix,"_{kT=%.2f}.dat",kT);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL, xi=%g, alpha=%g\n",xi,alpha_s);
  fprintf(out,"# columns: E/GeV,    Zpc,    Zpc w/ FCEL,   ratio\n");


  double res_outer, err;

  double as   = alpha_s*ootp,
         M2   = ( SQR(kT) + SQR(m) )/( xi*(1.-xi) ),
         Q2A  ;//= Qs2(L_eff(14.5),qhat,x2);

  double params[4]; //= {as,sqrt(Q2A/M2),sqrt(L_eff(14.5)/L_p),xi};
  params[0] = as;
  params[2] = sqrt(L_eff(14.5)/L_p);
  params[3] = .5;

  // Here are some parameters that can be changed:
  N_E=100; 
  Emin=1e2;
  Emax=1e10;
  // don't change anything after that.

  printf("Settings: kT = %2.2f, x2 = %g \n",kT,x2); 
  double frac;

  E = Emin;
  step=pow(Emax/Emin,1./(double)(N_E-1));

  for (int i=0;i<N_E;i++) { 
    frac = (double)i/(double)(N_E-1);

    printf(" E = %.5e , [%2.2f%]\n", E, 100.*frac); 

    x2 = M2/(4.*mp*E); // proper def. (?)
    Q2A  = Qs2(L_eff(14.5),qhat,x2);
    params[1] = sqrt(Q2A/M2);
    Z_orig = Z_sum_(-E,dsig,phi,params);
    Z_fcel = Z_sum_(+E,dsig,phi,params);

    //printf( "%.8e   %.8e   %.8e\n", E, Z_orig, Z_fcel);
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", E, Z_orig, Z_fcel, Z_fcel/Z_orig);

    E *= step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

void r_errors(double ga, char *prefix) {
  g=ga;
  int N_E;
  double E, Emin, Emax, step,
         R, Rmin, Rmax;

  double alpha_s = .5;

  //char *prefix=(char*)"r_FCEL2_E";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,".dat");
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# FCEL scaling, gamma=%g, [Hessian Method]\n",ga);
  fprintf(out,"# columns: E/GeV, R, R_min, R_max, R(qhat-), R(qhat+), R(xi-), R(xi+), R(kt-), R(kt+), R(m-), R(m+) \n");

  double res_outer, err;

  // Here are some parameters that can be changed:
  N_E=100; 
  Emin=1e2;
  Emax=1e10;
  // don't change anything after that.

  printf("Settings: gamma = %2.2f \n",ga); 
  double frac;

  E = Emin;
  step=pow(Emax/Emin,1./(double)(N_E-1));

  for (int i=0;i<N_E;i++) { 
    frac = (double)i/(double)(N_E-1);

    printf(" E = %.5e , [%2.2f%]\n", E, 100.*frac); 

    R_limits(E,&R,&Rmin,&Rmax); // calculation

    //printf( "%.8e   %.8e   %.8e\n", E, Z_orig, Z_fcel);
    fprintf( out, "%.8e   %.8e   %.8e   %.8e", E, R, Rmin, Rmax );

    for (int j=0; j<4; j++) {
      S[j] += dS[j];
      Rmin = R_nu(E,S); 
      S[j] -= 2.*dS[j];
      Rmax = R_nu(E,S); 
      S[j] += dS[j];
      fprintf( out, "   %.8e   %.8e", Rmin, Rmax);
    }
    fprintf( out, "\n");


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

