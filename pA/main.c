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
double alpha_s = 0.5; // coupling

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

  double qhat     = ((double *)params)[0]; // transport param.
  double xi       = ((double *)params)[1]; // prob. fraction
  double z        = ((double *)params)[2]; // fragmentation var
  double n        = ((double *)params)[3]; // pp exponent
  double m        = ((double *)params)[4]; // meson mass (?)

  double res_outer, err;

  double kT = pT/z; // partonic momenta rescaled by z
  double mT = sqrt( SQR(kT) + SQR(m) );

  double rs     = SQRTS*1e3, // root s
         x2     = ( mT/((1.-xi)*rs) )*exp(-y),
         y_max  = log(rs/mT),
         x_max  = fmin(1.,xi*exp(y_max-y)-1.),
         kappa  = 2.*mT*z/(.5*rs);

  double as   = alpha_s/(2.*M_PI),
         M2   = mT*mT/( xi*(1.-xi) ),
         Q2p  = Qs2(L_p,qhat,x2),
         Q2A  = Qs2(L_eff(A),qhat,x2);

  double _integrand(double x, void *params_out) {
    double A  = ((double *)params_out)[0];
    double B  = ((double *)params_out)[1];
    double C  = ((double *)params_out)[2];

    double  dy = SGN(Fc)*log(1.+x),
            P  = phat(x,A,B,C);

    double in[2] = {n,kappa}; // ~ (1.-kappa*cosh(y))^n
    return P*sig(y+dy,in)/( sig(y,in)*pow(1.+x,SGN(Fc)) );
  };

  double out[3] = {Q2A/M2,Q2p/M2,fabs(Fc)*as};
  integrator(0.,x_max,_integrand,out,&res_outer,&err);
  *res = res_outer ;

}

/*--------------------------------------------------------------------*/
// colour-bilities: g, g -> g, g

double rho_27(double xi, double *Fc) { *Fc = 2.*(Nc+1.);
  double xibar = 1.-xi,
         denom = 1.+SQR(xi)+SQR(xibar); denom*=(Nc+1.);
  return .5*(Nc+3.)/denom;
}

double rho_1(double xi, double *Fc) { *Fc = 0.;
  double dud = 0;
  return rho_27(xi,&dud)*4./( (Nc+3.)*(Nc-1.) );
}

double rho_8(double xi, double *Fc) { *Fc = Nc;
  double dud = 0;
  return 1.-rho_27(xi,&dud)*2.*(Nc+1.)/(Nc+3.);
}

typedef double (*rho)();

#define IRREPS 3
rho reps[IRREPS] = {&rho_1,&rho_8,&rho_27};

/*--------------------------------------------------------------------*/
// Hessian method
//
#define PARAMS 5  // q0, xi, z, n,  m_{Meson?}
double dS[PARAMS] = {.02,.25,.2,2.5,0.0};
double  S[PARAMS] = {.07,.50,.6,7.5,0.0};

double R_sum_(double pT, double y, void *params) {
  double Fc, temp, prob, res=.0,
         xi  = ((double *)params)[1];

  for (int i=0;i<IRREPS;i++) {
    prob = reps[i](xi,&Fc);
    //printf("%g\n",Fc);
    RpA_FCEL(Fc,pT,y,dsig,params,&temp);
    res+= prob*temp;
  };

  return res;
}

void R_limits(double pT, double y, double *R_ave, double *R_min, double *R_max) {
  double Sp[PARAMS], Sm[PARAMS], res1=0.,res2=0.,
         RSp,        RSm,        RS = R_sum_(pT,y,S);

  memcpy(Sp,S,sizeof(double)*PARAMS);
  memcpy(Sm,S,sizeof(double)*PARAMS);

  for (int k=0; k<PARAMS; k++) {

    // modify running arrays ...
    Sp[k]+=dS[k]; Sm[k]-=dS[k];

    RSp = R_sum_(pT,y,Sp); RSm = R_sum_(pT,y,Sm);

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
  double R, Rmin, Rmax;

  char *prefix=(char*)"R_pT6GeV";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,".dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");
  //fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  fprintf(out,"# columns: y, R_ave, R_min, R_max\n");

  double frac;
  double y=-6.;
  double pT=6.;

  while (y<6.) {
    frac = (y+6.)/12.;

    R_limits(pT,y,&R,&Rmin,&Rmax);

    //printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", y, R, Rmin, Rmax );
    y += .1;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

  return 0;
}

/*--------------------------------------------------------------------*/
int oldmain() {

  double R[3]; // 1;8;27
  double xi=S[1], prob;
  double res;
  double Fc;

  // starting params
  //RpA_FCEL(3.,2.,-6.,dsig,params,&res);
  //printf ("res = % .10f\n", res);

  char *prefix=(char*)"R_gg_gg_pT2GeV";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,".dat");
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  fprintf(out,"# columns: y, R(1), R(8), R(27), R_ave\n");

  double frac;
  double y=-6.;
  double pT=2.;

  while (y<6.) {
    frac = (y+6.)/12.;
    res = 0.;

    for (int i=0;i<3;i++) {
      prob = reps[i](xi,&Fc);
      RpA_FCEL(Fc,pT,y,dsig,S,&R[i]);
      res+= prob*R[i];
    }

    printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e   %.8e\n", y, R[0], R[1], R[2], res );
    y += .1;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

  return 0;
}
