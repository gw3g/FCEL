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

//FILE *in;
FILE *out;
void R_reps(double);   // pT input
void R_scan_y(double) ;// pT fixed
void R_scan_pT(double);// y  fixed

/*--------------------------------------------------------------------*/

#include <gsl/gsl_integration.h>
size_t calls=1e5; double tol=1e-2;

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
  //printf("test: %.2f\n",pT);
  integrator(0.,x_max,_integrand,out,&res_outer,&err);
  *res = res_outer ;

}

/*--------------------------------------------------------------------*/
/*typedef double (*rho)();

// colour-bilities: g, g -> g, g

#define C27  2.*(Nc+1.)
double rho_27(double xi, double *Fc) { *Fc = C27;
  double xibar = 1.-xi,
         denom = 1.+SQR(xi)+SQR(xibar); denom*=(Nc+1.);
  return .5*(Nc+3.)/denom;
}

double rho_1(double xi, double *Fc) { *Fc = 0.;
  double dud = 0;
  return rho_27(xi,&dud)*4./( (Nc+3.)*(Nc-1.) );
}

double rho_8(double xi, double *Fc) { *Fc = Ca;
  double dud = 0;
  return 1.-rho_27(xi,&dud)*2.*(Nc+1.)/(Nc+3.);
}
// */
// colour-bilities: q, g -> q, g
/*
#define C15  .5*(Nc+1.)*(3.*Nc-1.)/Nc
#define C6   .5*(Nc-1.)*(3.*Nc+1.)/Nc

double rho_15(double xi, double *Fc) { *Fc = Cf + C15 - Ca;
  double xibar = 1.-xi,
         denom = Cf*SQR(xi)+Nc*xibar; denom*=4.*(Nc+1.);
  return Nc*(Nc+2.)/denom;
}

double rho_6(double xi, double *Fc) { *Fc = Cf + C6 - Ca;
  double dud   = 0,
         denom = (Nc-1.)*(Nc+2.);
  return rho_15(xi,&dud)*(Nc+1.)*(Nc-2.)/denom;
}

double rho_3(double xi, double *Fc) { *Fc = 2.*Cf - Ca ;
  double dud   = 0,
         denom = (Nc+2.)*Nc;
  return rho_15(xi,&dud)*4.*(Nc+1.)*Cf*SQR(xi-.5*Nc/Cf)/denom;
}
// */
// colour-bilities: g, g -> q, \bar{q}
/*
double rho_1(double xi, double *Fc) { *Fc = 0;
  double xibar = 1.-xi,
         denom = SQR(Nc*xi) + SQR(Nc*xibar) - 1.;
  return 1./denom;
}

double rho_8(double xi, double *Fc) { *Fc = Ca;
  double dud   = 0;
  return 1.-rho_1(xi,&dud);
}
// */

//#define IRREPS 3
//rho reps[IRREPS] = {&rho_1,&rho_8,&rho_27};
//rho reps[IRREPS] = {&rho_3,&rho_6,&rho_15};
//rho reps[IRREPS] = {&rho_1,&rho_8};

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

   R_reps(2.);
   //R_scan_pT(0.);
   R_scan_pT(-2.);
   R_scan_y(2.);
   //R_scan_pT(4.);
   return 0;
}

/*--------------------------------------------------------------------*/
void R_reps(double pT) {
  int N_y;
  double y, y_min, y_max, step;

  double xi=S[1], R[IRREPS];
  double prob;
  double res;
  double Fc;

  char *prefix=(char*)"out/R_pT";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"%gGeV.dat",pT);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  fprintf(out,"# columns: y, {R1,R2,...}, R_ave\n");

  // Here are some parameters that can be changed:
  N_y=50; 
  y_min=-6.;
  y_max=+6.;
  // don't change anything after that.

  step=(y_max-y_min)/((double) N_y-1);
  y=y_min;

  printf(" Settings: pT=%g, with y_min=%g, y_max=%g\n",pT,y_min,y_max); 
  double frac;

  for (int i=0; i<N_y; i++) {
    frac = (double)i/(double)(N_y-1);
    res = 0.;

    fprintf( out, "%.8e", y);

    for (int j=0;j<IRREPS;j++) {
      prob = reps[j](xi,&Fc);
      RpA_FCEL(Fc,pT,y,dsig,S,&R[j]);
      res+= prob*R[j];
      fprintf( out, "   %.8e", R[j]);
    }

    printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); 
    fprintf( out, "   %.8e\n", res );
    y += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

void R_scan_y(double pT) {
  int N_y;
  double R, Rmin, Rmax, y, y_min, y_max, step;

  char *prefix=(char*)"out/RpA_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{rs=%.2f,pT=%.1f}.dat",SQRTS,pT);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pA, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# columns: y, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  N_y=50; 
  y_min=-2.;
  y_max=+2.;
  // don't change anything after that.

  step=(y_max-y_min)/((double) N_y-1);
  y=y_min;

  printf(" Settings: pT=%g, with y_min=%g, y_max=%g\n",pT,y_min,y_max); 
  double frac;

  for (int i=0; i<N_y; i++) {
    frac = (double)i/(double)(N_y-1);

    R_limits(pT,y,&R,&Rmin,&Rmax);

    printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", y, R, Rmin, Rmax );
    y += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


void R_scan_pT(double y) {
  int N_pT;
  double R, Rmin, Rmax, pT, pT_min, pT_max, step;

  char *prefix=(char*)"out/RpA_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{rs=%.2f,y=%.1f}.dat",SQRTS,y);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pA, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# columns: pT, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  N_pT=50; 
  pT_min=1.;
  pT_max=10.;
  // don't change anything after that.

  step=(pT_max-pT_min)/((double) N_pT-1);
  pT=pT_min;

  printf(" Settings: y=%g, with pT_min=%g, pT_max=%g\n",y,pT_min,pT_max); 
  double frac;

  for (int i=0; i<N_pT; i++) {
    frac = (double)i/(double)(N_pT-1);

    R_limits(pT,y,&R,&Rmin,&Rmax);

    printf(" pT = %.5e , [%2.2f%]\n", y , 100.*frac); 
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", pT, R, Rmin, Rmax );
    pT += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

