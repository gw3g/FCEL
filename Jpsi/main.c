#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
/*
 *    code for LHC predictions, RpA
 *    @author  GJ
 *    @version 0.?
 *
 *    COMPILE: gcc -lm -lgsl -lgslcblas main.c
 */


/*--------------------------------------------------------------------*/

#include "../fcel.h"
#include "sigma.h"
double SQRTS = .0189;   // TeV
double A = 184.;       // Pb
double alpha_s = 0.5;  // coupling

//FILE *in;
FILE *out;
void R_reps(double,char);   // pT input
void R_scan_full(char);
void R_scan_y(double,char) ;// pT fixed
void R_scan_pT(double,char);// y  fixed
void R_scan_yintegrated(double,double,char);
void R_FB(double);     // y  fixed
int  progress = 1;
void R_Casimir(double,double);// (pT,y)

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

void RpA_FCEL(Fc,xF,sig,params,res)
  double Fc;                                // Ca + CR - Cb
  double xF;                                 // rapidity
  double (*sig)(double,void*);              // l2_perp(Lp,q0,x2)
  void *params;                             // {u,v,w,x,y,z} = (see below)
  double *res;
{
  /*
   *  nuclear mod factor due to FCEL
   */
  if (fabs(Fc)<tol) {*res=1.; return;} // no mods.
  //printf("TEST - 1\n"); 

  double qhat     = ((double *)params)[0]; // transport param.
  double pT       = ((double *)params)[1]; // transverse mom.
  double nn       = ((double *)params)[2]; // pp exponent
  double M        = ((double *)params)[3]; // quarkonia mass (?)

  double res_outer, err;

  double MT = sqrt( SQR(pT) + SQR(M) ), M2   = MT*MT;
  //printf("           pT = %.2f, n = %.2f, M = %.2f\n",pT,nn,M); 

  double rs     = SQRTS*1e3, // root s [GeV]
         s      = SQR(rs),
         Ep     = s/(2.*mp),
         E      = Ep*( .5*xF + sqrt( SQR(.5*xF) + M2/s ) ),
         xp     = sqrt( SQR(xF) + 4.*M2/s ),
         x1     = .5*(xp+xF),
         x2     = .5*(xp-xF),
         x_max  = fmin(Ep/E-1.,1.);//fmin(1.,(fmin(xi,1.-xi)*exp(y_max-y)-1.)), // xi = .5

  double as   = alpha_s/(2.*M_PI),
         Q2p  = Qs2(L_p,qhat,x2),
         Q2A  = Qs2(L_eff(A),qhat,x2);
         //Q2A  = Qs2(10.11,qhat,x2);
         //Q2A  = Qs2(3.24,qhat,x2); // Be
         //Q2A  = Qs2(9.35,qhat,x2); // W
         //
  //printf("TEST - 2\n"); 

  if (xF>1.-SQR(M)/s) { xF=1.-SQR(M)/s; }
  double _integrand(double x, void *params_out) {
    double A_  = ((double *)params_out)[0];
    double B_  = ((double *)params_out)[1];
    double C_  = ((double *)params_out)[2];
    //printf("TEST - 3b\n"); 

    double P  = phat(x,A_,B_,C_);
    //printf("TEST - 3c\n"); 

    double xFp = E*(1.+x)/Ep;
           xFp = xFp - M2/s/xFp;
    if (xFp>1.-SQR(M)/s) { xFp=1.-SQR(M)/s; }
    //printf("TEST - 3d: xF = %.2f, xF' = %.2f\n",xF,xFp); 
    //printf("           M2 = %.2f, s = %.2f\n",M2,s); 

    double in[2] = {nn,M2/s}; //
    if (xp>1.) { return 0.; };
    return P*sig(xFp,in)/( sig(xF,in)*pow(1.+x,SGN(Fc)) );
  };

  //double out[3] = {Q2A/M2,Q2p/M2,fabs(Fc)*as};
  double out[3] = {Q2A/M2,fmax(Q2p,SQR(.25))/M2,fabs(Fc)*as};
  integrator(0.,x_max,_integrand,out,&res_outer,&err); // do integral

  //printf("TEST - 4\n"); 
  *res = res_outer ;
}

/*
double XS_pp(double pT, double y, double (*sig)(double,void*), void *params) {

  double z        = ((double *)params)[2]; // fragmentation var
  double n        = ((double *)params)[3]; // pp exponent

  double rs     = SQRTS*1e3, // root s
         kappa  = 4.*sqrt( SQR(pT) + SQR(MU) )/(rs); // meson transverse mass

    double in[3] = {n,kappa,-3.}; // ~ (1.-xF)^n
    return sig(y,in);
}//*/


/*--------------------------------------------------------------------*/
// Hessian method
//
/*#define PARAMS 5  //  q0, xi, z, n, m_Q
double dS[PARAMS] = {.02,.25,.2,1.,0.1};
double  S[PARAMS] = {.07,.50,.8,4.,1.5};

double R_sum_(double pT, double y, void *params) {
  double Fc, temp, prob, res=.0,
         xi  = ((double *)params)[1];

  for (int i=0;i<IRREPS;i++) {
    prob = reps[i](xi,&Fc);
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

void double_R(double pT, double y, double *R_ave, double *R_min, double *R_max) {
  double Sp[PARAMS], Sm[PARAMS], res1=0.,res2=0.,
         RSp,        RSm,        RS;// = R_sum_(pT,y,S);

  // apply Hessian method to the double ratio RpA(\sqrt{s}=8.16TeV)/RpA(\sqrt{s}=5.02TeV)
  SQRTS = 8.16;
  RS = R_sum_(pT,y,S);
  SQRTS = 5.02;
  RS = RS/R_sum_(pT,y,S);

  memcpy(Sp,S,sizeof(double)*PARAMS);
  memcpy(Sm,S,sizeof(double)*PARAMS);

  for (int k=0; k<PARAMS; k++) {

    // modify running arrays ...
    Sp[k]+=dS[k]; Sm[k]-=dS[k];

    SQRTS = 8.16;
    RSp = R_sum_(pT,y,Sp); 
    RSm = R_sum_(pT,y,Sm);
    SQRTS = 5.02;
    RSp = RSp/R_sum_(pT,y,Sp);
    RSm = RSm/R_sum_(pT,y,Sm);

    res1+= SQR( fmax( fmax( RSp-RS, RSm-RS ), 0. ) );
    res2+= SQR( fmax( fmax( RS-RSp, RS-RSm ), 0. ) );

    // ... undo those mods.
    Sp[k]-=dS[k]; Sm[k]+=dS[k];
  };

  *R_ave = RS;
  *R_max = RS + sqrt(res1);
  *R_min = RS - sqrt(res2);
}

void R_integrated(double pT, double y_min, double y_max, 
                  double *R_ave, double *R_min, double *R_max) {

  double Sp[PARAMS], Sm[PARAMS], res1=0.,res2=0.,
         numerator_Sp,    numerator_Sm,    numerator_S,
         denominator_Sp,  denominator_Sm,  denominator_S,
         RS, RSp, RSm, err;// = R_sum_(pT,y,S)*XS_pp(pT,y,dsig,S);
                           //
  memcpy(Sp,S,sizeof(double)*PARAMS);
  memcpy(Sm,S,sizeof(double)*PARAMS);

  double _pA(double y, void *params) {
    return R_sum_(pT,y,params)*XS_pp(pT,y,dsig,params);
  };

  double _pp(double y, void *params) {
    return XS_pp(pT,y,dsig,params);
  };

  // Forward:
  integrator(y_min,y_max,_pA,S,&numerator_S,&err);
  integrator(y_min,y_max,_pp,S,&denominator_S,&err);
  RS = numerator_S/denominator_S;

  // Backward:
  integrator(-y_max,-y_min,_pA,S,&numerator_S,&err);
  integrator(-y_max,-y_min,_pp,S,&denominator_S,&err);
  RS = RS/(numerator_S/denominator_S);


  for (int k=0; k<PARAMS; k++) {

    // modify running arrays ...
    Sp[k]+=dS[k]; Sm[k]-=dS[k];

    // Forward:
    integrator(y_min,y_max,_pA,Sp,&numerator_Sp,&err);
    integrator(y_min,y_max,_pp,Sp,&denominator_Sp,&err);
    RSp = numerator_Sp/denominator_Sp;
    // Backward:
    integrator(-y_max,-y_min,_pA,Sp,&numerator_Sp,&err);
    integrator(-y_max,-y_min,_pp,Sp,&denominator_Sp,&err);
    RSp = RSp/(numerator_Sp/denominator_Sp);

    // Forward:
    integrator(y_min,y_max,_pA,Sm,&numerator_Sm,&err);
    integrator(y_min,y_max,_pp,Sm,&denominator_Sm,&err);
    RSm = numerator_Sm/denominator_Sm;
    // Forward:
    integrator(-y_max,-y_min,_pA,Sm,&numerator_Sm,&err);
    integrator(-y_max,-y_min,_pp,Sm,&denominator_Sm,&err);
    RSm = RSm/(numerator_Sm/denominator_Sm);

    res1+= SQR( fmax( fmax( RSp-RS, RSm-RS ), 0. ) );
    res2+= SQR( fmax( fmax( RS-RSp, RS-RSm ), 0. ) );

    // ... undo those mods.
    Sp[k]-=dS[k]; Sm[k]+=dS[k];
  };

  *R_ave = RS;
  *R_max = RS + sqrt(res1);
  *R_min = RS - sqrt(res2);
}
//*/

/*--------------------------------------------------------------------*/

void R_scan_xF(double A1, double A2, double nn, double p_T, char *tag) {
  int N_xF;
  double R1, R2, R, xF, xF_min, xF_max, step, dummy;

  char *prefix=(char*)"out/R_";
  //char *prefix=(char*)"out/RpA_GJ_";
  char  suffix[20];
  char  filename[50];

  double pT = 1.;
  // filename
  strcpy(filename,prefix);
  strcat(filename,tag);
  sprintf(suffix,"_{rs=%.1f,pt=%.1f}.dat",1e3*SQRTS,p_T);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_AB, A=%.1f, B=%.1f, alpha=%g\n",A1,A2,alpha_s);
  fprintf(out,"# columns: xF, RpA, RpB, R_AB\n");

  // Here are some parameters that can be changed:
  N_xF=220; 
  xF_min=-.2;
  xF_max=.9;
  //N_y=21; 
  //y_min=-5.;
  //y_max=+5.;
  // don't change anything after that.

  step=(xF_max-xF_min)/((double) N_xF-1);
  xF=xF_min;

  if (progress) { printf(" Settings: pT=%g, with xF_min=%g, xF_max=%g\n",pT,xF_min,xF_max); }
  double frac;
  double par[4] = {.075,p_T,nn,3.};

  for (int i=0; i<N_xF; i++) {
    frac = (double)i/(double)(N_xF-1);

    A = A1;
    RpA_FCEL(Nc,xF,dsig,par,&R1); // calculation
    A = A2;
    RpA_FCEL(Nc,xF,dsig,par,&R2); // calculation

    if (progress) { printf(" xF = %.5e , [%2.2f%]\n", xF , 100.*frac); }
    //fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", y, R, Rmin, Rmax);
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", xF, R1, R2, R1/R2);

    xF += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


int main() {
  double Fc; char meson;
  //gsl_set_error_handler_off(); // live on the edge
  //
  SQRTS = .0387;
  //A = 184.;
  R_scan_xF(184.,9.,4.5,1.,"WBe");
  R_scan_xF(56.,9.,4.5,1.,"FeBe");

  SQRTS = .0189;
  R_scan_xF(184.,27.,1.4,0.,"WAl");
  R_scan_xF(184.,27.,1.4,1.,"WAl");
  R_scan_xF(184.,27.,1.4,2.,"WAl");
  R_scan_xF(184.,27.,1.4,3.,"WAl");
  R_scan_xF(184.,27.,1.4,4.,"WAl");
  //A = 9.;
  //R_scan_xF();
   //R_scan_pT(3.);
   //R_scan_pT(5.);
   return 0;
}
