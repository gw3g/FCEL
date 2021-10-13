#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
/*
 *    code for LHC predictions, RpA
 *    @author  GJ
 *    @version 0.8
 *
 *    COMPILE: gcc -lm -lgsl -lgslcblas main.c
 */


/*--------------------------------------------------------------------*/

#include "../fcel.h"
#include "sigma.h"
double SQRTS = 8.16;   // TeV
double A = 208.;       // Pb
double alpha_s = 0.5;  // coupling
const double mu_D = 1.8;     // GeV (D-meson mass)
const double mu_B = 5.3;     // GeV (B-meson mass)
double MU   = mu_D;

//FILE *in;
FILE *out;
void R_reps(double,char);   // pT input
void R_scan_y(double,char) ;// pT fixed
void R_scan_pT(double,char);// y  fixed
void R_scan_yintegrated(double,double,char);
void R_FB(double);     // y  fixed
int  progress = 0;
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
         x2     = ( mT/((.5)*rs) )*exp(-y), // xi = .5
         y_max  = log(rs/mT),
         //x_max  = 1.,//fmin(1.,(.5*exp(y_max-y)-1.)), // xi = .5
         x_max  = 1.,//fmin(1.,(fmin(xi,1.-xi)*exp(y_max-y)-1.)), // xi = .5
         kappa  = 4.*sqrt( SQR(pT) + SQR(MU) )/(rs), // meson transverse mass
         eps_g  = -3.;

  double as   = alpha_s/(2.*M_PI),
         M2   = mT*mT/( xi*(1.-xi) ),
         //M2   = (SQR(mT)+ xi*SQR(2.*m))/( xi*(1.-xi) ),
         Q2p  = Qs2(L_p,qhat,x2),
         //Q2A  = Qs2(L_eff(A),qhat,x2);
         Q2A  = Qs2(10.11,qhat,x2);

  double _integrand(double x, void *params_out) {
    double A  = ((double *)params_out)[0];
    double B  = ((double *)params_out)[1];
    double C  = ((double *)params_out)[2];

    double  dy = SGN(Fc)*log(1.+x),
            P  = phat(x,A,B,C);

    double in[3] = {n,kappa,eps_g}; // ~ (1.-xF)^n
    return P*sig(y+dy,in)/( sig(y,in)*pow(1.+x,SGN(Fc)) );
  };

  double out[3] = {Q2A/M2,Q2p/M2,fabs(Fc)*as};
  integrator(0.,x_max,_integrand,out,&res_outer,&err); // do integral

  *res = res_outer ;
}

double XS_pp(double pT, double y, double (*sig)(double,void*), void *params) {

  double z        = ((double *)params)[2]; // fragmentation var
  double n        = ((double *)params)[3]; // pp exponent

  double rs     = SQRTS*1e3, // root s
         kappa  = 4.*sqrt( SQR(pT) + SQR(MU) )/(rs); // meson transverse mass

    double in[3] = {n,kappa,-3.}; // ~ (1.-xF)^n
    return sig(y,in);
}

void RFB_FCEL(Fc,pT,ymin,ymax,sig,params,res)
  double Fc;                                // Ca + CR - Cb
  double pT;                                // transverse mom
  double ymin,ymax;
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

  double res_outer, res_inner, err;

  double _outer(double y, void *params_in) {

  double kT = pT/z; // partonic momenta rescaled by z
  double mT = sqrt( SQR(kT) + SQR(m) );

  double rs     = SQRTS*1e3, // root s
         x2     = ( mT/((.5)*rs) )*exp(-y), // xi = .5
         y_max  = log(rs/mT),
         x_max  = 1.,//fmin(1.,(.5*exp(y_max-y)-1.)), // xi = .5
         kappa  = 4.*sqrt( SQR(pT) + SQR(MU) )/(rs), // meson transverse mass
         eps_g  = -3.;

  double as   = alpha_s/(2.*M_PI),
         M2   = mT*mT/( xi*(1.-xi) ),
         Q2p  = Qs2(L_p,qhat,x2),
         Q2A  = Qs2(L_eff(A),qhat,x2);

  double _integrand(double x, void *params_out) {
    double A  = ((double *)params_out)[0];
    double B  = ((double *)params_out)[1];
    double C  = ((double *)params_out)[2];

    double  dy = SGN(Fc)*log(1.+x),
            P  = phat(x,A,B,C),
            denom;


    double in[3] = {n,kappa,eps_g}; // ~ (1.-xF)^n
    integrator(ymin,ymax,sig,in,&denom,&err);
    return P*sig(y+dy,in)/( denom*pow(1.+x,SGN(Fc)) );
  };

  double out[3] = {Q2A/M2,Q2p/M2,fabs(Fc)*as};
  integrator(0.,x_max,_integrand,out,&res_inner,&err); // do integral
  return res_inner;
  }

  integrator(ymin,ymax,_outer,NULL,&res_outer,&err); // do integral
  *res = res_outer ;
}


/*--------------------------------------------------------------------*/
// Hessian method
//
#define PARAMS 5  //  q0, xi, z, n, m_Q
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

void R_integrated(double pT, double y_min, double y_max, 
                  double *R_ave, double *R_min, double *R_max) {

  double Sp[PARAMS], Sm[PARAMS], res1=0.,res2=0.,
         numerator_Sp,    numerator_Sm,    numerator_S,
         denominator_Sp,  denominator_Sm,  denominator_S,
         RS, RSp, RSm, err;// = R_sum_(pT,y,S)*XS_pp(pT,y,dsig,S);

  double _pA(double y, void *params) {
    return R_sum_(pT,y,params)*XS_pp(pT,y,dsig,params);
  };

  double _pp(double y, void *params) {
    return XS_pp(pT,y,dsig,params);
  };

  integrator(y_min,y_max,_pA,S,&numerator_S,&err);
  integrator(y_min,y_max,_pp,S,&denominator_S,&err);
  RS = numerator_S/denominator_S;

  memcpy(Sp,S,sizeof(double)*PARAMS);
  memcpy(Sm,S,sizeof(double)*PARAMS);

  for (int k=0; k<PARAMS; k++) {

    // modify running arrays ...
    Sp[k]+=dS[k]; Sm[k]-=dS[k];

    integrator(y_min,y_max,_pA,Sp,&numerator_Sp,&err);
    integrator(y_min,y_max,_pp,Sp,&denominator_Sp,&err);
    RSp = numerator_Sp/denominator_Sp;

    integrator(y_min,y_max,_pA,Sm,&numerator_Sm,&err);
    integrator(y_min,y_max,_pp,Sm,&denominator_Sm,&err);
    RSm = numerator_Sm/denominator_Sm;

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
#include "scanners.h"

int main() {
  double Fc; char meson;
  //gsl_set_error_handler_off(); // live on the edge

  //R_Casimir(0.,+3.5);
  //R_Casimir(0.,+5.5);
  //R_Casimir(0.,-5.5);
  //R_Casimir(2.,+3.5);
  //R_Casimir(2.,-3.5);

/*
  // Repeat figs for 2003.06337

  SQRTS = 8.16; // TeV
  A = 208.;     // Pb
  alpha_s = 0.5;// coupling
  MU = 0.;
  meson = 'h';

  S[2] = .6; dS[2] = .2; // z -- frag. variable
  S[3] = 6.; dS[3] = 2.; // n -- exponent
  S[4] = 0.; dS[4] = 0.; // eff. quark mass

  { // FIG 2,3
    channel(1); // g, g -> g, g
    R_reps(2.);
    R_scan_y(2.,meson);   R_scan_y(6.,meson);
    R_scan_pT(0.,meson);  R_scan_pT(3.,meson);    R_scan_pT(5.,meson);
  }

  SQRTS = 5.02; // TeV
  { // FIG 5
    channel(1); // g, g -> g, g
    R_scan_pT(0.,meson);
  }

  SQRTS = 8.16; // TeV
  { // FIG 6
    channel(2); // q, g -> q, g
    R_reps(2.);
    R_scan_y(2.,meson);   R_scan_y(6.,meson);
  }

  { // FIG 7
    channel(3); // g, q -> q, g
    R_reps(2.);
    R_scan_y(2.,meson);   R_scan_y(6.,meson);
  }

  { // FIG 8
    channel(4); // g, g -> q, q
    R_reps(2.);
    R_scan_y(2.,meson);   R_scan_y(6.,meson);
  } //*/

  // >> LCHb, 5 TeV D0 <<

  //SQRTS = 5.02; // TeV
  //SQRTS = 8.16; // TeV
  //A = 208.;     // Pb
  //alpha_s = 0.5;// coupling

  // OXYGEN
  SQRTS = 9.9; // TeV
  A = 16.;     // Pb
  alpha_s = 0.5;// coupling

  S[2] = 0.8;  dS[2] = 0.2;// z -- frag. variable
  S[3] = 4.0;  dS[3] = 1.0;// n -- exponent
  S[4] = 1.3;  dS[4] = 0.2;// eff. quark mass
  MU = mu_D;
  meson = 'D';

  //S[2] = 0.7;  dS[2] = 0.2;// z -- frag. variable
  //S[3] = 15.0;  dS[3] = 5.0;// n -- exponent
  //S[4] = 0.0;  dS[4] = 0.0;// eff. quark mass
  //MU = 0;
  //meson = 'h';


  //S[2] = 0.9;  dS[2] = 0.1;// z -- frag. variable
  //S[3] = 2.0;  dS[3] = 0.5;// n -- exponent
  //S[4] = 4.6;  dS[4] = 0.5;// eff. quark mass
  //MU = mu_B;
  //meson = 'B';

  {
    //channel(1); // g, g -> g, g
    channel(4); // g, g -> q, q

    R_scan_y(1.,meson);
    R_scan_y(2.,meson);
    R_scan_y(3.,meson);
    R_scan_y(4.,meson);
    R_scan_y(5.,meson);
    //R_scan_pT(3,meson);
    //R_scan_pT(-4,meson);
    //R_scan_pT(-2,meson);
    //R_scan_pT(+0,meson);
    //R_scan_pT(+0,meson);
    //R_scan_pT(+1,meson);
    //R_scan_pT(+2,meson);
    //R_scan_pT(+3,meson);
    //R_scan_pT(+4,meson);
    //R_scan_pT(+5,meson);
    //R_scan_pT(+6,meson);
    //R_scan_y(0.,meson);
    //R_scan_y(1.,meson);
    //R_scan_y(2.,meson);
    //R_scan_yintegrated(2.5,3.5,meson);
  }//*/


  /*{ //
    channel(1); // g, g -> g, G

    R_reps(.5);  
    R_reps(1.5);  
    R_reps(2.5);  
    R_reps(3.5);  
    R_reps(4.5);  
    R_reps(5.5);  
    R_reps(6.5);  
    R_reps(7.5);  
    R_reps(8.5);  
    R_reps(9.5);
  }//*/


  /*{ //
    channel(4); // g, g -> q, q

    R_scan_pT(+2.25,meson);  
    R_scan_pT(+2.75,meson);  
    R_scan_pT(+3.25,meson);  
    R_scan_pT(+3.75,meson);  

    R_scan_pT(-2.75,meson);  
    R_scan_pT(-3.25,meson);  
    R_scan_pT(-3.75,meson);  
    R_scan_pT(-4.25,meson);  

    R_scan_y(.5,meson);  
    R_scan_y(1.5,meson);  
    R_scan_y(2.5,meson);  
    R_scan_y(3.5,meson);  
    R_scan_y(4.5,meson);  
    R_scan_y(5.5,meson);  
    R_scan_y(6.5,meson);  
    R_scan_y(7.5,meson);  
    R_scan_y(8.5,meson);  
    R_scan_y(9.5,meson);
  }//*/

  /*{
    SQRTS = 8.16; // TeV
    A = 208.;     // Pb
    alpha_s = 0.5;// coupling

    S[2] = 0.8;  dS[2] = 0.2;// z -- frag. variable
    S[3] = 4.0;  dS[3] = 1.0;// n -- exponent
    S[4] = 1.5;  dS[4] = 0.5;// eff. quark mass

    channel(4); // g, g -> q, q

    R_FB(1.);
  }//*/

  /*{
    SQRTS = 8.16; // TeV
    A = 208.;     // Pb
    alpha_s = 0.5;// coupling

    S[2] = 0.9;  dS[2] = 0.1;// z -- frag. variable
    S[3] = 4.0;  dS[3] = 1.0;// n -- exponent
    //S[4] = 1.3;  dS[4] = 0.2;// eff. quark mass
    S[4] = 4.6;  dS[4] = 0.5;// eff. quark mass

    channel(4); // g, g -> q, q

    R_FB(1.);
  }//*/


   //R_scan_y(10.);

   //R_FB(2.5);
   //R_FB(4.0);

   //R_scan_pT(0.);
   //R_scan_pT(3.);
   //R_scan_pT(5.);
   return 0;
}
