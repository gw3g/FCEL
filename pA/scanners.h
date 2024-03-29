
/*--------------------------------------------------------------------*/
void R_reps(double pT, char H) {
  int N_y;
  double y, y_min, y_max, step;

  double xi=S[1], R[3];
  double prob;
  double res;
  double Fc;

  char *prefix=(char*)"out/R_reps_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{rs=%.2f,pT=%.1f}_%c.dat",SQRTS,pT,H);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pPb, z=%g, alpha=%g\n",S[2],alpha_s);
  fprintf(out,"# columns: y, {R1,R2,...}, R_ave\n");

  // Here are some parameters that can be changed:
  N_y=200;
  y_min=-6.;
  y_max=+6.;
  // don't change anything after that.

  step=(y_max-y_min)/((double) N_y-1);
  y=y_min;

  if (progress) { printf(" Settings: pT=%g, with y_min=%g, y_max=%g\n",pT,y_min,y_max); }
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

    if (progress) { printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); }
    fprintf( out, "   %.8e\n", res );
    y += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);

}

void R_scan_full(char H) {
  int N_y, N_pT;
  double R, Rmin, Rmax, 
         y, y_min, y_max, 
         pT, pT_min, pT_max,
         step_y, step_pT, dummy;

  //char *prefix=(char*)"out/RpA_";
  char *prefix=(char*)"out/RoR_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  //sprintf(suffix,"_{sqrs=%.2f}_%c.dat",SQRTS,H);
  sprintf(suffix,"_%c.dat",H);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  //fprintf(out,"# R_pA, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# R-over-R = R_pA(8.16TeV)/RpA(5.02TeV), A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# Columns:\n");
  fprintf(out,"# pT/GeV, y, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  N_pT    = 101; 
  pT_min  = 0.;
  pT_max  = 10.;

  N_y     = 241; 
  y_min   = -6.;
  y_max   = +6.;

  //N_y=21; 
  //y_min=-5.;
  //y_max=+5.;
  // don't change anything after that.

  step_y =(y_max-y_min)/((double) N_y-1);
  step_pT=(pT_max-pT_min)/((double) N_pT-1);

  y  = y_min;
  pT = pT_min;

  if (progress) { printf(" Settings: sqrs=%g, with y_min=%g,  y_max=%g\n",SQRTS,y_min,y_max); }
  if (progress) { printf("                         pT_min=%g, pT_max=%g\n\n",pT_min,pT_max); }
  double frac;

  for (int j=0; j<N_pT; j++) {
    for (int i=0; i<N_y; i++) {
      frac = ((double)j*(double)N_y+(double)i)/(double)(N_y*N_pT-1);
      if (progress) { printf(" pT = %.5e GeV, y = %.5e , [%2.2f%]\n", pT, y, 100.*frac); }

      //R_limits(pT,y,&R,&Rmin,&Rmax); // calculation
      double_R(pT,y,&R,&Rmin,&Rmax); // calculation

      fprintf( out, "%.8e   %.8e   %.8e   %.8e   %.8e\n", pT, y, R, Rmin, Rmax);

      y += step_y;
    }
    y = y_min;
    pT += step_pT;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


void R_scan_y(double pT, char H) {
  int N_y;
  double R, Rmin, Rmax, y, y_min, y_max, step, dummy;

  char *prefix=(char*)"out/RpO_";
  //char *prefix=(char*)"out/RpA_GJ_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{rs=%.2f,pT=%.1f}_%c.dat",SQRTS,pT,H);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pA, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# columns: y, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  N_y=200; 
  y_min=-6.;
  y_max=+6.;
  //N_y=21; 
  //y_min=-5.;
  //y_max=+5.;
  // don't change anything after that.

  step=(y_max-y_min)/((double) N_y-1);
  y=y_min;

  if (progress) { printf(" Settings: pT=%g, with y_min=%g, y_max=%g\n",pT,y_min,y_max); }
  double frac;

  for (int i=0; i<N_y; i++) {
    frac = (double)i/(double)(N_y-1);

    R_limits(pT,y,&R,&Rmin,&Rmax); // calculation

    if (progress) { printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); }
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", y, R, Rmin, Rmax);

    y += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


void R_scan_pT(double y, char H) {
  int N_pT;
  double R, Rmin, Rmax, pT, pT_min, pT_max, step, dummy;

  char *prefix=(char*)"out/RpA_GJ_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{rs=%.2f,y=%.1f}_%c.dat",SQRTS,y,H);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pA, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# columns: pT/GeV, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  //N_pT=100; 
  //pT_min=.1;
  //pT_max=10.;
  N_pT=21; 
  pT_min=.0;
  pT_max=10.;
  // don't change anything after that.

  step=(pT_max-pT_min)/((double) N_pT-1);
  pT=pT_min;

  if (progress) { printf(" Settings: y=%g, with pT_min=%g, pT_max=%g\n",y,pT_min,pT_max); }
  double frac;

  for (int i=0; i<N_pT; i++) {
    frac = (double)i/(double)(N_pT-1);

    R_limits(pT,y,&R,&Rmin,&Rmax); // calculation

    if (progress) { printf(" pT = %.5e , [%2.2f%]\n", y , 100.*frac); }
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", pT, R, Rmin, Rmax );

    pT += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

void R_scan_yintegrated(double y_min, double y_max, char H) {
  int N_pT;
  double R, Rmin, Rmax, pT, pT_min, pT_max, step, dummy;

  //char *prefix=(char*)"out/RpA_integrated_";
  char *prefix=(char*)"out/RFB_integrated_";
  char  suffix[70];
  char  filename[100];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{rs=%.2f,ymin=%.1f,ymax=%.1f}_%c.dat",SQRTS,y_min,y_max,H);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_FB, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# Columns:\n");
  fprintf(out,"# pT/GeV, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  //N_pT=100; 
  //pT_min=.1;
  //pT_max=20.;
  N_pT=101; 
  pT_min=.0;
  pT_max=10.;
  // don't change anything after that.

  step=(pT_max-pT_min)/((double) N_pT-1);
  pT=pT_min;

  if (progress) { printf(" Settings: y_min=%g, y_max=%g, with pT_min=%g, pT_max=%g\n",
                          y_min,y_max,pT_min,pT_max); }
  double frac;

  for (int i=0; i<N_pT; i++) {
    frac = (double)i/(double)(N_pT-1);

    R_integrated(pT,y_min,y_max,&R,&Rmin,&Rmax); // calculation

    if (progress) { printf(" pT = %.5e , [%2.2f%]\n", pT , 100.*frac); }
    fprintf( out, "%.8e   %.8e   %.8e   %.8e\n", pT, R, Rmin, Rmax );

    pT += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


void R_Casimir(double pT, double y) { // function of global final colour
  int N_CR;
  double  R, R_minus, R_plus, Fc,
         CR, CR_min, CR_max, step;

  char *prefix=(char*)"out/R_CR_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"{rs=%.2f,pT=%.1f,y=%.1f}.dat",SQRTS,pT,y);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# columns: C_R, R_ave, R_min, R_max\n");

  // Here are some parameters that can be changed:
  N_CR=80; 
  CR_min=-20.;
  CR_max=20.;
  // don't change anything after that.

  step=(CR_max-CR_min)/((double) N_CR-1);
  CR=CR_min;

  if (progress) { printf(" Settings: y=%g, pT=%g, with CR_min=%g, CR_max=%g\n",y,pT,CR_min,CR_max); }
  double frac;

  for (int i=0; i<N_CR; i++) {
    frac = (double)i/(double)(N_CR-1);
    Fc = Cf + CR - Cf;

    RpA_FCEL(Fc,pT,y-.5,dsig,S,&R_minus);
    RpA_FCEL(Fc,pT,y,dsig,S,&R);
    RpA_FCEL(Fc,pT,y+.5,dsig,S,&R_plus);

    if (progress) { printf(" pT = %.5e, y = %.5e , [%2.2f%]\n", pT, y , 100.*frac); }
    fprintf( out, "%.8e   %.8e   %.8e    %.8e\n", CR, R, R_minus, R_plus );
    CR += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}

/*--------------------------------------------------------------------*/

void R_rel_uncertainty(double pT, char H) {
  int N_y;
  double R, Rp, Rm, y, y_min, y_max, step, dummy;

  //double dS[5] = {.02,.25,.2,1.,0.2};
  //double  S[5] = {.07,.50,.8,4.,1.3};
  double dS[5] = {.02,.25,.2,5.,0.0};
  double  S[5] = {.07,.50,.6,15.,0.0};

  char *prefix=(char*)"out/RelError_pO_";
  //char *prefix=(char*)"out/RpA_GJ_";
  char  suffix[20];
  char  filename[50];

  // filename
  strcpy(filename,prefix);
  sprintf(suffix,"_{rs=%.2f,pT=%.1f}_%c.dat",SQRTS,pT,H);
  strcat(filename,reaction);
  strcat(filename,suffix);
  out=fopen(filename,"w");
  fprintf(out,"# R_pA, A=%.1f, alpha=%g\n",A,alpha_s);
  fprintf(out,"# columns: y, R, R(qhat-), R(qhat+), R(xi-), R(xi+), R(z-), R(z+), R(n-), R(n+), R(m-), R(m+) \n");

  // Here are some parameters that can be changed:
  N_y=200; 
  y_min=-6.;
  y_max=+6.;
  //N_y=21; 
  //y_min=-5.;
  //y_max=+5.;
  // don't change anything after that.

  step=(y_max-y_min)/((double) N_y-1);
  y=y_min;

  if (progress) { printf(" Settings: pT=%g, with y_min=%g, y_max=%g\n",pT,y_min,y_max); }
  double frac;

  for (int i=0; i<N_y; i++) {
    frac = (double)i/(double)(N_y-1);
    fprintf( out, "%.8e   ", y);

    R = R_sum_(pT,y,S); 
    fprintf( out, "%.8e", R);

    for (int j=0; j<5; j++) {
      S[j] += dS[j];
      Rp = R_sum_(pT,y,S); 
      S[j] -= 2.*dS[j];
      Rm = R_sum_(pT,y,S); 
      S[j] += dS[j];
      fprintf( out, "   %.8e   %.8e", Rm, Rp);
    }
    fprintf( out, "\n", Rm, Rp);

    if (progress) { printf(" y = %.5e , [%2.2f%]\n", y , 100.*frac); }

    y += step;
  }

  printf(" Saved to file ["); printf(filename); printf("]\n"); fclose(out);
}


