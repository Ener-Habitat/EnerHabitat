/**************** conduction.c para no aire acondicionado 1D  *********/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "arrays.h"             //includes the routine to make dynamic arrays
#include "new_input.h"              //includes the routines to start the executabla
#include "tools.h"
#include "sol.h"
/******************************************************************************************
*******************************************************************************************/


static char help[]="\
***************************************************************************\n\
conduction.c\n\
***************************************************************************\n\
\n\
One dimensional transient conduction in line\n\
***************************************************************************\n\
 ";

main(int argc,char *argv[]) 
{
  
  /* */
  /* */
   
  
  char inpfile[80];
  /* Declare variables defined by the user */
  double L1,L2,L3,L4,L5,L6,L7,La;
  double dt,ho,hi;
  double k1,k2,k3,k4,k5,k6,k7;
  double rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7;
  double t_Tamax,Tmax,Tmin,Ig,Id,Ib,beta,gamma,A;
  double Lat,Lon;
  double Fc;
  double a11,a12,a13;
  double a21,a22,a23;
  double e21,e22,e23;
  int miniDeltaTK;
  int month,nx,cont_herramienta;
  char user[15];
  
  if (argv[1] == NULL || argv[1][0] == '-') 
    {
      strcpy(inpfile, "muro.e");
      strcat(inpfile, ".inp");
    } 
  else 
    strcpy(inpfile, argv[1]);
  input_reset();
  input_insert("N", "Number of elements", &nx, 'i');
  input_insert("mDTK", "Mini DeltaTK", &miniDeltaTK, 'i');
  input_insert("tT", "Hour of maximum temperature ", &t_Tamax, 'd');
  input_insert("Tmax", "Maximum temperature", &Tmax, 'd');
  input_insert("Tmin", "Minimum temperature", &Tmin, 'd');
  input_insert("Ig", "Global solar radiation", &Ig, 'd');
  input_insert("Id", "Diffuse solar radiation", &Id, 'd');
  input_insert("Ib", "Beam solar radiation", &Ib, 'd');
  input_insert("beta", "Wall inclination [degree]", &beta, 'd');
  input_insert("gamma", "Wall orientation [degree]", &gamma, 'd');
  input_insert("a", "Absortance ", &A, 'd');
  input_insert("La", "Room lenght", &La, 'd');
  input_insert("Fc", "Correction factor(0,1]", &Fc, 'd');
  input_insert("a11", "thermal conductivity 1", &a11, 'd');
  input_insert("a12", "thermal conductivity 2", &a12, 'd');
  input_insert("a13", "thermal conductivity 2", &a13, 'd');
  input_insert("a21", "thermal conductivity 3", &a21, 'd');
  input_insert("a22", "thermal conductivity 3", &a22, 'd');
  input_insert("e21", "thermal conductivity 1", &e21, 'd');
  input_insert("e22", "thermal conductivity 2", &e22, 'd');
  input_insert("e23", "thermal conductivity 3", &e23, 'd');
  input_insert("k1", "thermal conductivity 1", &k1, 'd');
  input_insert("k2", "thermal conductivity 2", &k2, 'd');
  input_insert("k3", "thermal conductivity 3", &k3, 'd');
  input_insert("k4", "thermal conductivity 4", &k4, 'd');
  input_insert("rhoc1", "thermal capacity 1", &rhoc1, 'd');
  input_insert("rhoc2", "thermal capacity 2", &rhoc2, 'd');
  input_insert("rhoc3", "thermal capacity 3", &rhoc3, 'd');
  input_insert("rhoc4", "thermal capacity 3", &rhoc4, 'd');
  input_insert("L1", "L1 [m]", &L1, 'd');
  input_insert("L2", "L2 [m]", &L2, 'd');
  input_insert("L3", "L3 [m]", &L3, 'd');
  input_insert("L4", "L4 [m]", &L4, 'd');
  input_insert("ho", "h outside [W/m2 oC]", &ho, 'd');
  input_insert("hi", "h inside [W/m2 oC]", &hi, 'd');
  input_insert("user", "usuario", &user, 's');
  input_insert("mes", "mes", &month, 'i');
  input_insert("Lo", "Longitud", &Lon, 'd');
  input_insert("Lat", "Latitud", &Lat, 'd');
  input_insert("cont", "contador", &cont_herramienta, 'i');
  input_load(inpfile);
  if (input_options(argc,argv,help)) input();  
  input_save(inpfile);
  
  
  double dx,t,criterio,t_max,Tsa;
  int error,dias;
  double *a,*b,*c,*d,*T,*Tn,*To,*Tf,*P,*Q;
  double *rhoc,*k;
  int i,j,itera;
  double Hi;
  double Qin,rhoair,cair,Tint,Ta,Tc,DtaT;
  double Tintmax,Tintmin,Tsamax,Tsamin,Tintaverage;
  double Tsa1,Tsa0;
  double TPIcold,TPIhot,NumCold,NumHot,DenCold,DenHot;
  double DDHcold,DDHhot;
  double t_Tintmax,t_Tsamax;
  FILE *f_in;

  int *juliano;
  double *delta,*orto, *ocaso,*Ho, tau[13];
  
  L5 = L6 = L7 = 0.;
  rhoc5 = rhoc6 = rhoc7 = 0.;
  k5 = k6 = k7 = 0.;

  

  dt = 1.;
  t_max  = 86400.; 
  rhoair = 1.205; 
  cair = 1005.;  



  double L11,L22,L33;  //geometr'ia de la cavidad hueca
  double Ceq,C2,RR,Hr,Hc,Aa,Ac;
  
  a12 = a12/2.;
  L11 = L1 = e21;
  L22 = L2 = e22;
  L33 = L3 = e23;
  Aa = a21/(a11 + a21 + a12 );
  Ac = (a11 + a12)/(a11 + a21 + a12);
  Ceq = L1*rhoc1 + L2*(Aa*rhoair*cair+Ac*rhoc1)+L3*rhoc1;
  C2 = (Ceq - (L11 + L33)*rhoc1)/L22;
  rhoc2 = C2;
  dx = (L1+L2+L3+L4+L5+L6+L7)/nx;
  printf("Aa = %f\n",Aa);
  printf("Ra = %f\n",a21/e22); 
  

  printf("L1 = %f\n",e21);
  printf("L2 = %f\n",e22);
  printf("L3 = %f\n",e23);
  printf("rhoc1 = %f\n",rhoc1*e21);
  printf("rhoc2 = %f\n",rhoc1*e22*.21 + rhoair*cair*e22*.79);
  printf("rhoc3 = %f\n",rhoc1*e23);
  printf("rhoca = %f\n",rhoair*cair);
  printf("rhoc_total = %f\n",rhoc1*e21 + rhoc1*.21*e22 + rhoair*cair*.79*e22 + rhoc1*e23);
  printf("Ceq = %f\n",Ceq);//getchar();
  
  T  = one_d_double_array(nx);
  Tn = one_d_double_array(nx);
  To = one_d_double_array(nx);
  Tf = one_d_double_array(nx);
  k = one_d_double_array(nx);
  rhoc = one_d_double_array(nx);
  a = one_d_double_array(nx);
  b = one_d_double_array(nx);
  c = one_d_double_array(nx);
  d = one_d_double_array(nx);
  P = one_d_double_array(nx);
  Q = one_d_double_array(nx);

  juliano = one_d_int_array(13);
  delta = one_d_double_array(13);
  orto = one_d_double_array(13);
  ocaso = one_d_double_array(13);
  Ho = one_d_double_array(13);


  double Econd,Econv,Erad;


  double Tintmaxo,Tintmino,Tintmaxf,Tintminf;
  dia_juliano(juliano,15,2011,12);   //dia, anio, mes
  calculo_declinacion_delta(delta,juliano,12);
  calculo_orto_ocaso(orto, ocaso, Lat, delta, 12);
  calculo_duracion_dia(tau, 12, orto);
  calculo_hora_orto(Ho, tau, 12);
  set_krhoc(k,rhoc,dx,L1,L2,L3,L4,L5,L6,L7,k1,k2,k3,k4,k5,k6,k7,rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7);
  Ta_Tc_DtaT(Ig,ho,Tmax,Tmin,t_Tamax,Ho[month],&Ta,&Tc,&DtaT);
  Tint = Ta;
  initial_conditions(T,Tn,nx,Ta); 
  error = 1;
  criterio = 1e-3;
  dias = 0;
  Tintmaxo = Tintmaxf = -100.;
  Tintmino = Tintminf =  100.;
  double Ein;
  double T1,T2,Tave,Keq;
  double Tavemax= -1000.,Tavemin=1000.;
  double **TK,DeltaTK;
  int nT = 10,*contadorTK,n;
  double contador_total = 0.;

  TK = two_d_double_array(nT,2);
  contadorTK = one_d_int_array(nT);
  for (n = 0; n < nT; ++n){
    TK[n][0] = 0.;
    TK[n][1] = 0.;
    contadorTK[n] = 0;
  }
  int Days = 15;
  
  while (dias < Days)  {
    ++dias;
    Qin = Tintaverage  = 0.;    Tintmax = Tsamax = -100.;  Tintmin = Tsamin = 100.;
    TPIhot = TPIcold = NumHot = NumCold = DenHot = DenCold = DDHhot = DDHcold = 0.;
    abrefile_muro(A,gamma,beta,&f_in,Lat,month);
    Econd = Econv = Erad = 0.;
    for (t = 0; t <= t_max; t += dt) { 
      Tsa = time_evolution_Tsa(&Ta,t,Ig,A,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta);
      Tsa1 = time_evolution_Tsa(&Ta,t,Ig,1.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta);
      Tsa0 = time_evolution_Tsa(&Ta,t,Ig,0.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta);
      convective_coefficients(&hi,beta,T[nx-1],Tint);
      //recalculate k for second layer
      set_krhocv3(k,rhoc,dx,L1,L2,L3,k1,k2,T,rhoc1,rhoc2,&RR,&Hr,&Hc,Aa,Ac,a11,a12,a21,Fc,&Econd,&Econv,&Erad,&T1,&T2,&Tave,&Keq);
      calculate_coefficients(a,b,c,d,dt,dx,k,nx,rhoc,T,Tsa,ho,Tint,hi);
      Ein = 0.;
      solve_PQ(a,b,c,d,P,Q,Tn,T,nx,&Tint,hi,rhoair,cair,La,&Qin,dt,&Tintaverage,&Ein);
      max_min(&Tsamax,&Tsamin,&Tintmax,&Tintmin,&t_Tintmax,&t_Tsamax,t,Tsa,Tint);
      discomfort_degree_hours(Tint,Tc,dt,&DDHhot,&DDHcold);
      if (Tint < Tc) NumCold += Tc - Tint; if (Tsa0 < Tc) DenCold += Tc - Tsa0; 
      if (Tint > Tc) NumHot  += Tint - Tc; if (Tsa1 > Tc) DenHot  += Tsa1 - Tc; 
      if (fmod(t,10.)<dt) {
	//rintf(f_in,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%f\t%f\t%         1       2    3  4   5  6  7  8  9    10
	fprintf(f_in,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%f\t%f\t%f\t%f\n",t/3600.,Tint,Ta,Tsa,Hr,Hc,T1,T2,Tave,Keq);
      }
      if (t == 0)  	for (i = 0; i < nx; ++i)      To[i] = T[i];
      if (t == t_max) 	for (i = 0; i < nx; ++i)      Tf[i] = T[i];
      if (dias ==Days-1) {
	//localiza m'aximo y min de Tave
	if (Tave > Tavemax) Tavemax = Tave;
	if (Tave < Tavemin) Tavemin = Tave;
	//printf("aqu'i?\n");
      }
      if (dias == Days) { 
	DeltaTK = (Tavemax - Tavemin)/(nT-1);
	if (Tave == Tavemin) {
	  TK[0][0] = Tavemin;
	  TK[0][1] = Keq;
	}
	for (n = 0; n < nT; ++n) {
	  if      (Tave >= ( Tavemin + DeltaTK*n - (DeltaTK/miniDeltaTK) ) && Tave <= ( Tavemin + DeltaTK*n + (DeltaTK/miniDeltaTK) ) ) {
	    TK[n][0] = Tave;
	    TK[n][1] += Keq;
	    ++contadorTK[n];
	    //printf("contadorTK[%d] = %d\n",n,contadorTK[n]);getchar();
	  }
	}
	if (Tave == Tavemax) {
	  /* TK[nT-1][0] = Tavemax; */
	  /* TK[nT-1][1] = Keq; */
	}
      } 
    }  

    
    Tintminf = Tintmin;
    Tintmaxf = Tintmax;
    error = error_maxmin(Tintmaxo,Tintmino,Tintmaxf,Tintminf,criterio);
    Tintmino = Tintmin;
    Tintmaxo = Tintmax;
    
    fclose(f_in);
    TPIhot = (1. - NumHot/DenHot)*100.;    TPIcold = (1. - NumCold/DenCold)*100.;
    //error = calculate_error(To,Tf,nx,criterio);
    //++error;
  }

  /* for (n = 0; n < nT; ++n)  */
  /*   printf("%f\t%f\t%d\n",TK[n][0],TK[n][1],contadorTK[n]); */
  
  /* for (n = 0; n < nT; ++n)  */
  /*   TK[n][1] = TK[n][1]/ (double)contadorTK[n]; */
  for (n = 0; n < nT; ++n) {
    TK[n][1] = TK[n][1]/ (double)contadorTK[n];
    contador_total += contadorTK[n];
  }
  
  
  printf("1D  Aa = %f\n\n",Aa); 
  printf("Rr = %f\n\n",RR);
  printf("DF\tLT\tET\n");
  printf("%.2f\t",(Tintmax-Tintmin)/(Tsamax-Tsamin));                               
  printf("%.1f\t",(t_Tintmax-t_Tsamax)/3600.);
  printf("%f\t",Qin/3600.);
  printf("%f\t%f\t%f\n",Econd/(Econd+Econv+Erad),Econv/(Econd+Econv+Erad),Erad/(Econd+Econv+Erad));

  double keq_promedio=0., keq_pesada = 0.;
  
  char file_TK[180];
  char orientation[140];
  for (n = 0; n < nT; ++n) {
    keq_promedio += TK[n][1]/nT;
    keq_pesada   += TK[n][1] * contadorTK[n]/contador_total;
  }


  if (gamma == 0. )   sprintf(orientation,"sur");
  if (gamma == 180. ) sprintf(orientation,"norte");
  if (gamma == -90. ) sprintf(orientation,"este");
  if (gamma == 90. )  sprintf(orientation,"oeste");
 
  printf("mes = %d\n",month);
  
  sprintf(file_TK,"./dat/TK_muro_Lat%.2f_a%.1f_mes%d_%s.csv",
	  Lat,A,month,orientation);
  
  f_in = fopen(file_TK,"w"); 
  printf(file_TK);
  printf("\n");
  for (n = 0; n < nT; ++n) 
    //fprintf(f_in,"%f\t%f\t%d\n",TK[n][0],TK[n][1],contadorTK[n]);
    //fprintf(f_in,"%f\n%f\n",TK[n][0],TK[n][1]);
    //fprintf(f_in,"%f\t%f\t%f\n%f\t%f\t%f\n",TK[n][0],keq_promedio,keq_pesada,TK[n][1],keq_promedio,keq_pesada);
    fprintf(f_in,"%f\t%f\n",TK[n][0],TK[n][1]);
}
 


