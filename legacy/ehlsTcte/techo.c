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
  double kbeam,khollow,ktop,kfinish;
  double rhocbeam,rhochollow,rhoctop,rhocfinish;
  double t_Tamax,Tmax,Tmin,Ig,Id,Ib,beta,gamma,A;
  double Lat,Lon;
  double Fc;
  double d1,d2,d3;
  int miniDeltaTK;
  int month,nx,cont_herramienta;
  char user[15];
  
  if (argv[1] == NULL || argv[1][0] == '-') 
    {
      strcpy(inpfile, "techo.e");
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
  //input_insert("Ig", "Global solar radiation", &Ig, 'd');
  input_insert("Id", "Diffuse solar radiation", &Id, 'd');
  input_insert("Ib", "Beam solar radiation", &Ib, 'd');
  input_insert("beta", "Wall inclination [degree]", &beta, 'd');
  input_insert("gamma", "Wall orientation [degree]", &gamma, 'd');
  input_insert("a", "Absortance ", &A, 'd');
  input_insert("La", "Room lenght", &La, 'd');
  input_insert("Fc", "Correction factor(0,1]", &Fc, 'd');
  input_insert("d1", "d1 [m]", &d1, 'd');
  input_insert("d2", "d2 [m]", &d2, 'd');
  input_insert("d3", "d3 [m]", &d3, 'd');
  input_insert("kbeam",   "thermal conductivity T-beam       [W/mK]", &kbeam,   'd');
  input_insert("khollow", "thermal conductivity hollow-block [W/mK]", &khollow,  'd');
  input_insert("ktop",    "thermal conductivity top-concrete [W/mK]", &ktop,    'd');
  input_insert("kfinish", "thermal conductivity plaster      [W/mK]", &kfinish, 'd');
  input_insert("rhocbeam",   "thermal capacitance T-beam       []", &rhocbeam,   'd');
  input_insert("rhochollow", "thermal capacitance hollow-block []", &rhochollow,  'd');
  input_insert("rhoctop",    "thermal capacitance top-concrete []", &rhoctop,    'd');
  input_insert("rhocfinish", "thermal capacitance plaster      []", &rhocfinish, 'd');
  input_insert("L1", "L1 [m]", &L1, 'd');
  input_insert("L2", "L2 [m]", &L2, 'd');
  input_insert("L3", "L3 [m]", &L3, 'd');
  input_insert("L4", "L4 [m]", &L4, 'd');
  input_insert("L5", "L5 [m]", &L5, 'd');
  input_insert("L6", "L6 [m]", &L6, 'd');
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
   double Qin,rhoair,cair,Tint,Ta,Tc,DtaT;

   dt = 1.; 
   t_max  = 86400.; 
   rhoair = 1.205; 
   cair = 1005.;  
   dx = (L1+L2+L3+L4+L5+L6)/nx;
   //printf("dx = %f\n",dx);
   
   
  int error,dias;
  double *a,*b,*c,*d,*T,*Tn,*To,*Tf,*P,*Q;
  double *rhoc,*k;
  int i,j,itera;
  double Hi;
  double Tintmax,Tintmin,Tsamax,Tsamin,Tintaverage;
  double Tsa1,Tsa0;
  double TPIcold,TPIhot,NumCold,NumHot,DenCold,DenHot;
  double DDHcold,DDHhot;
  double t_Tintmax,t_Tsamax;
  FILE *f_in;

  int *juliano;
  double *delta,*orto, *ocaso,*Ho, tau[13];
  
  L7    = 0.;
  



  //DEFINIENDO LA GEMETR'IA DEL TECHO
  double Ceq,RR,Hr,Hc,Aa,Ab,At;
  double d4,d5,d6,d7,d8,d9,D;
  double Ltotal;

  d5 = d7 = d3;
  d4 = d6 = d8 = d2;
  d9 = d1;
  D = d1 + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9;
  Ltotal = L1 + L2 + L3 + L4 + L5 + L6 + L7;

  Aa = (d3 + d5 + d7 )/D;
  Ab = (d2 + d4 + d6 + d8 )/D;
  At = (d1 + d9)/D;
    
  printf("Ltotal = %f\n",Ltotal);
  printf("D      = %f\n",D);
  printf("\n");
  printf("Aa     = %f\n",Aa);
  printf("Ab     = %f\n",Ab);
  printf("At     = %f\n",At);
  //CALCULATION OF THE EQUIVALENT THERMAL CAPACITANCE FOR LAYER 3,4,5,6
  //Debe haber Ceq para las capas 3,4,5 y 6
  double Ceq3,Ceq4,Ceq5,Ceq6;
  Ceq3 = (d1/2. + d9/2.)/D*rhocbeam + (d1/2. + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9/2.)/D*rhoctop;
  Ceq4 = (d1/2. + d9/2.)/D*rhocbeam + (d1/2. + d2 + d3 + d4 + d5 + d6 + d7 + d8 + d9/2.)/D*rhochollow;
  Ceq5 = (d1/2. + d9/2.)/D*rhocbeam + (d1/2. + d2 + d4 + d6 +d8 + d9/2. )/D*rhochollow
       + (d3 + d5 + d7)/D *rhoair*cair;
  Ceq6 = (d1 + d9)/D*rhocbeam + (d2 + d3 + d4 + d5 + d6 + d7 + d8)/D*rhochollow;
  printf("\n");
  printf("Ceq3 = %.3e\n",Ceq3);
  printf("Ceq4 = %.3e\n",Ceq4);
  printf("Ceq5 = %.3e\n",Ceq5);
  printf("Ceq6 = %.3e\n",Ceq6);
  //CALCULATION OF THE EQUIVALENT THERMAL RESISTANCE FOR LAYERS 3,4,6 NOT 5
  double Req3,Req4,Req6;
  Req3 = pow( 1./(L3/(d1/(2.*D))/kbeam) + 1./(L3/( (d1/2.+d2+d3+d4+d5+d6+d7+d8+d9/2.)/D)/ktop)
	      + 1./(L3/(d9/(2.*D))/kbeam),-1.);
  Req4 = pow( 1./(2.*L4*D/d1/kbeam) + 1./(L4*D/(d1/2.+d2+d3+d4+d5+d6+d7+d8+d9/2.)/khollow)
	      +1./(2.*L4*D/d9/kbeam), -1.);
  Req6 = pow( 1./(L6*D/d1/kbeam) + 1./(L6*D/(d2 + d3 + d4 + d5 + d6 + d7 + d8)/khollow)
	      + 1./(L6*D/d9/kbeam),-1.);
  //CALCULATION OF THE EQUIVALENT THERMAL CONDUCTIVITY FOR LAYERS 3,4,6 NOT 5
  double keq3,keq4,keq6;
  keq3 = L3/Req3;
  keq4 = L4/Req4;
  keq6 = L6/Req6;
  printf("\n");
  printf("keq3 = %f\n",keq3);
  printf("keq4 = %f\n",keq4);
  printf("keq6 = %f\n",keq6);
  

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
  dia_juliano(juliano,15,2011,12);   //dia, anio, mes
  calculo_declinacion_delta(delta,juliano,12);
  calculo_orto_ocaso(orto, ocaso, Lat, delta, 12);
  calculo_duracion_dia(tau, 12, orto);
  calculo_hora_orto(Ho, tau, 12);
  double k1,k2,k3,k4,k5,k6,k7;
  k1 = kfinish;
  k2 = ktop;
  k3 = keq3;
  k4 = keq4;
  k5 = 2.;
  k6 = keq6;
  k7 = 0.;
  double rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7;
  rhoc1 = rhocfinish;
  rhoc2 = rhoctop;
  rhoc3 = Ceq3;
  rhoc4 = Ceq4;
  rhoc5 = Ceq5;
  rhoc6 = Ceq6;
  rhoc7 = 0.;
  set_krhoc(k,rhoc,dx,L1,L2,L3,L4,L5,L6,L7,k1,k2,k3,k4,k5,k6,k7,rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7); 


  double Econd,Econv,Erad;


  double Tintmaxo,Tintmino,Tintmaxf,Tintminf;


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
  int Days = 20;
  
  while (dias < Days)  {
    ++dias;
    Qin = Tintaverage  = 0.;    Tintmax = Tsamax = -100.;  Tintmin = Tsamin = 100.;
    TPIhot = TPIcold = NumHot = NumCold = DenHot = DenCold = DDHhot = DDHcold = 0.;
    abrefile_techo(A,gamma,beta,&f_in,Lat,month);
    Econd = Econv = Erad = 0.;
    
    for (n = 0; n < nT; ++n){
      TK[n][0] = 0.;
      TK[n][1] = 0.;
      contadorTK[n] = 0;
    }
    for (t = 0; t <= t_max; t += dt) {
      Tsa = time_evolution_Tsa(&Ta,t,Ig,A,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta);
      Tsa1 = time_evolution_Tsa(&Ta,t,Ig,1.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta);
      Tsa0 = time_evolution_Tsa(&Ta,t,Ig,0.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta);
      convective_coefficients(&hi,beta,T[nx-1],Tint);
      //recalculate k for second layer
      set_krhocv3_techo(k,rhoc,dx,
			L1,L2,L3,L4,L5,L6,
			k1,k2,k3,k4,k5,k6,
			kfinish,ktop,kbeam,khollow,
			T,
			rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,
			&RR,&Hr,&Hc,
			d1,d2,d3,d4,d5,d6,d7,d8,d9,
			(d3+d5+d7)/D,(d1+d2+d4+d6+d8+d9)/D,
			Fc,&Econd,&Econv,&Erad,&T1,&T2,&Tave,&Keq);
      calculate_coefficients(a,b,c,d,dt,dx,k,nx,rhoc,T,Tsa,ho,Tint,hi);
      Ein = 0.;
      solve_PQ(a,b,c,d,P,Q,Tn,T,nx,&Tint,hi,rhoair,cair,La,&Qin,dt,&Tintaverage,&Ein);
      max_min(&Tsamax,&Tsamin,&Tintmax,&Tintmin,&t_Tintmax,&t_Tsamax,t,Tsa,Tint);
      discomfort_degree_hours(Tint,Tc,dt,&DDHhot,&DDHcold);
      if (Tint < Tc) NumCold += Tc - Tint; if (Tsa0 < Tc) DenCold += Tc - Tsa0;
      if (Tint > Tc) NumHot  += Tint - Tc; if (Tsa1 > Tc) DenHot  += Tsa1 - Tc;
      if (fmod(t,10.)<dt) {
	//                                                                     1      2  3  4   5  6  7  8    9   10     
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
  	for (n = 0 ; n <= nT; ++n) {
  	  if      (Tave > ( Tavemin + DeltaTK*n - (DeltaTK/miniDeltaTK) ) && Tave < ( Tavemin + DeltaTK*n + (DeltaTK/miniDeltaTK) ) ) {
  	    TK[n][0] = Tave;
  	    TK[n][1] += Keq;
  	    ++contadorTK[n];
  	    //printf("contadorTK[%d] = %d\n",n,contadorTK[n]);getchar();
  	  }
  	}
  	if (Tave == Tavemax) {
	  // TK[nT][0] = Tavemax;
  	  //TK[nT][1] = Keq;
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
  }
  
  for (n = 0; n < nT; ++n) {
    TK[n][1] = TK[n][1]/ (double)contadorTK[n];
    contador_total += contadorTK[n];
  }

  //printf("contador_total = %f\n",contador_total);
  
  printf("1D  Aa = %f\n\n",Aa);
  printf("Rr = %f\n\n",RR);
  printf("DF\tLT\tET\n");
  printf("%.2f\t",(Tintmax-Tintmin)/(Tsamax-Tsamin));
  printf("%.1f\t",(t_Tintmax-t_Tsamax)/3600.);
  printf("%f\n",Qin/3600.);

  double keq_promedio=0., keq_pesada = 0.;
  
  for (n = 0; n < nT; ++n) {
    keq_promedio += TK[n][1]/nT;
    keq_pesada   += TK[n][1] * contadorTK[n]/contador_total;
  }


  
  char file_TK[180];
  sprintf(file_TK,"./dat/TK_techo_Lat%.2f_a%.1f_mes%d.csv",
	  Lat,A,month);

  
  f_in = fopen(file_TK,"w");
  printf(file_TK);
  printf("\n");
  fprintf(f_in,"#Ktabla\tKeq_promedio\tKeq_pesada\n");
  for (n = 0; n < nT; ++n)
    //fprintf(f_in,"%f,%f,%f,%f,%f\n",TK[n][0],TK[n][1],contadorTK[n]/contador_total,keq_promedio,keq_pesada);
    fprintf(f_in,"%f\t%f\t%f\n%f\t%f\t%f\n",TK[n][0],keq_promedio,keq_pesada,TK[n][1],keq_promedio,keq_pesada);
  
}
 


