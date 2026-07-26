/**************** conduction.c para no aire acondicionado 1D  *********/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <libpq-fe.h>

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
  double f;
  int month,nx,cont_herramienta;
  char user[15];
  char baseDatos[15];
  
  if (argv[1] == NULL || argv[1][0] == '-') 
    {
      strcpy(inpfile, "conduction.e");
      strcat(inpfile, ".inp");
    } 
  else 
    strcpy(inpfile, argv[1]);
  input_reset();
  //input_insert("dt", "delta  t [s]", &dt, 'd');
  input_insert("N", "Number of elements", &nx, 'i');
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
  input_insert("f", "factor time (0,1]", &f, 'd');
  input_insert("k1", "thermal conductivity 1", &k1, 'd');
  input_insert("k2", "thermal conductivity 2", &k2, 'd');
  input_insert("k3", "thermal conductivity 3", &k3, 'd');
  input_insert("k4", "thermal conductivity 4", &k4, 'd');
  input_insert("k5", "thermal conductivity 5", &k5, 'd');
  input_insert("k6", "thermal conductivity 6", &k6, 'd');
  input_insert("k7", "thermal conductivity 7", &k7, 'd');
  input_insert("rhoc1", "thermal capacity 1", &rhoc1, 'd');
  input_insert("rhoc2", "thermal capacity 2", &rhoc2, 'd');
  input_insert("rhoc3", "thermal capacity 3", &rhoc3, 'd');
  input_insert("rhoc4", "thermal capacity 4", &rhoc4, 'd');
  input_insert("rhoc5", "thermal capacity 5", &rhoc5, 'd');
  input_insert("rhoc6", "thermal capacity 6", &rhoc6, 'd');
  input_insert("rhoc7", "thermal capacity 7", &rhoc7, 'd');
  input_insert("L1", "L1 [m]", &L1, 'd');
  input_insert("L2", "L2 [m]", &L2, 'd');
  input_insert("L3", "L3 [m]", &L3, 'd');
  input_insert("L4", "L4 [m]", &L4, 'd');
  input_insert("L5", "L5 [m]", &L5, 'd');
  input_insert("L6", "L6 [m]", &L6, 'd');
  input_insert("L7", "L7 [m]", &L7, 'd');
  input_insert("ho", "h outside [W/m2 oC]", &ho, 'd');
  input_insert("hi", "h inside [W/m2 oC]", &hi, 'd');
  input_insert("user", "usuario", &user, 's');
  input_insert("mes", "mes", &month, 'i');
  input_insert("Lo", "Longitud", &Lon, 'd');
  input_insert("Lat", "Latitud", &Lat, 'd');
  input_insert("cont", "contador", &cont_herramienta, 'i');
  input_insert("bd", "baseDatos", &baseDatos, 's');
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
  double Tsa1,Tsa0,Is;
  double TPIcold,TPIhot,NumCold,NumHot,DenCold,DenHot;
  double DDHcold,DDHhot;
  double t_Tintmax,t_Tsamax;
  FILE *f_in;

  int *juliano;
  double *delta,*orto, *ocaso,*Ho, tau[13];
  
  
  
  database_begin(user,cont_herramienta,month,baseDatos);

  

  dt = 1.;
  t_max  = 86400.;
  rhoair = 1.1797660470258469;
  cair = 1005.458757;  
  dx = (L1+L2+L3+L4+L5+L6+L7)/nx;
  
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

  while (error > 0)  {
    ++dias;
    Qin = Tintaverage  = 0.;    Tintmax = Tsamax = -100.;  Tintmin = Tsamin = 100.;
    TPIhot = TPIcold = NumHot = NumCold = DenHot = DenCold = DDHhot = DDHcold = 0.;
    abrefile(user,month,cont_herramienta,&f_in);
    for (t = 0; t <= t_max; t += dt) { 
      Tsa = time_evolution_Tsa(&Ta,t,Ig,A,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta,&Is);
      Tsa1 = time_evolution_Tsa(&Ta,t,Ig,1.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta,&Is);
      Tsa0 = time_evolution_Tsa(&Ta,t,Ig,0.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta,&Is);
      convective_coefficients(&hi,beta,T[nx-1],Tint);
      calculate_coefficients(a,b,c,d,dt,dx,k,nx,rhoc,T,Tsa,ho,Tint,hi);
      solve_PQ(a,b,c,d,P,Q,Tn,T,nx,&Tint,hi,rhoair,cair,La,&Qin,dt,&Tintaverage);
      max_min(&Tsamax,&Tsamin,&Tintmax,&Tintmin,&t_Tintmax,&t_Tsamax,t,Tsa,Tint);
      discomfort_degree_hours(Tint,Tc,dt,&DDHhot,&DDHcold);
      if (Tint < Tc) NumCold += Tc - Tint; if (Tsa0 < Tc) DenCold += Tc - Tsa0; 
      if (Tint > Tc) NumHot  += Tint - Tc; if (Tsa1 > Tc) DenHot  += Tsa1 - Tc; 
      if (fmod(t,600.)<dt) 
	fprintf(f_in,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",Is,Tsa,Ta,T[0],T[nx-1],Tint,Tc,DtaT);
     
      if (t == 0)  	for (i = 0; i < nx; ++i)      To[i] = T[i];
      if (t == t_max) 	for (i = 0; i < nx; ++i)      Tf[i] = T[i];
      
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
  
  
  
  //  abrefile_indice(L1,L2,L3,k1,k2,k3,rhoc1,rhoc2,rhoc3,La,user,month,cont_herramienta,&f_in,Tsamax,Tsamin,Tintmax,Tintmin,t_Tintmax,t_Tsamax,Qin,Tintaverage,dt,TPIhot,TPIcold,DDHhot,DDHcold);

  abrefile_indice(user,cont_herramienta,&f_in,Tsamax,Tsamin,Tintmax,Tintmin,t_Tintmax,t_Tsamax,Qin,Tintaverage,dt,TPIhot,TPIcold,DDHhot,DDHcold,month);


  database_end(user,cont_herramienta,month,baseDatos);
}



