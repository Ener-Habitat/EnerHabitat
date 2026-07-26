/*********** conduction.c para no aire acondicionado 2D  **********/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <libpq-fe.h>
/*********************************/

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
  double E;
  double a11,a12,a13,a14;
  double a21,a22,a23;
  double e21,e22,e23;
  int layer;
  double ho,hi;
  double k1,k2,k3,k4,k5,k6,k7,kr;
  double rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7,rhocr;
  double t_Tamax,Tmax,Tmin,Ig,Id,Ib,beta,gamma,A;
  double Lat,Lon;
  int month,ny,nx,cont_herramienta,tipo;
  char user[15],mail[200];
  char baseDatos[15];
  
  if (argv[1] == NULL || argv[1][0] == '-') 
    {
      strcpy(inpfile, "conduction.e");
      strcat(inpfile, ".inp");
    } 
  else 
    strcpy(inpfile, argv[1]);
  input_reset();
  input_insert("ny", "Number of elements", &ny, 'i');
  input_insert("nx", "Number of elements", &nx, 'i');
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
  input_insert("k1", "thermal conductivity 1", &k1, 'd');
  input_insert("k2", "thermal conductivity 2", &k2, 'd');
  input_insert("k3", "thermal conductivity 3", &k3, 'd');
  input_insert("k4", "thermal conductivity 4", &k4, 'd');
  input_insert("k5", "thermal conductivity 5", &k5, 'd');
  input_insert("k6", "thermal conductivity 6", &k6, 'd');
  input_insert("k7", "thermal conductivity 7", &k7, 'd');
  input_insert("kr", "thermal conductivity fill", &kr, 'd');
  input_insert("rhoc1", "thermal capacity 1", &rhoc1, 'd');
  input_insert("rhoc2", "thermal capacity 2", &rhoc2, 'd');
  input_insert("rhoc3", "thermal capacity 3", &rhoc3, 'd');
  input_insert("rhoc4", "thermal capacity 4", &rhoc4, 'd');
  input_insert("rhoc5", "thermal capacity 5", &rhoc5, 'd');
  input_insert("rhoc6", "thermal capacity 6", &rhoc6, 'd');
  input_insert("rhoc7", "thermal capacity 7", &rhoc7, 'd');
  input_insert("rhocr", "thermal capacity fill", &rhocr, 'd');
  input_insert("L1", "L1 [m]", &L1, 'd');
  input_insert("L2", "L2 [m]", &L2, 'd');
  input_insert("L3", "L3 [m]", &L3, 'd');
  input_insert("L4", "L4 [m]", &L4, 'd');
  input_insert("L5", "L5 [m]", &L5, 'd');
  input_insert("L6", "L6 [m]", &L6, 'd'); 
  input_insert("L7", "L7 [m]", &L7, 'd');
  input_insert("layer", "#layer from outside", &layer, 'i');
  input_insert("e", "Emissivity", &E, 'd');
  input_insert("a11", "a11 [m]", &a11, 'd');
  input_insert("a12", "a12 [m]", &a12, 'd');
  input_insert("a13", "a13 [m]", &a13, 'd');
  input_insert("a14", "a14 [m]", &a14, 'd');
  input_insert("a21", "a21 [m]", &a21, 'd');
  input_insert("a22", "a22 [m]", &a22, 'd');
  input_insert("a23", "a23 [m]", &a23, 'd');
  input_insert("e21", "e21 [m]", &e21, 'd');
  input_insert("e22", "e22 [m]", &e22, 'd');
  input_insert("e23", "e23 [m]", &e23, 'd');
  input_insert("ho", "h outside [W/m2 oC]", &ho, 'd');
  input_insert("hi", "h inside [W/m2 oC]", &hi, 'd');
  input_insert("user", "usuario", &user, 's');
  input_insert("mes", "mes", &month, 'i');
  input_insert("Lo", "Longitud", &Lon, 'd');
  input_insert("Lat", "Latitud", &Lat, 'd');
  input_insert("cont", "contador", &cont_herramienta, 'i');
  input_insert("mail", "mail", &mail, 's');  
  input_insert("tipo", "tipo", &tipo, 'i');  
  input_insert("bd", "baseDatos", &baseDatos, 's'); 
  input_load(inpfile);
  if (input_options(argc,argv,help)) input();  
  input_save(inpfile);
  

  double dt,dx,dy,t,error,t_max,Tsa,Is,Tso;
  double **a,**b,**c,**d,**T,**To,**Tn,**Terror,*P,*Q;
  int **NT;
  double **rhoc,**k;
  int i,j,itera;
  double Hi;
  double Qin,rhoair,cair,Tint,Ta,Tc,DtaT,Thueco,hh,kair;
  double Qrup,Qrdown;
  double Tintmax,Tintmin,Tsamax,Tsamin,Tintaverage,Tsi;
  double Tsa1,Tsa0,Tsurface;
  double TPIcold,TPIhot,NumCold,NumHot,DenCold,DenHot;
  double DDHcold,DDHhot;
  double t_Tintmax,t_Tsamax;
  double X,Y;
  int *juliano;
  double *delta,*orto, *ocaso,*Ho, tau[13];  
  FILE *f_in;
  int i1,j1,i2,j2;
  int animation;
  double YY[8];
  double y1,y2;
  double Nur;

  //VARIABLES PARA VERIFICAR EL RAYLEIGH
  double Tarriba,Tabajo,Ra,gr,Beta,nu,alphaair;
  Tarriba = Tabajo = 0.;
  gr = 9.81;
  Beta = 1./300.;
  nu = 1.11e-5;
  rhoair = 1.668790;
  kair = 0.0262;
  cair = 1000.; 
   
  database_begin(user,cont_herramienta,month,baseDatos);
  
  YY[0] = 0.;
  YY[1] = L1;
  YY[2] = L2;
  YY[3] = L3;
  YY[4] = L4;
  YY[5] = L5;
  YY[6] = L6;
  YY[7] = L7;
  YY[layer] = e21 + e22 + e23;
  
  if (a14 == 0.) {   //THIS MEANS ONE HOLLOW
    X = a21 + a11 + a12/2.;
  }
  else {
    X = a11 + a21 + a12 + a22 + a13 + a23 + a14;
  }
  if (tipo == 4)
	X = 2.* a11 + a21;
	
  
  dt = 1.;
  t_max  = 86400.;
  rhoair = 1.1797660470258469;
  cair = 1005.458757; 
  dx = X/nx;
  Y = 0.;
  for (j = 0; j <= 7; ++j ) 
    Y += YY[j];
  dy = Y/ny;
  
  y1 = 0.;
  for (j = 0; j < layer ; ++j)
    y1 += YY[j];
  
  y2 = (y1 + YY[layer])/dy + 0.5;
  y1 = y1/dy + 0.5;

  i1 = a11/dx + 0.5 ;
  //j1 = e21/dy + 0.5 + y1 ;
  j1 = e21/dy  + (int) y1 ;
  i2 = (a11+a21)/dx + 0.5;
  //j2 = (e21+e22)/dy + 0.5 + y1;
  j2 = (e21+e22)/dy  + (int) y1;

  NT = two_d_int_array(nx,ny);
  T  = two_d_double_array(nx,ny);
  To  = two_d_double_array(nx,ny);
  Tn = two_d_double_array(nx,ny);
  Terror = two_d_double_array(nx,ny);
  k = two_d_double_array(nx,ny);
  rhoc = two_d_double_array(nx,ny);
  a = two_d_double_array(nx,ny);
  b = two_d_double_array(nx,ny);
  c = two_d_double_array(nx,ny);
  d = two_d_double_array(nx,ny);
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

  if (tipo == 1) {
    //printf("vigueta y bovedilla hueca\n");
    draw_viguetabovedilla2hueca(nx,ny,NT,i1,j1,i2,j2); 
    set_krhoc(k,rhoc,dx,dy,nx,ny,L1,L2,L3,L4,L5,L6,L7,k1,k2,k3,k4,k5,k6,k7,rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7);
  }

  if (tipo == 2) {
    //printf("vigueta y bovedilla rellena\n");
    draw_viguetabovedilla2rellena(nx,ny,NT,i1,j1,i2,j2);
    set_krhocrelleno(k,rhoc,dx,dy,nx,ny,L1,L2,L3,L4,L5,L6,L7,k1,k2,k3,k4,k5,k6,k7,kr,rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7,rhocr,NT);
  }

  if(tipo==4) {
    i1 = a11/dx/2. + 0.5 ;
    j1 = e21/dy  + (int) y1 ;
    i2 = (a11/2.+a21)/dx + 0.5;
    j2 = (e21+e22)/dy  + (int) y1;
    draw_viguetabovedillarellena(nx,ny,NT,i1,i2,j1,j2); 
    set_krhoc_viguetabovedillarellena(k,rhoc,dx,dy,nx,ny,L1,L2,L3,L4,L5,L6,L7,k1,k2,k3,k4,k5,k6,k7,kr,rhoc1,rhoc2,rhoc3,rhoc4,rhoc5,rhoc6,rhoc7,rhocr,NT);
  } 



  Ta_Tc_DtaT(Ig,ho,Tmax,Tmin,t_Tamax,Ho[month],&Ta,&Tc,&DtaT);
  /* printf("DtaT = %.2f\n",DtaT); */
  /* printf("Tc = %.2f\n",Tc); */
  error = 10.;
  Tint = Tc+DtaT;
  Thueco = Tint;
  initial_conditions(T,Tn,Terror,To,nx,ny,Tint); 
  
  while (error > 1e-5)  { 
    Qin = Tintaverage  = Qrup = Qrdown = 0.;    Tintmax = Tsamax = -100.;  Tintmin = Tsamin = 100.; animation = 0;
    TPIhot = TPIcold = NumHot = NumCold = DenHot = DenCold  = DDHhot = DDHcold = 0.;
    abrefile(user,month,cont_herramienta,&f_in);
    for (t = 0; t <= t_max; t += dt) {
      Tsa = time_evolution_Tsa(&Ta,t,Ig,A,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta,&Is);
      Tsa1 = time_evolution_Tsa(&Ta,t,Ig,1.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta,&Is);
      Tsa0 = time_evolution_Tsa(&Ta,t,Ig,0.,ho,Tmax,Tmin,t_Tamax,Ho[month],tau[month]*2.,delta[month],Ib,Id,gamma,Lat,beta,&Is);
      interchange(To,T,nx,ny);
      Tsi = Tsint (T,nx,ny);
      Tso = Tsout (T,nx,ny);
      convective_coefficients(&hi,beta,Tsurface,Tint);
      solve_PQ(a,b,c,d,P,Q,Tn,T,nx,ny,&Tint,hi,rhoair,cair,La,&Qin,dt,dx,dy,k,rhoc,NT,Tsa,ho,To,X,t,&Tintaverage,&Thueco,&hh,i1,j1,i2,j2,a21,e22,E,&Qrup,&Qrdown,tipo,beta,&Tarriba,&Tabajo,&Nur);
      max_min(&Tsamax,&Tsamin,&Tintmax,&Tintmin,&t_Tintmax,&t_Tsamax,t,Tsa,Tint,&Tsi,T,nx,ny);
      discomfort_degree_hours(Tint,Tc,dt,&DDHhot,&DDHcold);
      if (Tint < Tc) NumCold += Tc - Tint;  if (Tsa0 < Tc) DenCold += Tc - Tsa0;    
      if (Tint > Tc) NumHot  += Tint - Tc;  if (Tsa1 > Tc) DenHot  += Tsa1 - Tc;  
      if (fmod(t,600.)<dt)  {
	fprintf(f_in,"%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",Is,Tsa,Ta,Tso,Tsi,Tint,Tc,DtaT);
	//printf("%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",t/3600.,Is,Tsa,Ta,Tso,Tsi,Tint,Tc,DtaT);
      }
    }
    
    fclose(f_in);
    //printf("ya cerro el archivo\n");getchar();
    TPIhot = (1. - NumHot/DenHot)*100.;    TPIcold = (1. - NumCold/DenCold)*100.;	
    error = calculate_error(T,Terror,nx,ny);
    interchange(Terror,T,nx,ny);
  }
  printf("TPIhot  = %f\n",TPIhot);
  printf("TPIcold = %f\n",TPIcold);
  printf("DDHcold = %f\n",DDHcold);
  printf("DDHhot = %f\n",DDHhot);
  
  abrefile_indice(user,cont_herramienta,month,&f_in,
		  Tsamax,Tsamin,Tintmax,Tintmin,t_Tintmax,
		  t_Tsamax,Qin,Tintaverage,dt,X,TPIhot,TPIcold,DDHhot,DDHcold);
  //animation_gnuplot(animation,X,Y,Tintmax,Tintmin); 
  
  
 database_end(user,cont_herramienta,month,baseDatos);
  char orden[200];
  
  sprintf(orden,"java -classpath /var/lib/tomcat6/webapps/Cie/WEB-INF/lib/mail.jar:. Mailparametros %s",mail);
  system(orden);
}




