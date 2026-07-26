/*8 de septiembre del 2014, Temixco, Morelos, M'exico       */
/*                                          */
#define _XOPEN_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>


//define pointers used to monitor values in the monitor.h file
FILE *f_in;
FILE *f_in2;
FILE *f_in3;
 
double DF_ave,LT_ave,DF_sigma,LT_sigma,Tave,Tin_sigma;
#include "arrays.h"            // includes the routine to make dynamic arrays
#include "tools.h"             // tools or rutines for the code
#include "new_input.h"         // includes the routines to start the executable

/******************************************************************************************
*******************************************************************************************/


static char help[]="\
***************************************************************************\n\
datos.c\n\
***************************************************************************\n\
\n\
Calcula el retraso A, tiempos de retraso m'aximos y minimos TR \n\
***************************************************************************\n\
 ";

main(int argc,char *argv[]) 
{
  

  /* Routine to save initial values in file "program.e.inp" */   
  /**********************************************************/
  char inpfile[120],file1[140];
  int m1,m2;
  int d1,d2;
  double Lon,Lat,msnm;
  if (argv[1] == NULL || argv[1][0] == '-')
    {
      strcpy(inpfile, "datos.e");
      strcat(inpfile, ".inp");
    } 
  else 
    strcpy(inpfile, argv[1]);
  /************************************************************************/
  input_reset();
  input_insert("file", "archivo datos", file1, 's');
  input_insert("Lon", "Longitud", &Lon, 'd');
  input_insert("Lat", "Latitud", &Lat, 'd');
  input_insert("msnm", "Altura snm", &msnm, 'd');
  /************************************************************************/
  input_load(inpfile);
  if (input_options(argc,argv,help)) input();  
  input_save(inpfile);
  /* Ends routine to save initial values in file "program.e.inp"  */
  /****************************************************************/
  
  if  ((f_in=fopen(file1,"r"))== NULL) {
    printf("Can't open file\n");
    exit(-1);
  }
  
  FILE    *F;
  struct tm f,f2,f3,ftmp;
  time_t t,t2,t3;
  int lines,i,j;
  double pi;

  pi =acos(-1.);
  
  //NUEVAS VARIABLES
  double **Ig,**Ib,**Id;
  double *Tmax,*Tmin,*tTmax,*tTmin;
  
  Ig    = two_d_double_array(12,24);
  Ib    = two_d_double_array(12,24);
  Id    = two_d_double_array(12,24);
  Tmax  = one_d_double_array(12);
  tTmax = one_d_double_array(12);
  Tmin  = one_d_double_array(12);
  tTmin = one_d_double_array(12);
  
  
  for (i = 0; i < 12; ++i) {
    Tmax[i] = tTmax[i] = Tmin[i] = tTmin[i] = 0.;
    for (j = 0; j < 24; ++j) 
      Ig[i][j] = Ib[i][j] = Id[i][j] = 0.;
  }
  
  //VARIABLES PARA EVITAR LOS PRIMEROS RENGLONES
  char renglon[1000],espacio[1000],espacio2[1000];
 
  double year,mes,dia,hora,Ta,Tamax,Tamin,tTamax,tTamin;
  double G,I,D;
  double *counter_month;
  
  Tamax = -1000.;
  Tamin =  1000.;
  
  counter_month = one_d_double_array(12);
  for (i = 0; i < 12; ++i) 
    counter_month[i] = 0.;
  
  int dtmp,position,contador=0;
  //ABRE EL ARCHIVO CON DATOS DEL INTERIOR Y EXTERIOR PARA E+
  if((F=fopen(file1,"r")) != NULL) {
    for (i = 0; i < 8; ++i) { 
      if ( fgets (renglon , 500 , F) != NULL );
    }
    while (!feof(F)) { 
      //A~no,mes,dia,hora,minuto,*?,
      fscanf(F,"%lf,%lf,%lf,%lf,%*lf,%44s,%lf,%*lf,%*lf,%*lf,%*lf,%*lf,%*lf,%lf,%lf,%lf,%500s"
	     ,&year,&mes,&dia,&hora,espacio,&Ta,&G,&I,&D,espacio2);
      m2 = mes;
      d2 = dia;
      //fscanf(F,"%lf,%lf,%lf,%lf,%*lf,%s",&year,&mes,&dia,&hora,espacio);
      //printf("%.0f/%.0f/%.0f\t%.0f\t%f\t%f\t%f\n",year,mes,dia,hora,G,I,D);//getchar();
      //printf("contador_mes[%d] = %f\n",m2-1,counter_month[m2-1]);getchar();
      ++counter_month[m2-1];
      memset(&f, 0, sizeof(f));   //llena de ceros las estructuras
      f.tm_mday  = dia;
      f.tm_mon   = mes-1;
      f.tm_year = year - 1900; 
      t=mktime(&f);
      dtmp = f.tm_yday;
      position = hora -  1;
      Ig[m2-1][position] += G; 
      Ib[m2-1][position] += G-D;
      Id[m2-1][position] += D;
      if (Ta > Tamax) { Tamax = Ta;   tTamax = hora; }
      if (Ta < Tamin) { Tamin = Ta;   tTamin = hora; }
      if (hora == 24.) { 
	Tmax[m2-1]  = Tmax[m2-1]  + Tamax;
	Tmin[m2-1]  = Tmin[m2-1]  + Tamin;
	tTmax[m2-1] = tTmax[m2-1] + tTamax;
	tTmin[m2-1] = tTmin[m2-1] + tTamin;
	//printf("Tmax = %f\n",Tamax);getchar();
	Tamax =-1000.; Tamin = 1000.;
      }
      
    }
  }
  fclose(F);




  for (i = 0 ; i< 12; ++i) {
    counter_month[i]  /= 24.;
    //printf("counter = %f\n",counter_month[i]);getchar();
    for (j = 0; j < 24; ++j) {
      Ig[i][j] /= counter_month[i];
      Ib[i][j] /= counter_month[i];
      Id[i][j] /= counter_month[i];
    }
  }
  
  double *Igmes,*Ibmes,*Idmes;

  Igmes = one_d_double_array(12);
  Ibmes = one_d_double_array(12);
  Idmes = one_d_double_array(12);
  for (i = 0; i < 12; ++i)      Igmes[i] = Ibmes[i] =  Idmes[i] = 0.;
  

  for (i = 0; i < 12; ++i) {
    Tmax[i]  /= counter_month[i];
    Tmin[i]  /= counter_month[i];
    tTmax[i] /= counter_month[i];
    tTmin[i] /= counter_month[i];

    /* printf("Tmax[%d] = %f\n",i,Tmax[i]); */
    /* printf("Tmin[%d] = %f\n",i,Tmin[i]); */
    /* printf("tTmax[%d] = %f\n",i,tTmax[i]); */
    /* printf("tTmin[%d] = %f\n",i,tTmin[i]); */
    /* getchar(); */
  }

  for (i = 0; i < 24; ++i)  {
    //printf("%d",i);
    for(j = 0; j <  12; ++j)  {
      Igmes[j] += Ig[j][i];
      Idmes[j] += Id[j][i];
      Ibmes[j] += Ib[j][i];
      //printf("\t%f",Ib[j][i]);
    }
    //printf("\n");
  }

  double *Igmax,*Ibmax,*Idmax;

  Igmax = one_d_double_array(12);
  Ibmax = one_d_double_array(12);
  Idmax = one_d_double_array(12);
  for (i = 0; i < 12; ++i)   Igmax[i] = Ibmax[i] = Idmax[i] = 0.;
  

  //printf("%s\n",file1);
  char f_out[180];
  char f_out2[180];
  char m[50];
  
  sprintf(f_out,"./dat/Radiacion%s.dat",file1);
  sprintf(f_out2,"./dat/Temperaturas%s.dat",file1);
  printf(f_out);
  printf("\n");
  printf(f_out2);
  printf("\n");
  f_in2 = fopen(f_out,"w");
  f_in3 = fopen(f_out2,"w");



  fprintf(f_in2,"##############\n");
  fprintf(f_in2,"lugar\t%s\n",file1);
  for (i = 0; i < 12; ++i) {
    if (i==0)  fprintf(f_in2,"mes\tEnero\n");
    if (i==1)  fprintf(f_in2,"mes\tFebrero\n");
    if (i==2)  fprintf(f_in2,"mes\tMarzo\n");
    if (i==3)  fprintf(f_in2,"mes\tAbril\n");
    if (i==4)  fprintf(f_in2,"mes\tMayo\n");
    if (i==5)  fprintf(f_in2,"mes\tJunio\n");
    if (i==6)  fprintf(f_in2,"mes\tJulio\n");
    if (i==7)  fprintf(f_in2,"mes\tAgosto\n");
    if (i==8)  fprintf(f_in2,"mes\tSeptiembre\n");
    if (i==9)  fprintf(f_in2,"mes\tOctubre\n");
    if (i==10) fprintf(f_in2,"mes\tNoviembre\n");
    if (i==11) fprintf(f_in2,"mes\tDiciembre\n");
    Igmax[i] = Igmes[i]*pi/24.;
    Ibmax[i] = Ibmes[i]*pi/24.;
    Idmax[i] = Idmes[i]*pi/24.;
    //printf("Ig= %f\t",i,Igmax[i]);
    fprintf(f_in2,"Ib\t%.2f\n",Ibmax[i]);
    fprintf(f_in2,"Id\t%.2f\n",Idmax[i]);
  }
  




  fprintf(f_in3,"#Archivo con datos clim'aticos #\n");
  fprintf(f_in3,"Localidad\t%s\n",file1);
  fprintf(f_in3,"Longitud\t%.2f\n",Lon);
  fprintf(f_in3,"Latitud\t%.2f\n",Lat);
  fprintf(f_in3,"Altura_msn\t%.0f\n",msnm);


  for (i = 0; i < 12; ++i) {
    if (i==0)  fprintf(f_in3,"mes\tEnero\n");
    if (i==1)  fprintf(f_in3,"mes\tFebrero\n");
    if (i==2)  fprintf(f_in3,"mes\tMarzo\n");
    if (i==3)  fprintf(f_in3,"mes\tAbril\n");
    if (i==4)  fprintf(f_in3,"mes\tMayo\n");
    if (i==5)  fprintf(f_in3,"mes\tJunio\n");
    if (i==6)  fprintf(f_in3,"mes\tJulio\n");
    if (i==7)  fprintf(f_in3,"mes\tAgosto\n");
    if (i==8)  fprintf(f_in3,"mes\tSeptiembre\n");
    if (i==9)  fprintf(f_in3,"mes\tOctubre\n");
    if (i==10) fprintf(f_in3,"mes\tNoviembre\n");
    if (i==11) fprintf(f_in3,"mes\tDiciembre\n");
    //printf("Ig= %f\t",i,Igmax[i]);
    fprintf(f_in3,"tTmax\t%.2f\n",tTmax[i]);
    fprintf(f_in3,"Tmax\t%.2f\n",Tmax[i]);
    fprintf(f_in3,"Tmin\t%.2f\n",Tmin[i]);
    fprintf(f_in3,"Imax\t%.2f\n",Ibmax[i]+Idmax[i]);
  }















}
 
   


