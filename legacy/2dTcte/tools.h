/************** tools.h  *************/



double Tsout(double **T,int nx,int ny) {

  double tmp;
  int i;
  
  tmp = 0.;

  for (i = 0; i < nx ; ++i) { 
    tmp += T[i][0];
  }

  return tmp/nx;
}

double Tsin(double **T,int nx,int ny) {

  double tmp;
  int i;
  
  tmp = 0.;

  for (i = 0; i < nx ; ++i) { 
    tmp += T[i][ny-1];
  }

  return tmp/nx;
}


convective_coefficients(double *hi,double beta,double Ts,double Tint) {

  if (beta <= 45. ) { //Esta condici'on s'olo se aplica a un techo
    if (Tint>Ts) *hi = 9.4;
    else *hi = 6.6;
  }
  else *hi = 8.1; //Coeficiente convectivo interior para cualquier muro  

}

void draw_viguetabovedilla2hueca(int nx,int ny,int **NT,int i1,int j1,int i2,int j2) {
  int i,j;

  //esquina superior izquierda
  NT[0][0] = 1;
  //esquina superior derecha
  NT[nx-1][0] = 2;
  //esquina inferior izquierda
  NT[0][ny-1] = 3;
  //esquina inferior derecha
  NT[nx-1][ny-1] = 4;
  
  //laterales adiab'aticos
  for (j = 1; j <= ny-2; ++j) {
    NT[0][j] = 6;
    NT[nx-1][j] = 7;
  }

  for (i = 1; i <=nx-2; ++i) {
    //frontera exterior convectiva
    NT[i][0] = 5;
    //frontera interior convectiva
    NT[i][ny-1] = 8;
  }
  //nodos interiores
  for (i = 1; i <= nx-2; ++i)
    for (j = 1; j <= ny-2; ++j)
      NT[i][j] = 13;
  
  //upper wall interior convection
  j = j1-1;
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 9;
    
  
  //lower wall interior convection
  j = j2;
  for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 10;

  //left wall convective interior
  i = i1-1;
  for (j = j1; j < j2 ; ++j)
    NT[i][j] = 11;
  //left wall convective interior
  i = i2;
  for (j = j1; j < j2 ; ++j)
    NT[i][j] = 12;
  
   for (j = j1; j < j2 ; ++j)
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 0;
  
  

  /* //imprime esquema de materiales */
  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d,%d\t",i,j); */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */
  //imprime esquema de materiales


  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d\t",NT[i][j]); */
  /*   printf("\n"); */
  /* } */

  //getchar();

  
}

void draw_viguetabovedilla2rellena(int nx,int ny,int **NT,int i1,int j1,int i2,int j2) {
  int i,j;

  //esquina superior izquierda
  NT[0][0] = 1;
  //esquina superior derecha
  NT[nx-1][0] = 2;
  //esquina inferior izquierda
  NT[0][ny-1] = 3;
  //esquina inferior derecha
  NT[nx-1][ny-1] = 4;
  
  //laterales adiab'aticos
  for (j = 1; j <= ny-2; ++j) {
    NT[0][j] = 6;
    NT[nx-1][j] = 7;
  }

  for (i = 1; i <=nx-2; ++i) {
    //frontera exterior convectiva
    NT[i][0] = 5;
    //frontera interior convectiva
    NT[i][ny-1] = 8;
  }
  //nodos interiores
  for (i = 1; i <= nx-2; ++i)
    for (j = 1; j <= ny-2; ++j)
      NT[i][j] = 13;
  
  //upper wall interior convection
  j = j1-1;
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 13 ;
    
  
  //lower wall interior convection
  j = j2;
  for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 13;

  //left wall convective interior
  i = i1-1;
  for (j = j1; j < j2 ; ++j)
    NT[i][j] = 13;
  //left wall convective interior
  i = i2;
  for (j = j1; j < j2 ; ++j)
    NT[i][j] = 13;
  
   for (j = j1; j < j2 ; ++j)
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 14;  //para establecer que es relleno
   
  

  /* //imprime esquema de materiales */
  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d,%d\t",i,j); */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */
  //imprime esquema de materiales


  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d\t",NT[i][j]); */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */

  
}



void set_krhoc_viguetabovedillarellena(double **k,double **rhoc,double dx,double dy,int nx, int ny,
				       double L1,double L2,double L3,double L4,double L5,double L6,double L7,
				       double k1,double k2,double k3,double k4,double k5,double k6,double k7,double kr,
				       double rhoc1,double rhoc2,double rhoc3,double rhoc4,double rhoc5,double rhoc6,double rhoc7,double rhocr,
		      int **NT)
{
  int i,j;
  for (i = 0; i < nx; ++i) {
    for (j = 0; j < (L1/dy) ; ++j) {
      k[i][j] = k1;
      rhoc[i][j] = rhoc1;
    }
    for (; j <  (L1 + L2)/dy  ; ++j) {
      k[i][j] = k2;
      rhoc[i][j] = rhoc2;
    }
    for (; j <  (L1 + L2 + L3)/dy   ; ++j) {
      k[i][j] = k3;
      rhoc[i][j] = rhoc3;
    }
    for (; j <  (L1 + L2 + L3 + L4)/dy   ; ++j) {
      k[i][j] = k4;
      rhoc[i][j] = rhoc4;
    }
    
    for (; j <  (L1 + L2 + L3 + L4 + L5)/dy   ; ++j) {
      k[i][j] = k5;
      rhoc[i][j] = rhoc5;
    }
    
    for (; j <  (L1 + L2 + L3 + L4 + L5 + L6)/dy   ; ++j) {
      k[i][j] = k6;
      rhoc[i][j] = rhoc6;
    }

    for (; j <  (L1 + L2 + L3 + L4 + L5 + L6 + L7)/dy  ; ++j) {
      k[i][j] = k7;
      rhoc[i][j] = rhoc7;
    }
  }


  for (i = 0; i < nx; ++i)
    for (j=0; j < ny; ++j) 
      if (NT[i][j] == 14) {
	k[i][j] = kr;
	rhoc[i][j] = rhocr;
	NT[i][j] = 13;
      }

  for (j = ny-1, i = 1; i < nx-1; ++i) 
    NT[i][j] = 8;
  /*
  for (j = 0; j < ny; ++j ) {
    for (i = 0; i < nx; ++i)
      printf("%d\t",NT[i][j]);
    printf("\n");
  }

  getchar();
  */


  /* for (j=0; j < ny; ++j) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%.2f\t",k[i][j]); */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */

}




void  draw_viguetabovedillarellena(int nx,int ny,int **NT,int i1,int i2,int j1,int j2) {
  int i,j;

  //esquina superior izquierda
  NT[0][0] = 1;
  //esquina superior derecha
  NT[nx-1][0] = 2;
  //esquina inferior izquierda
  NT[0][ny-1] = 3;
  //esquina inferior derecha
  NT[nx-1][ny-1] = 4;
  
  //laterales adiab'aticos
  for (j = 1; j <= ny-2; ++j) {
    NT[0][j] = 6;
    NT[nx-1][j] = 7;
  }

  for (i = 1; i <=nx-2; ++i) {
    //frontera exterior convectiva
    NT[i][0] = 5;
    //frontera interior convectiva
    NT[i][ny-1] = 8;
  }
  //nodos interiores
  for (i = 1; i <= nx-2; ++i)
    for (j = 1; j <= ny-2; ++j)
      NT[i][j] = 13;
  
  //PARA ESTABLECER QUE ES RELLENO Y LUEGO NT[][] = 13
   for (j = j1; j <= j2 ; ++j)
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 14;  

  /* //imprime esquema de materiales */
  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d,%d\t",i,j); */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */
  //imprime esquema de materiales

   /*
  for (j = 0; j < ny; ++j ) {
    for (i = 0; i < nx; ++i)
      printf("%d\t",NT[i][j]);
    printf("\n");
  }

  getchar();
   */
  
}






void set_krhocrelleno(double **k,double **rhoc,double dx,double dy,int nx, int ny,
	  double L1,double L2,double L3,double L4,double L5,double L6,double L7,
		      double k1,double k2,double k3,double k4,double k5,double k6,double k7,double kr,
		      double rhoc1,double rhoc2,double rhoc3,double rhoc4,double rhoc5,double rhoc6,double rhoc7,double rhocr,
		      int **NT)
{
  int i,j;
  for (i = 0; i < nx; ++i) {
    for (j = 0; j < (L1/dy) ; ++j) {
      k[i][j] = k1;
      rhoc[i][j] = rhoc1;
    }
    for (; j <  (L1 + L2)/dy  ; ++j) {
      k[i][j] = k2;
      rhoc[i][j] = rhoc2;
    }
    for (; j <  (L1 + L2 + L3)/dy   ; ++j) {
      k[i][j] = k3;
      rhoc[i][j] = rhoc3;
    }
    for (; j <  (L1 + L2 + L3 + L4)/dy   ; ++j) {
      k[i][j] = k4;
      rhoc[i][j] = rhoc4;
    }
    
    for (; j <  (L1 + L2 + L3 + L4 + L5)/dy   ; ++j) {
      k[i][j] = k5;
      rhoc[i][j] = rhoc5;
    }
    
    for (; j <  (L1 + L2 + L3 + L4 + L5 + L6)/dy   ; ++j) {
      k[i][j] = k6;
      rhoc[i][j] = rhoc6;
    }

    for (; j <  (L1 + L2 + L3 + L4 + L5 + L6 + L7)/dy  ; ++j) {
      k[i][j] = k7;
      rhoc[i][j] = rhoc7;
    }
  }


  for (i = 0; i < nx; ++i)
    for (j=0; j < ny; ++j) 
      if (NT[i][j] == 14) {
	k[i][j] = kr;
	rhoc[i][j] = rhocr;
	NT[i][j] = 13;
      }

}









void animation_gnuplot(int animation,double X,double Y,double Tsamax,double Tsamin) {
  int i;

  FILE *f_in;
  char file_animation[180];
  sprintf(file_animation,"./dat/animation.gnp");
  // archivo datos temporales NOMBRE
  f_in = fopen(file_animation,"w");
  printf(file_animation);
  printf("\n");
  
  
  fprintf(f_in,"set pm3d at b\n");
  fprintf(f_in,"unset surface\n");
  fprintf(f_in,"set view map\n");
  fprintf(f_in,"set size ratio %f\n",X/Y);
  fprintf(f_in,"set cbrange [%f:%f]\n",Tsamin,Tsamax);
  
  
  fprintf(f_in,"set term gif animate delay 50\n");
  fprintf(f_in,"set out 'animation.gif'\n");
  for (i = 0; i < animation; ++i) {
    fprintf(f_in,"set title '%.2f h'\n",(double) i/144.*24.);
    fprintf(f_in,"sp [][%f:0]'./dat/matrix_%d.dat' using 1:2:3 t ''\n",Y,i);
  }
  
  fprintf(f_in,"set term X'\n");
  fclose(f_in);

}
void field_temperature(double **T,double Thueco,int nx,int ny,int **NT,int animation,double X,double Y) {
  int i,j;
  FILE *f_in;
  char file_matrix[180];

  sprintf(file_matrix,"./dat/matrix_%d.dat",animation);
  // archivo datos temporales NOMBRE
  f_in = fopen(file_matrix,"w");
  printf(file_matrix);
  printf("\n");
     

  
  for (j = 0; j < ny; ++j )  {
    for (i = 0; i < nx; ++i) {
      if (NT[i][j]!=0) 
	fprintf(f_in,"%f\t%f\t%f\n",X*i/nx,Y*j/ny,T[i][j]);
      else 
	fprintf(f_in,"%f\t%f\t%f\n",X*i/nx,Y*j/ny,Thueco);
    }
    fprintf(f_in,"\n");
  }
   
  fclose(f_in);
}

void draw_viguetabovedilla2(int nx,int ny,int **NT,int i1,int j1,int i2,int j2) {
  int i,j;

  //esquina superior izquierda
  NT[0][0] = 1;
  //esquina superior derecha
  NT[nx-1][0] = 2;
  //esquina inferior izquierda
  NT[0][ny-1] = 3;
  //esquina inferior derecha
  NT[nx-1][ny-1] = 4;
  
  //laterales adiab'aticos
  for (j = 1; j <= ny-2; ++j) {
    NT[0][j] = 6;
    NT[nx-1][j] = 7;
  }

  for (i = 1; i <=nx-2; ++i) {
    //frontera exterior convectiva
    NT[i][0] = 5;
    //frontera interior convectiva
    NT[i][ny-1] = 8;
  }
  //nodos interiores
  for (i = 1; i <= nx-2; ++i)
    for (j = 1; j <= ny-2; ++j)
      NT[i][j] = 13;
  
  //upper wall interior convection
  j = j1-1;
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 9;
    
  
  //lower wall interior convection
  j = j2;
  for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 10;

  //left wall convective interior
  i = i1-1;
  for (j = j1; j < j2 ; ++j)
    NT[i][j] = 11;
  //left wall convective interior
  i = i2;
  for (j = j1; j < j2 ; ++j)
    NT[i][j] = 12;
  
   for (j = j1; j < j2 ; ++j)
    for ( i = i1; i < i2 ; ++i)
      NT[i][j] = 0;
  
  

  /* //imprime esquema de materiales */
  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d,%d\t",i,j); */
  /*   printf("\n"); */
  /* } */
  /* getchar(); */
  //imprime esquema de materiales


  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d\t",NT[i][j]); */
  /*   printf("\n"); */
  /* } */

  //getchar();

  
}
void interchange(double **T1, double **T2,int nx,int ny) {
  int i,j;

    for (i = 0; i < nx; ++i)
      for (j = 0; j <ny; ++j)
	T1[i][j] = T2[i][j];
}

abrefile_indice(char *user,int month, int cont_herramienta,FILE **f_in,
		double Tsemax,double Tsemin,double Tsintmax,double Tsintmin,
		double Tint,double dt,double Qcool,double Qheat,double X) {
  char file_indices[180];
  //double contador;

  //contador = 86400/ dt;

  sprintf(file_indices,"./dat/indice_%s_%d_%d.csv",user,month,cont_herramienta);
  // archivo datos temporales NOMBRE
  *f_in = fopen(file_indices,"w");
  printf(file_indices);

  fprintf(*f_in,"%.3f\t",Qcool/3600./X);
  fprintf(*f_in,"%.3f\t",Qheat/3600./X);
  fprintf(*f_in,"%.3f\t",(Qcool+Qheat)/3600./X);
  fprintf(*f_in,"%.2f\t",(Tsintmax-Tsintmin)/(Tsemax-Tsemin));
  fprintf(*f_in,"%.1f\t",Tint);
  fclose(*f_in);
  printf("\n");
}

void abrefile(char *user,int month,int cont_herramienta,FILE **f_in) {
  char s_file_capas[140];
  

  sprintf(s_file_capas,"./dat/%s_%d_%d.csv",user,month,cont_herramienta);
  printf(s_file_capas);
  printf("\n");


  // archivo datos temporales NOMBRE
  *f_in = fopen(s_file_capas,"w");
  fprintf(*f_in,"Is\tTsa\tTa\tT_paredext\tTparedint\tTint\tT_n\tDeltaT_n\n");
  fprintf(*f_in,"[W/m2]\t[oC]\t[oC]\t[oC]\t[oC]\t[oC]\t[oC]\t[oC]\n");

}
void set_nodetype(int nx, int ny,int **NT)
{

  int i,j;
 
  //esquina superior izquierda
  NT[0][0] = 1;
  //esquina superior derecha
  NT[nx-1][0] = 2;
  //esquina inferior izquierda
  NT[0][ny-1] = 3;
  //esquina inferior derecha
  NT[nx-1][ny-1] = 4;
  
  //laterales adiab'aticos
  for (j = 1; j <= ny-2; ++j) {
    NT[0][j] = 6;
    NT[nx-1][j] = 7;
  }

  for (i = 1; i <=nx-2; ++i) {
    //frontera exterior convectiva
    NT[i][0] = 5;
    //frontera interior convectiva
    NT[i][ny-1] = 8;
    }

  //nodos interiores
  for (i = 1; i <= nx-2; ++i)
    for (j = 1; j <= ny-2; ++j)
      NT[i][j] = 13;

  //imprime esquema de materiales
  /* for (j = 0; j < ny; ++j ) { */
  /*   for (i = 0; i < nx; ++i) */
  /*     printf("%d\t",NT[i][j]); */
  /*   printf("\n"); */
  /* } */

  /* getchar(); */


}
void set_krhoc(double **k,double **rhoc,double dx,double dy,int nx, int ny,
	  double L1,double L2,double L3,double L4,double L5,double L6,double L7,
	  double k1,double k2,double k3,double k4,double k5,double k6,double k7,
	  double rhoc1,double rhoc2,double rhoc3,double rhoc4,double rhoc5,double rhoc6,double rhoc7)
{
  int i,j;
  for (i = 0; i < nx; ++i) {
    for (j = 0; j < (L1/dy) ; ++j) {
      k[i][j] = k1;
      rhoc[i][j] = rhoc1;
    }
    for (; j <  (L1 + L2)/dy  ; ++j) {
      k[i][j] = k2;
      rhoc[i][j] = rhoc2;
    }
    for (; j <  (L1 + L2 + L3)/dy   ; ++j) {
      k[i][j] = k3;
      rhoc[i][j] = rhoc3;
    }
    for (; j <  (L1 + L2 + L3 + L4)/dy   ; ++j) {
      k[i][j] = k4;
      rhoc[i][j] = rhoc4;
    }
    
    for (; j <  (L1 + L2 + L3 + L4 + L5)/dy   ; ++j) {
      k[i][j] = k5;
      rhoc[i][j] = rhoc5;
    }
    
    for (; j <  (L1 + L2 + L3 + L4 + L5 + L6)/dy   ; ++j) {
      k[i][j] = k6;
      rhoc[i][j] = rhoc6;
    }

    for (; j <  (L1 + L2 + L3 + L4 + L5 + L6 + L7)/dy  ; ++j) {
      k[i][j] = k7;
      rhoc[i][j] = rhoc7;
    }
  }


}
void initial_conditions (double **T,double **Tn,double **Terror,double **To,int nx,int ny,double Tint) {
  int i,j;
  for (i = 0; i < nx; ++i ) 
    for (j = 0; j < ny; ++j) { 
    T[i][j] = Tint;
    Tn[i][j] = Tint;
    Terror[i][j] = Tint;
    To[i][j] = Tint;
  }

}


void Ta_Tc_DtaT(double Ig,double ho,double Tmax,double Tmin,double Hi,
      double Ho,double *Ta,double *Tn,double *DtaT) {
  
  double tm2,t,y,pi;
  
  *Ta = 0.;
  pi = acos(-1.);
  y = 0.;

  for (t = 0; t <= 86400; ++t) {
    if (t/3600. <= Ho) 
      y = (cos(pi*(Ho-t/3600.)/(24.+Ho-Hi) )+ 1.)/2. ;
    
    if ( (t/3600.>Ho) && (t/3600.<=Hi) ) 
      y = (cos(pi*(t/3600.-Ho)/ (Hi  - Ho) ) +1.)/2.;
    
    if (t/3600. > Hi) 
      y = (cos(pi*(24. + Ho-t/3600.)/(24. + Ho - Hi) ) + 1.)/2. ;
    
    tm2= y*Tmin + (1.-y)*Tmax;
    *Ta += tm2;
  }
  
  *Ta /= 86400.;
  *Tn = 0.54*(*Ta) + 13.5;  
  
  double Delta,tmp2;
 
  
  Delta = ( Tmax - Tmin );
  if (Delta < 13.) tmp2 = 2.5/2.;
  if (Delta >=13. && Delta < 16.) tmp2 = 3.0/2.;
  if (Delta >=16. && Delta < 19.) tmp2 = 3.5/2.;
  if (Delta >=19. && Delta < 24.) tmp2 = 4.0/2.;
  if (Delta >=24. && Delta < 28.) tmp2 = 4.5/2.;
  if (Delta >=28. && Delta < 33.) tmp2 = 5./2.;
  if (Delta >=33. && Delta < 38.) tmp2 = 5.5/2.;
  if (Delta >=38. && Delta < 45.) tmp2 = 6./2.;
  if (Delta >=45. && Delta < 52.) tmp2 = 6.5/2.;
  if (Delta >=52.) tmp2 = 7./2.;
  *DtaT = tmp2;    
}



 
double time_evolution_Tsa(double *Ta,double t,double Ig,
			  double a,double ho,double Tmax,double Tmin,double Hi,
			  double Ho,double tau,double delta,double Ib,double Id,double gamma,
			  double phi,double beta,double *Is) {
  double tm2,tm3;
  double y,pi,CF;
  double omega,X,theta,thetaz;
  double I; 
  double Idtheta,Ibtheta;
  
  pi = acos(-1.);
  X = pi/180.;
  Ibtheta = Idtheta = 0.;


  // Varies lineraly the CF factor with slope of the wall, >90 CF=0.
  CF = 3.9*(1. - beta/90. );
  if (beta >90.) CF =0. ;
  
    
  if (t/3600. <= Ho) 
    y = (cos(pi*(Ho-t/3600.)/(24.+Ho-Hi) )+ 1.)/2. ;
  
  if ( (t/3600.>Ho) && (t/3600.<=Hi) ) 
    y = (cos(pi*(t/3600.-Ho)/ (Hi  - Ho) ) +1.)/2.;
  
  if (t/3600. > Hi) 
    y = (cos(pi*(24. + Ho-t/3600.)/(24. + Ho - Hi) ) + 1.)/2. ;
  
  tm2 = y*Tmin + (1.-y)*Tmax;
  *Ta =tm2;
  
  if ( sin(2.*M_PI*( t/3600./tau - Ho/tau))>0.) {
    omega = t/3600*15. -180.;
    //printf("omega = %f\n",omega);
    thetaz = cos(phi*X)*cos(delta*X)*cos(omega*X)
      + sin(phi*X)*sin(delta*X);
    theta =  sin(delta*X)*sin(phi*X)*cos(beta*X) 
      - sin (delta*X)*cos(phi*X)*sin(beta*X)*cos(gamma*X)
      + cos(delta*X)*cos(phi*X)*cos(beta*X)*cos(omega*X)
      + cos(delta*X)*sin(phi*X)*sin(beta*X)*cos(gamma*X)*cos(omega*X)
      + cos(delta*X)*sin(beta*X)*sin(gamma*X)*sin(omega*X);/*
      if (fmod(t,600.)<dtat) 	
      printf("%f\t%f\n",(omega+180.)/15.,acos(theta)/X);//getchar();*/
      theta = acos(theta);
      thetaz = acos(thetaz);
      Ibtheta = Ib* sin(2.*M_PI*(t/tau/3600.-Ho/tau))/cos(thetaz)*cos(theta);
      Idtheta = Id*sin(2.*M_PI*(t/tau/3600.-Ho/tau))*(1.-beta/180.);
      if (Ibtheta <0.) Ibtheta = 0.;
      //printf("relaci'on = %f\n",1.-beta/180.);getchar();
      I = Ig*sin(2.*M_PI*(t/tau/3600.-Ho/tau));
      /*if (fmod(t,60.)<dtat) 	{
	printf("%f\t%f\t%f\t%f\n",(omega+180.)/15.,I,Ibtheta,Idtheta);
	getchar();
	}*/
  }
  else
    I =0.;
  
  //printf("Mpi=%f\n",M_PI);getchar(); 
  tm3 = tm2 + a*(Idtheta+Ibtheta)/ho - CF;
  *Is = Idtheta+Ibtheta;

  /*
  printf("Radiaci'on = %f\n",Idtheta+Ibtheta);
  printf("Tm3 = %f\tCF=%f\tt=%f\n",tm3,CF,t);
  */

  return tm3;

}

void max_min(double *Tsemax,double *Tsemin,double *Tsintmax,double *Tsintmin,
        double t,double Tint,double **T,int nx,int ny,double *Tsi) {

  double Tavese,Tavesint;
  int i;
  Tavese = Tavesint = 0.;
  
  for (i = 0; i < nx; ++i) {
    Tavese   += T[i][0];
    Tavesint += T[i][ny-1];
  }

  Tavese   = Tavese   / (nx-1);
  Tavesint = Tavesint / (nx-1); 

 

  if (Tavese > *Tsemax)
    *Tsemax  = Tavese;
  if (Tavese  < *Tsemin)
    *Tsemin  = Tavese;

  *Tsi = Tavesint;

  if (Tavesint < *Tsintmin)
    *Tsintmin = Tavesint;

  if (Tavesint > *Tsintmax)
    *Tsintmax = Tavesint;


}


double calculate_error (double **T,double **Terror,int nx,int ny) 
{
  int i,j;
  double error;
  error = 0.;
  
  for (i=0; i < nx; ++i) 
    for (j = 0; j < ny; ++j)
      error += fabs(T[i][j] - Terror[i][j]);

  return error/nx/ny;
  
}


void calculate_coefficients(double **a, double **b, double **c, double **d, double dt, double dx, double dy,
			    double **k, int nx, int ny, int **NT, double **rhoc, double **T, 
			    double Tsa,double ho,double Tint,double hi,
			    double **To,double *Thueco,double *hhueco,
			    double Qud,double Qul,double Qur,
			    double Qru,double Qrd,double Qrl,
			    double Qdu,double Qdr,double Qdl,
			    double Qlu,double Qlr,double Qld)
{



   int i,j;
   double an,as,ae,aw,apo;
   double Ti,Th,hh;
   Ti = Tint;
   Th = *Thueco;
   hh = *hhueco;

   for (i = 0; i < nx; ++i) 
     for (j = 0; j < ny; ++j) {
       //printf("NT[%d][%d] = %d\n",i,j,NT[i][j]);getchar();
       switch (NT[i][j]) {
       case 1: //superior left node with convective upper and adiabatic left
	 /* printf("case 1 in\n");getchar(); */
	 /* printf("[%d][%d] = \n",i,j);getchar(); */
	 apo = rhoc[i][j] * dx * dy / dt;
	 /* printf("case 1 middle\n");getchar(); */
	 an = 0.;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = 0.;
	 
	 //	 printf("case 1 middle\n");getchar();
	 a[i][j] =  apo + ae + ho*dx + as;  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  ho*dx* Tsa + as * T[i][j+1] + apo * To[i][j];
	 //printf("case 1\n");getchar();
	 break;
       case 2: // superior right node with convective upper and adiabatic right
	 //printf("case 2 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = 0.;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = 0;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + ho*dx + as + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  ho*dx* Tsa + as * T[i][j+1] + apo * To[i][j];
	 //printf("case 2\n");getchar();
	 break;
       case 3: // left lower node adiabatic left convective down
	 //printf("case 3 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = 0.;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = 0.;
	 
	 a[i][j] =  apo + an + hi*dx + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + hi*dx * Ti + apo * To[i][j];

	 //	  printf("case 3\n");getchar();
	 break;
       case 4:  
	 //	 printf("case 4 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = 0.;
	 ae = 0.;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + an + hi*dx + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + hi*dx * Ti + apo * To[i][j]; 
	 
	 //	 printf("case 4\n");getchar();
	 break;
       case 5://convective boundary in the upper side
	 //printf("case 5 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = 0.;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + ho*dx + as + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  ho*dx* Tsa + as * T[i][j+1] + apo * To[i][j];
	 
	 //printf("case 5\n");getchar();
	 break;      
       case 6:// adiabatic boundary at left side
	 //printf("case 6 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = 0.;
	 
	 a[i][j] =  apo + an + as + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + as * T[i][j+1] + apo * To[i][j];
	 
	 //printf("case 6\n");getchar();
	 break;
       case 7: //adiabatic boundary at right side
	 //printf("case 7 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = 0.;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + an + as + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + as * T[i][j+1] + apo * To[i][j];
	 
	 //printf("case 7\n");getchar();
	 break;
       case 8://convective boundary in the lower limit
	 //printf("case 8 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = 0.;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + an + hi*dx + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + hi* dx * Ti + apo * To[i][j];
	 //printf("case 8\n");getchar();
	 break;

       case 9:  //upper hole convection
	 //printf("case 9 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = 0.;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + an + hh*dx + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + hh* dx * Th + apo * To[i][j] - Qur - Qud - Qul;
	 
	 //printf("case 9\n");getchar();
	 break;

       case 10://lower hole convection
	 //printf("case 10 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = 0.;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + hh*dx + as + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  hh*dx* Th + as * T[i][j+1] + apo * To[i][j] - Qdl - Qdu - Qdr;
	 break;      


       case 11:  //left hole convective wall

	 //printf("case 11 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = 0.;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + an + as + hh*dy + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + as * T[i][j+1] + apo * To[i][j]  + hh * dy * Th - Qlu -  Qlr - Qld;
	 //printf("d[%d][%d] = %f\n",i,j,d[i][j]);//getchar();
	 //printf("case 11\n");getchar();
	 break;


       case 12:  //right hole convective wall
	 // printf("case 12 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = 0.;
	
	 a[i][j] =  apo + an + as + ae + hh*dy;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + as * T[i][j+1] + apo * To[i][j] + hh * dy *Th - Qrd - Qrl - Qru;

	 //printf("case 12\n");getchar();
	 break;

       case 13:  //interior nodes

	 //printf("case 13 in\n");getchar();
	 apo = rhoc[i][j] * dx * dy / dt;
	 an = (2.* k[i][j-1]*k[i][j]/(k[i][j]+k[i][j-1]) ) * dx / dy;
	 as = (2.* k[i][j+1]*k[i][j]/(k[i][j]+k[i][j+1]) ) * dx / dy;
	 ae = (2.* k[i+1][j] * k[i][j]/(k[i+1][j] + k[i][j])) * dy / dx;
	 aw = (2.* k[i-1][j] * k[i][j]/(k[i-1][j] + k[i][j])) * dy / dx;
	 
	 a[i][j] =  apo + an + as + ae + aw;                  
	 b[i][j] = ae;           
	 c[i][j] = aw;          
	 d[i][j] =  an* T[i][j-1] + as * T[i][j+1] + apo * To[i][j];
	 //printf("case 13\n");getchar();
	 break;
      
       case 0:
	 break;
       }
     }

   //printf("despu'es de los case\n");getchar();
						
   for (i = 0; i < nx; ++i)
     for (j = 0; j<ny;++j)
       if (NT[i][j] == 0)
	 {
	 a[i][j] =  1;                  
	 b[i][j] = 0.;           
	 c[i][j] = 0.;          
	 d[i][j] =  Th;
	 //printf("Th = %f\n",Th);getchar();
	 }
 
 }


  
  
void solve_PQ (double **a,double **b,double **c,double **d,double *P,double *Q,double **Tn,double **T,int nx,int ny,
	       double Tint,double hi,double rhoair,double cair,double La,double *Qin,double dt,
	       double dx,double dy,double **k,double **rhoc,int **NT,double Tsa,double ho,double **To,double X,double t,
	       double *Thueco, double *hh, 
	       int i1, int j1, int i2, int j2,double a11,double a21,double a12,double e22,
	       double E,double *Qrup,double *Qrdown,double *Qcool,double *Qheat,
	       double tipo,double beta, double *Tarriba,double *Tabajo,double *Nur) {

   int i,j;
  double Ti,error,Qh,Qhueco,Th;
  double Tup,Tdown,Tupave,Tdownave,Tleftave,Trightave;
  double cont_up,cont_down,cont_left;
  double H,sigma,F12;
  double Qrtmp,dot1,dot2; 
  double dot11,dot22;
  double Ra,gr,Beta,nu,alphaair,kair;
  double Qupdown,Qupright,Qupleft;
  double Qrightup,Qrightdown,Qrightleft;
  double Qdownup,Qdownright,Qdownleft;
  double Qleftup,Qleftright,Qleftdown;
  double Fur,Fud,Ful;
  double Frd,Frl,Fru;
  double Fdl,Fdu,Fdr;
  double Flu,Flr,Fld;
  double h,l,wi,wj;

  gr = 9.81;
  Beta = 1./300.;
  nu = 1.11e-5;  
  kair = 0.0262;
  alphaair = kair/rhoair/cair;
  
  //printf("alphaair in = %e\n",alphaair);getchar();
  Ti = Tint;
  Th = *Thueco;
  error = 0;
  sigma = 5.6704e-8;
  h = e22;
  l = a21;
  wj = h;  wi = l;
  Fur = 0.5* (1. + wj/wi - pow(1.+(wj*wj/wi/wi),0.5));
  Ful = Fur;  Fud = 1. - 2.*Fur;
  wj = l;  wi = h;
  Fru =  0.5* (1. + wj/wi - pow(1.+(wj*wj/wi/wi),0.5));
  Frd = Fru;  Frl = 1. - 2.*Fru;
  Fdl = Ful;  Fdr = Fur;  Fdu = Fud;
  Flu = Fru;  Flr = Frl;  Fld = Frd;


  Qh = Qhueco = 0.;
  do  {
    Tupave = Tdownave = Trightave = Tleftave = 0.; 
    error = cont_up = cont_left = 0.;
    *hh = 0.4005*pow(fabs(Tupave-Tdownave),0.3033)/pow(e22,0.0901);
    if ( tipo==1)  {   //bloque simetrico hueco aire
      for ( i = i1,j = j1-1; i < i2 ; ++i)  {   //calculate average up and down temperatures
	Tupave += T[i][j1-1];
	Tdownave += T[i][j2];
	++cont_up;
      }
      for (j = j1; j < j2; ++j) {
	Tleftave += T[i1-1][j];
	Trightave += T[i2][j];
	++cont_left;
      }
      Tupave    =     Tupave/cont_up;
      Tdownave  =   Tdownave/cont_up;
      Tleftave  =  Tleftave/cont_left;
      Trightave = Trightave/cont_left; 
      *Tarriba = Tupave;
      *Tabajo  = Tdownave;
      
      
      double C = 1.0;
      Qupdown     =  dx*E*sigma*(pow(Tupave+273.15,4.)-pow(Tdownave+273.15,4.))*Fud;
      Qupleft     =  dx*E*sigma*(pow(Tupave+273.15,4.)-pow(Tleftave+273.15,4.))*Ful*C;
      Qupright    =  dx*E*sigma*(pow(Tupave+273.15,4.)-pow(Trightave+273.15,4.))*Fur*C;
      
      Qrightup    =  dy*E*sigma*(pow(Trightave+273.15,4.)-pow(Tupave+273.15,4.))*Fru*C;
      Qrightdown  =  dy*E*sigma*(pow(Trightave+273.15,4.)-pow(Tdownave+273.15,4.))*Frd*C;
      Qrightleft  =  dy*E*sigma*(pow(Trightave+273.15,4.)-pow(Tleftave+273.15,4.))*Frl*C;
      
      Qdownup     =  dx*E*sigma*(pow(Tdownave+273.15,4.)-pow(Tupave+273.15,4.))*Fdu;
      Qdownright  =  dx*E*sigma*(pow(Tdownave+273.15,4.)-pow(Trightave+273.15,4.))*Fdr*C;
      Qdownleft   =  dx*E*sigma*(pow(Tdownave+273.15,4.)-pow(Tleftave+273.15,4.))*Fdl;
      
      Qleftup     =  dy*E*sigma*(pow(Tleftave+273.15,4.)-pow(Tupave+273.15,4.))*Flu*C;
      Qleftright  =  dy*E*sigma*(pow(Tleftave+273.15,4.)-pow(Trightave+273.15,4.))*Flr*C;
      Qleftdown   =  dy*E*sigma*(pow(Tleftave+273.15,4.)-pow(Tdownave+273.15,4.))*Fld*C;
          
      *Nur = sigma*(pow(Tupave+273.15,4.)-pow(Tdownave+273.15,4.))*Fud*E*h/kair/(Tupave-Tdownave);
      double Qradi = 0.;
      Qradi = -Qdownup/dx - Qdownleft/dx - Qdownright/dx; 
      *Nur = Qradi/kair/(Tupave-Tdownave)*h; 
      if (beta ==90.) { //es un muro 
	*hh = 0.4005*pow(fabs(Tupave-Tdownave),0.3033)/pow(e22,0.0901);
      }
      if (beta ==0) { // un techo 
	Ra = gr*Beta*(Tdownave-Tupave)*pow(e22,3.)/nu/alphaair;
	dot11 = 1. - 1708./Ra;
	dot22 = pow(Ra/5830.,1./3.)-1.;
	if ( dot11<0.) dot11 = 0.;
	if  (dot22<0.) dot22 = 0.;
	*hh = kair/e22*(1. + 1.44*dot11 + dot22);
	if (Tdownave<Tupave) {
	  *hh = kair/e22;
	}
      }
    }
    if (tipo ==2) {
      *Qrup   = 0.;
      *Qrdown = 0.;
      Qrtmp = 0.; 
      Tupave = Tdownave = 0.;
      *hh = 1.;
    }
    
    
    calculate_coefficients(a,b,c,d,dt,dx,dy,k,nx,ny,NT,rhoc,T,Tsa,ho,Tint,hi,To,Thueco,hh,
			   Qupdown,Qupleft,Qupright,
			   Qrightup,Qrightdown,Qrightleft,
			   Qdownup,Qdownright,Qdownleft,
			   Qleftup,Qleftright,Qleftdown);
    for (j = 0; j < ny ; ++j) {
      P[0] = b[0][j]/a[0][j];
      Q[0] = d[0][j]/a[0][j];
      for (i = 1 ; i < nx; ++i) {
	P[i] = b[i][j] / ( a[i][j] - c[i][j] * P[i-1] ); 
	Q[i] = ( d[i][j] + c[i][j] * Q[i-1] ) / (a[i][j] - c[i][j] * P[i-1] );
      }
      Tn[nx-1][j] = Q[nx-1] ;
      for (i = nx-2; i > -1; --i)  {
	Tn[i][j] = P[i] * Tn[i+1][j] + Q[i];
      }
    }
    
    //evaluate error
    for (j = 0; j < ny; ++j) 
      for (i = 0; i < nx; ++i) 
	error += (T[i][j] - Tn[i][j] )/T[i][j]/nx/ny; 
  
    interchange(T,Tn,nx,ny);
    
  }  while (fabs(error) >1e-10 );
 

  //Evaluate Qin and total heat transfer by convection
  for (i = 0; i < nx; ++i) {
    if (T[i][ny-1] > Tint)
      *Qcool += hi*dt*dx*(T[i][ny-1]-Tint);
    if (T[i][ny-1] < Tint)
      *Qheat += hi*dt*dx*(Tint - T[i][ny-1]);
  } 
  
  cont_up = cont_down = 0.;
  Tup = Tdown = 0.;

  //calculate average temperatures
  //upper wall 
  j = j1-1;
  for ( i = i1; i < i2 ; ++i)  {
    Tup += T[i][j];
    ++cont_up;
  }
  //lower wall
  j = j2;
  for ( i = i1; i < i2 ; ++i) {
    Tdown += T[i][j];
    ++cont_down;
  }
  Tup = Tup/cont_up;
  Tdown = Tdown/cont_down;
  


  //calculate heat flux  
  //upper wall 
  j = j1-1;
  for ( i = i1; i < i2 ; ++i)  
    Qhueco += *hh * dx * (T[i][j] - Th);
  //lower wall
  j = j2;
  for ( i = i1; i < i2 ; ++i)
    Qhueco += *hh * dx * (T[i][j] - Th);
  //left wall convective interior
  i = i1-1;
  for (j = j1; j < j2 ; ++j)
    Qhueco += *hh * dy * (T[i][j] - Th);
  //left wall convective interior
  i = i2;
  for (j = j1; j < j2 ; ++j)
    Qhueco += *hh * dy * (T[i][j] - Th);
  
     /* if (T[i][ny-1] > Tinn)  */
     /*    *Qin += hi*dt*dx*(T[i][ny-1]-Tinn);  */
  
  // printf("Tw[%d][%d] = %f\t Tw[%d][%d] = %f\n",0,0,Tn[0][0],0,ny-1,Tn[0][ny-1]); //getchar();
  //*Tint = ( Qh + (rhoair*cair*La*X/dt)*Ti)*dt/(rhoair*cair*La*X);
  *Thueco = ( Qhueco + (rhoair*cair*a21*e22/dt)*Th)*dt/(rhoair*cair*a21*e22);
  //printf("Tint = %f\n",Tinn);

  //*Tintaverage += *Tint;
  
}



database_begin(char *user,int cont_herramienta,int month, char *baseDatos) {

  int pid;
  PGconn *conn;
  PGresult *activa;
  char orden[250];

char pghost[] = "localhost";
char pgport[] = "5432";
char pguser[] = "postgres";
char pgpass[] = "";

  pid = getpid();  
  //conn = PQsetdbLogin("localhost","5432",NULL,NULL,"usuarios","postgres","");
  //conn = PQsetdbLogin("localhost","5432","NULL","NULL",baseDatos,"postgres","");
  conn = PQsetdbLogin(pghost,pgport,NULL,NULL,baseDatos,pguser,pgpass);
  //sprintf(conn,"PQsetdbLogin('localhost','5432','NULL','NULL','%s','postgres','')",baseDatos);
  sprintf(orden,"update anuales set pid='%d',activo='t' where login='%s' and sc='%d' and mes='%d'",pid,user,cont_herramienta,month);

  if (PQstatus(conn) != CONNECTION_BAD) {
  	activa = PQexec(conn, orden);
  	//printf("NOMBRE DE LA BASE DE DATOS DESPUES DE EJECUTAR: %s\n",baseDatos);
	  
  }

  else{
  	printf("NO se conecto");
	fprintf(stderr, "Connection to database failed: %s\n",
            PQerrorMessage(conn));
        
  }

  PQclear(activa);
  PQfinish(conn);


}


database_end(char *user,int cont_herramienta,int month, char *baseDatos) {

  int pid;
  PGconn *conn;
  PGresult *res;
  PGresult *activa;
  PGresult *desactiva;
  char orden2[230];
  char orden3[230]; 

  char pghost[] = "localhost";
  char pgport[] = "5432";
  char pguser[] = "postgres";
  char pgpass[] = "";
 
  
  pid = getpid();
  
  //conn = PQsetdbLogin("localhost","5432",NULL,NULL,"usuarios","postgres","");
  //conn = PQsetdbLogin("localhost","5432","NULL","NULL",baseDatos,"postgres","");
  conn = PQsetdbLogin(pghost,pgport,NULL,NULL,baseDatos,pguser,pgpass);
  //sprintf(conn,"PQsetdbLogin('localhost','5432','NULL','NULL','%s','postgres','')",baseDatos);
  sprintf(orden3,"update sesiones set estado=estado+1 where login='%s'",user);
  //res = PQexec(conn, orden3);
  sprintf(orden2,"update anuales set pid='%d',activo='f' where login='%s' and sc='%d' and mes='%d'",pid,user,cont_herramienta,month);
  //desactiva = PQexec(conn, orden2);


 if (PQstatus(conn) != CONNECTION_BAD) {
  	res = PQexec(conn, orden3);
        desactiva = PQexec(conn, orden2);  	
  }

  else{
  	printf("NO se conecto");
	fprintf(stderr, "Connection to database failed: %s\n",
            PQerrorMessage(conn));
        
  }


  
  PQclear(desactiva);
  PQclear(res);
  PQfinish(conn);
}
