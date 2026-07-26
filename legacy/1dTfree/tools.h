/**************** tools.h  *********/
set_krhoc(double *k,double *rhoc,double dx,
	  double L1,double L2,double L3,double L4,double L5,double L6,double L7,
	  double k1,double k2,double k3,double k4,double k5,double k6,double k7,
	  double rhoc1,double rhoc2,double rhoc3,double rhoc4,double rhoc5,double rhoc6,double rhoc7)
{
  int i;
  for (i = 0; i < (L1/dx) ; ++i) {
    k[i] = k1;
    rhoc[i] = rhoc1;
    //printf("i = %d\n",i);getchar();
  }
  
  for (; i <  (L1 + L2)/dx  ; ++i) {
    k[i] = k2;
    rhoc[i] = rhoc2;
    //printf("i = %d\n",i);getchar();
  }
    
  
  for (; i <  (L1 + L2 + L3)/dx   ; ++i) {
    k[i] = k3;
    rhoc[i] = rhoc3;
    //printf("i = %d\n",i);getchar();
  }
  for (; i <  (L1 + L2 + L3 + L4)/dx   ; ++i) {
    k[i] = k4;
    rhoc[i] = rhoc4;
    //printf("i = %d\n",i);getchar();
  }

  for (; i <  (L1 + L2 + L3 + L4 + L5)/dx   ; ++i) {
    k[i] = k5;
    rhoc[i] = rhoc5;
    //printf("i = %d\n",i);getchar();
  }

  for (; i <  (L1 + L2 + L3 + L4 + L5 + L6)/dx   ; ++i) {
    k[i] = k6;
    rhoc[i] = rhoc6;
    //printf("i = %d\n",i);getchar();
  }

  for (; i <  (L1 + L2 + L3 + L4 + L5 + L6 + L7)/dx   ; ++i) {
    k[i] = k7;
    rhoc[i] = rhoc7;
    //printf("i = %d\n",i);getchar();
    }
}

convective_coefficients(double *hi,double beta,double Ts,double Tint) {

  if (beta <= 45. ) { //Esta condici'on s'olo se aplica a un techo
    if (Tint>Ts) *hi = 9.4;
    else *hi = 6.6;
  }
  else *hi = 8.1; //Coeficiente convectivo interior para cualquier muro  

}


calculate_coefficients(double *a, double *b, double *c, double *d, double dt, double dx, 
		       double *k, int nx, double *rhoc, double *T, double To,double ho,
		       double Ti,double hi)
 {
   int i,j;
   
   


   b[0] = (2.*k[0]*k[1]) / (k[0]+k[1]) / dx;
   c[0] = 0.;
   d[0] = rhoc[0] * dx / dt * T[0] + ho*To;
   a[0] = rhoc[0] * dx / dt  + ho + b[0] ;
   
   for (i = 1; i <= nx-2; ++i)  {
     b[i] = ( 2.*k[i]  *k[i+1] ) / ( k[i]+k[i+1]) / dx; 
     c[i] = ( 2.*k[i-1]*k[i]   ) / ( k[i]+k[i-1]) / dx;
     d[i] = rhoc[i] * dx / dt * T[i] ;
     a[i] = rhoc[i] * dx / dt + b[i] + c[i] ;
     }

   i = nx-1;
   b[i] = 0.;
   c[i] = (2.*k[i-1]*k[i]) / (k[i]+k[i-1]) / dx;
   d[i] = rhoc[i] * dx / dt * T[i] + hi*Ti ;
   a[i] = rhoc[i] * dx / dt  +  c[i] + hi;

}



Ta_Tc_DtaT(double Ig,double ho,double Tmax,double Tmin,double Hi,
      double Ho,double *Ta,double *Tn,double *DtaT) {
  
  double tm2,tm3,t,y,pi;
  
  *Ta = 0.;
  pi = acos(-1.);
  
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

max_min(double *Tsamax,double *Tsamin,double *Tintmax,double *Tintmin,
	double *t_Tintmax, double *t_Tsamax,double t,double Tsa,double Tint) {

  
  if (Tsa  > *Tsamax)  { 
    *Tsamax  = Tsa;
    *t_Tsamax = t;
  }
  if (Tsa  < *Tsamin) {
    *Tsamin  = Tsa;
  }
  
  if (Tint < *Tintmin) {
    *Tintmin = Tint;
  }
  if (Tint > *Tintmax) {
    *Tintmax = Tint;
    *t_Tintmax = t;
  }


}


double time_evolution_Tsa(double *Ta,double t,double Ig,
			  double a,double ho,double Tmax,double Tmin,double Hi,
			  double Ho,double tau,double delta,double Ib,double Id,double gamma,
			  double phi,double beta,double *Is) {
  double tm2,tm3;
  double y,pi,CF;
  double omega,X,theta,thetaz;
  double Ibh,Ibp,I; 
  double Igtheta,Idtheta,Ibtheta;
  
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
      + cos(delta*X)*sin(beta*X)*sin(gamma*X)*sin(omega*X);
      theta = acos(theta);
      thetaz = acos(thetaz);
      Ibtheta = Ib* sin(2.*M_PI*(t/tau/3600.-Ho/tau))/cos(thetaz)*cos(theta);
      Idtheta = Id*sin(2.*M_PI*(t/tau/3600.-Ho/tau))*(1.-beta/180.);
      if (Ibtheta <0.) Ibtheta = 0.;
      //printf("relaci'on = %f\n",1.-beta/180.);getchar();
      I = Ig*sin(2.*M_PI*(t/tau/3600.-Ho/tau));
      
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
  /*  if (fmod(t,600.)<1.) 	{
    printf("%f\t%f\t%f\t%f\n",(omega+180.)/15.,I,Ibtheta,Idtheta);
    //getchar();
    }*/
  
  return tm3;

}

initial_conditions (double *T,double *Tn,int nx,double Tint) {
  int i;
  for (i = 0; i < nx; ++i ) {
    T[i] = Tint;
    Tn[i] = Tint + 1.;
  }
  //  printf("Tint = %f\n",Tint);getchar();
}
int calculate_error (double *To,double *Tf,int nx,double criterio) 
{
  int i,tmp;
  tmp = 0;
 
  for (i=0; i < nx; ++i)  
    if ( fabs(To[i] - Tf[i]) > criterio)
      ++tmp;

  return tmp;
  
}

int error_maxmin(double Tintmaxo,double Tintmino,double Tintmaxf,double Tintminf,double criterio) {

  int tmp;
  tmp = 0;
  if (  (fabs(Tintmaxf-Tintmaxo)> criterio ) &&  (fabs(Tintminf-Tintmino)> criterio ) )
    ++tmp;

  return tmp;

}
solve_PQ (double *a,double *b,double *c,double *d,double *P,double *Q,double *Tn,double *T,int nx,
	  double *Tint,double hi,double rhoair,double cair,double La,double *Qin,double dt,
	  double *Tintaverage) {
  int i;
  double Tinn,Qtmp;
  
  Tinn = *Tint;
  P[0] = b[0]/a[0];
  Q[0] = d[0]/a[0];
  for (i = 1 ; i < nx; ++i) {
    P[i] = b[i] / ( a[i] - c[i] * P[i-1] ); 
    Q[i] = ( d[i] + c[i] * Q[i-1] ) / (a[i] - c[i] * P[i-1] );
  }
  Tn[nx-1] = Q[nx-1] ;
  for (i = nx-2; i != -1; --i) 
    Tn[i] = P[i] * Tn[i+1] + Q[i];
  
  for (i = 0; i < nx; ++i) 
  

    for (i = 0; i < nx; ++i) 
    T[i] = Tn[i];  

  //*Tint = (hi*(T[nx-1]-Tinn) + (rhoair*cair*La/dt)*(Tinn))*dt/(rhoair*cair*La);
  *Tint += hi*dt/(rhoair*cair*La)*(T[nx-1]-Tinn)  ;
  *Tintaverage += *Tint;
  if (T[nx-1] > Tinn)
    *Qin += hi*dt*(T[nx-1]-Tinn);

}

double  discomfort_degree_hours(double Tin,double Tc,double dt,double *DDHhot,double *DDHcold) {
  if (Tin < Tc) 
    *DDHcold += (Tc - Tin)*dt/3600.; 
  if (Tin >  Tc)
    *DDHhot  += (Tin - Tc)*dt/3600.; 

}

abrefile(char *user,int month,int cont_herramienta,FILE **f_in) {
  char s_file_capas[140];
  

  sprintf(s_file_capas,"./dat/%s_%d_%d.csv",
	  user,month,cont_herramienta);
  printf(s_file_capas);
  printf("\n");
  // archivo datos temporales NOMBRE
  *f_in = fopen(s_file_capas,"w");
  fprintf(*f_in,"Is\tTsa\tTa\tT_paredext\tTparedint\tTint\tT_n\tDeltaT_n\n");
  fprintf(*f_in,"[W/m2]\t[oC]\t[oC]\t[oC]\t[oC]\t[oC]\t[oC]\t[oC]\n");
}


abrefile_indice(char *user,int cont_herramienta,FILE **f_in,
		double Tsamax,double Tsamin,double Tintmax,double Tintmin,double t_Tintmax,
		double t_Tsamax,double Qin,double Tintaverage,double dt,
		double TPIhot,double TPIcold,double DDHhot,double DDHcold,int month) {
  char file_indices[180];
  double contador;

  contador = 86400/ dt;

  sprintf(file_indices,"./dat/indice_%s_%d_%d.csv",user,month,cont_herramienta);
  // archivo datos temporales NOMBRE
  *f_in = fopen(file_indices,"w");
  printf(file_indices);
  
  fprintf(*f_in,"%.3f\t",Qin/3600.);                  
  fprintf(*f_in,"%.2f\t",(Tintmax-Tintmin)/(Tsamax-Tsamin));                               
  fprintf(*f_in,"%.1f\t",(t_Tintmax-t_Tsamax)/3600.);            
  fprintf(*f_in,"%.1f\t",Tintaverage/contador);                                                     
  fprintf(*f_in,"%.1f\t",Tintmin);                                                                   
  fprintf(*f_in,"%.1f\t",Tintmax);                                                                   
  fprintf(*f_in,"%.1f\t",TPIhot);                                                              
  fprintf(*f_in,"%.1f\t",TPIcold);                                                             
  fprintf(*f_in,"%.1f\t",DDHhot);                                                             
  fprintf(*f_in,"%.1f\t",DDHcold);                                                                 
  fclose(*f_in);
  printf("\n");
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
