/**************** sol.h  *********/
dia_juliano(int *juliano, int dia, int anio, int mes)

{
	int m;
	struct tm f;
	time_t t;
	
	
	memset(&f, 0, sizeof(f));   //llena de ceros las estructuras
	
	//	printf("Cada d'ia = %d  de mes corresponde al calendario juliano a:\n", dia);
	
	for ( m = 1 ; m <= mes ; ++m){
		
		
		f.tm_mday  = dia; 
		f.tm_mon   = m-1;
		f.tm_year = anio - 1900; 
		
		
		t=mktime(&f);
		juliano[m]=f.tm_yday+1;
		
		
	}
	
}

calculo_declinacion_delta(double *delta, int *juliano, int mes)
{
	int m;
	double radianes,pi;
	pi = acos(-1.);
	radianes = pi / 180.;

	for ( m = 1; m <= mes; ++m){
	  delta[m]=23.45*sin(radianes*((360./365.)*(284.+juliano[m])));
	}
}

/****************************************************************************************/
calculo_orto_ocaso(double *orto, double *ocaso, 
				   double phi, double *delta, 
				   int mes)
{
	
  double grados,radianes,pi;
  pi = acos(-1.);
  radianes = pi / 180.;
  grados=180./pi;
  int m;
  
  for ( m=1 ; m<=mes; ++m){
    
    orto[m]=grados*(acos(((-1.)*tan(radianes*phi))*(tan(radianes*delta[m]))));
    ocaso[m]=(-1.)*(orto[m]);
    
  }
}

/****************************************************************************************/
calculo_duracion_dia(double *D_dia, int mes, double *orto)
{
	int m;
	
	for ( m=1 ; m<=mes; ++m){
		
		D_dia[m] = (2. * orto[m])/15.;
		
		//printf("mes= %d\torto = %lf\tduraci'on del d'ia=%lf\n",m,orto[m],D_dia[m]);
	}
	
}

/****************************************************************************************/
calculo_hora_orto(double *t_min, double *D_dia, int mes)
{
	int m;
	
	for ( m=1 ; m<=mes; ++m){
		
		t_min[m] = 12.-(D_dia[m]/2.);
		
		//printf("mes= %d\thora del amanecer=%lf\n",m,t_min[m]);
	}
	
	//printf("--------------------------------------------------\n");
	
}
