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




/**************** tools.h  *********/
set_krhocv3(double *k,double *rhoc,double dx,
	   double L1,double L2,double L3,
	   double k1,double k2,double *T,
	   double rhoc1,double rhoc2,
	    double *RR,double *Hr,double *Hc,double Aa,double Ac,double a11,double a12,double a21,double Fc,
	    double *Econd,double *Econv, double *Erad,
	    double *T1,double *T2,double *Tave,double *Keq)
{
  int i;
  double Tup,Tdown;
  double hr,hc,E,sigma,Nu; 
  double kair,rhoair,cair;
  double Req1,Req3,Rc;
  double keq;
  double Ra,Nu1,Nu2,Nu3;
  double Beta,nu,alphaair,gr;
  double e22,h,l;
  double wj,wi,tmp1;
  double Fud;
  double Lcond,Lconv;
  
  e22 =h= L2;
  l = 0.164;
  kair    = 0.0262;
  rhoair  = 1.1797660470258469;
  cair    = 1005.458757;  
  sigma = 5.6704e-8;
  E = 0.9;
  gr = 9.81;
  Beta = 1./300.;
  nu = 1.11e-5;  
  alphaair = kair/rhoair/cair;
  Lcond = a11+a12;
  Lconv = a21;
  
  wj = h;  wi = l;
  tmp1 = 0.5* (1. + wj/wi - pow(1.+(wj*wj/wi/wi),0.5));
  Fud = 1. - 2.*(tmp1); 
  //printf("Fud = %f\n",Fud);getchar();

  for (i = 0; i < (L1/dx) ; ++i) {
  }
  Tup = T[i];
  *T1 = Tup;
  *Tave = 0.;
  int cont = 0;
  for (; i <  (L1 + L2)/dx  ; ++i) {
    *Tave += T[i];
    ++cont;
  }
  Tdown = T[i];
  *T2 = Tdown;
  *Tave = *Tave/cont;
  //CALCULO DEL COEFICIENTE CONVECTIVO
  Ra = gr*Beta*fabs(Tup-Tdown)*pow(e22,3.)/nu/alphaair;
  Nu1 = 0.0605*pow(Ra,1./3.);
  Nu2 = pow(  1. + pow(0.104*pow(Ra,0.293)/(1.+pow(6310./Ra,1.36)),3),1./3.);
  Nu3 = 0.242*pow(Ra*h/l,.272);
  if (Nu1>Nu2) Nu = Nu1;
  else Nu = Nu2;
  if (Nu3> Nu) Nu = Nu3;
  hc  = Nu*kair/e22;

  //CALCULO DEL COEFICIENTE RADIATIVO
  //Fc = Aa;
  //Fc = 1.;
  hr = E*sigma*Fud*Fc*(pow(Tup+273.,4.) - pow(Tdown+273.,4.))/(Tup-Tdown);
  if ( (Tup-Tdown) == 0) { hr =0.; }
  *Hr = hr;
  *Hc = hc;

  //MODELO DE BORB'ON 
  double R1,R3,Rr,Rconv,Reqcr,Req2;
  Rconv = 1./(Lconv*hc);
  if (hr == 0.) 
    Reqcr = Rconv;
  else    {    
      Rr = 1./(Lconv*hr);
      Reqcr =  1./(1./Rconv + 1./Rr ); 
  }
  
  R1 = L1/(k1*(Lcond+Lconv));
  R3 = L3/(k1*(Lcond+Lconv));
  Rc = L2/(k1*Lcond);
  //Calculo de las contribuciones en la energ'ia
  if ( (Tup-Tdown) > 0.) {
    *Econd += (Tup-Tdown)/Rc; 
    *Econv += (Tup-Tdown)/Rconv;
    *Erad  += (Tup-Tdown)/Rr;
  }
  Req2 =  1./( 1./Rc + 1./Reqcr);
  keq = L2/(Req2*(Lcond+Lconv));
  *RR = (R1 + Req2 + R3)*(Lcond+Lconv);
  for (i = 0; i < (L1/dx) ; ++i) {   }
  
  for (; i <  (L1 + L2)/dx  ; ++i) {
    k[i] = keq;
    rhoc[i] = rhoc2;
  }
  *Keq = keq;
  //  printf("%f\t%f\t%f\t%f\n",keq,*T1,*T2,*Tave);
  
  
}



/**************** tools.h  *********/
set_krhocv3_techo(double *k,double *rhoc,double dx,
		  double L1,double L2,double L3, double L4,double L5,double L6,
		  double k1,double k2,double k3,double k4,double k5,double k6,
		  double kfinish,double ktop,double kbeam,double khollow,
		  double *T,
		  double rhoc1,double rhoc2,double rhoc3,double rhoc4,double rhoc5,double rhoc6,
		  double *RR,double *Hr,double *Hc,
		  double d1,double d2,double d3,double d4,double d5,double d6,double d7,double d8,double d9,
		  double Aa,double Ac,
		  double Fc,
		  double *Econd,double *Econv, double *Erad,
		  double *T1,double *T2,double *Tave,double *Keq)
{
  int i;
  double Tup,Tdown;
  double hr,hc,E,sigma,Nu; 
  double kair,rhoair,cair;
  double Req1,Req3,Rc;
  double keq;
  double Ra;
  double Beta,nu,alphaair,gr;
  double e22,h,l;
  double wj,wi,tmp1;
  double D;
  double Fud;
  double Req,Rrad;
  
  h = L5;   //Altura de la cavidad
  l = d3;   //Ancho de la cavidad
  D = d1 + d2 +d3+d4+d5+d6+d7+d8+d9;
  kair    = 0.0262;
  rhoair  = 1.1797660470258469;
  cair    = 1005.458757;  
  sigma = 5.6704e-8;
  E = 0.9;
  gr = 9.81;
  Beta = 1./300.;
  nu = 1.11e-5;  
  alphaair = kair/rhoair/cair;
  
  //printf("Fud = %f\n",Fud);getchar();

  for (i = 0; i < ( L1 + L2 + L3 + L4 )/dx ; ++i) {
  }
  Tup = T[i];
  *T1 = Tup;
  *Tave = 0.;
  int cont = 0;
  for (; i <  (L1 + L2 + L3 + L4 + L5)/dx  ; ++i) {
    *Tave += T[i];
    ++cont;
  }
  Tdown = T[i];
  *T2 = Tdown;
  *Tave = *Tave/cont;
  //CALCULO DEL COEFICIENTE CONVECTIVO
  double dot11,dot22;
  Ra = gr*Beta*(Tup-Tdown)*pow(L5,3.)/nu/alphaair;
  if (Ra  <1e-5) {
    Ra = 0.;
    dot11 = 0.;
  }
  else dot11 = 1. - 1708./Ra;
  dot22 = pow(Ra/5830.,1./3.)-1.;
  if ( dot11<0.) dot11 = 0.;
  if  (dot22<0.) dot22 = 0.;
  hc  = kair/L5*(1. + 1.44*dot11 + dot22);
  if ( (Tdown  <Tup ) || (Ra ==0.) )
    hc = kair/L5;
  double Rconv,Rcond;
  Rconv  = 1./(hc*(d3 +d5 +d7 )/D);
  
  //Radiative coefficient for the roof
  wj = h;  wi = l;
  tmp1 = 0.5* (1. + wj/wi - pow(1.+(wj*wj/wi/wi),0.5));
  Fud = 1. - 2.*(tmp1); 
  Rcond = pow( 1./ (L5/(kbeam * d1)/(2.*D)) + 1./(L5/(khollow*(d1/2.+d2 +d4 +d6 +d8 +d9/2.)/D)) + 1./ (L5/(kbeam * d9)/(2.*D)) , -1.);
  //CALCULO DEL COEFICIENTE RADIATIVO
  if  (  fabs(Tup-Tdown) <1e-5 ) hr = 0.;
  else  hr = E*sigma*Fud*Fc*(pow(Tup+273.,4.) - pow(Tdown+273.,4.))/(Tup-Tdown);
  if (hr  == 0.) {
    hr = 0.;
    // printf("Rcon_r  if = %f\n",Rcon_r);
    // printf("Rconv_r if = %f\n",Rconv_r);getchar();
    Req = pow( 1./Rcond + 1./Rconv, -1.);
    //printf("Req_r i= %f\n",Req_r);getchar();
  }
  else {
    //printf("Tup_r - Tdown_r = %f\n",Tup_r-Tdown_r);
    Rrad = 1./(hr*(d3+d5+d7)/D);
    Req = pow( 1./Rcond + 1./Rconv + 1./Rrad , -1.);
  }
  
  keq = L5/Req;
  
  
  for (i = 0; i < (L1+L2+L3)/dx ; ++i) {   }
  
  for (; i <  (L1 + L2 + L3 + L4 )/dx  ; ++i) {
    k[i] = keq;
    rhoc[i] = rhoc2;
  }
  *Keq = keq;
  //  printf("%f\t%f\t%f\t%f\n",keq,*T1,*T2,*Tave);
  
  
}



/**************** tools.h  *********/
set_krhoc2(double *k,double *rhoc,double dx,
	   double L1,double L2,double L3,
	   double k1,double k2,double *T,
	   double rhoc1,double rhoc2,
	   double *RR,double *Hr,double *Hc,double Aa,double Ac)
{
  int i;
  double Tup,Tdown;
  double hr,hc,E,sigma,Nu; 
  double kair,rhoair,cair;
  double Req1,Req2,Req3,Rc;
  double keq;
  double Ra,Nu1,Nu2,Nu3;
  double Beta,nu,alphaair,gr;
  double e22,h,l;
  double Fud,wi,wj,tmp1;
   
  
  e22 =h= 0.07;
  l = 0.164;
  kair    = 0.0262;
  rhoair  = 1.1797660470258469;
  cair    = 1005.458757;  
  sigma = 5.6704e-8;
  E = 0.9;
  gr = 9.81;
  Beta = 1./300.;
  nu = 1.11e-5;  
  alphaair = kair/rhoair/cair;

  wj = h;  wi = l;
  tmp1 = 0.5* (1. + wj/wi - pow(1.+(wj*wj/wi/wi),0.5));
  Fud = 1. - 2.*(tmp1); 



  for (i = 0; i < (L1/dx) ; ++i) {
  }
  Tup = T[i];
  for (; i <  (L1 + L2)/dx  ; ++i) {
  }
  Tdown = T[i];
  

  Ra = gr*Beta*fabs(Tup-Tdown)*pow(e22,3.)/nu/alphaair;
  Nu1 = 0.0605*pow(Ra,1./3.);
  Nu2 = pow(  1. + pow(0.104*pow(Ra,0.293)/(1.+pow(6310./Ra,1.36))  ,3) ,1./3.);
  Nu3 = 0.242*pow(Ra*h/l,.272);
  if (Nu1>Nu2) Nu = Nu1;
  else Nu = Nu2;
  if (Nu3> Nu) Nu = Nu3;
  hc  = Nu*kair/e22;
  Nu = hc*L2/kair;
  if (Nu<1.) hc = kair/L2;


  double Fc;
  //printf("Fud = %f\n",Fud);getchar();
  Fc = Aa;
  hr = E*sigma*Fud*Fc*(pow(Tup+273.15,4.) - pow(Tdown+273.15,4.))/(Tup-Tdown);
  if ( (Tup-Tdown) == 0) { hr =0.; }
  *Hr = hr;
  *Hc = hc;
  double R1,R3,Rr,Rconv;
  //MODELO DE BORB'ON
  /* printf("Aa = %f\n",Aa); */
  /* printf("Ac = %f\n",Ac);//getchar(); */
  Rr = 1./(Aa*hr);
  Rconv = 1./(Aa*hc);
  if (hr == 0.)   Req1 = Rconv;
  else Req1 =  1./(1./Rconv + 1./Rr );
  R1 = L1/(k1*Aa);
  R3 = L3/(k1*Aa);
  Rc = (L1+L2+L3)/(k1*Ac);
  Req2 = R1 + Req1 + R3 ;
  Req3 = 1./ (1./Rc  + 1./Req2);
  // printf("Req3 = %f\n",Req3);getchar();




   keq = L2/(Req3 - (L1+L3)/k1);
   *RR = Req3;

  for (i = 0; i < (L1/dx) ; ++i) {
  }
  
  for (; i <  (L1 + L2)/dx  ; ++i) {
    k[i] = keq;
    rhoc[i] = rhoc2;
    //printf("rhoc2 = %f\n",rhoc2);getchar();
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

max_min(double *Tsemax,double *Tsemin,double *Tsintmax,double *Tsintmin,
	double t,double Tint,double *T,int nx) {

  
  if (T[0] > *Tsemax)  
    *Tsemax  = T[0];
  if (T[0]  < *Tsemin) 
    *Tsemin  = T[0];
  
  
  if (T[nx-1] < *Tsintmin) 
    *Tsintmin = T[nx-1];
  
  if (T[nx-1] > *Tsintmax) 
    *Tsintmax = T[nx-1];


}


double time_evolution_Tsa(double *Ta,double t,double Ig,
			  double a,double ho,double Tmax,double Tmin,double Hi,
			  double Ho,double tau,double delta,double Ib,double Id,double gamma,
			  double phi,double beta) {
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
	  double Tint,double hi,double rhoair,double cair,double La,double *Qin,double dt,
	  double *Tintaverage,double *Qcool,double *Qheat) {
  int i;
  double Tinn,Qtmp;
  
  Tinn = Tint;
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
  
  
  if (T[nx-1] > Tint)
    *Qcool += hi*dt*(T[nx-1]-Tint);
  
  if (T[nx-1] < Tint)
    *Qheat += hi*dt*(Tint - T[nx-1]);
}

double  discomfort_degree_hours(double Tin,double Tc,double dt,double *DDHhot,double *DDHcold) {
  if (Tin < Tc) 
    *DDHcold += (Tc - Tin)*dt/3600.; 
  if (Tin >  Tc)
    *DDHhot  += (Tin - Tc)*dt/3600.; 

}

void abrefile_techo(double a,double gamma,double beta,FILE **f_in,double Lat,int month) {
  char s_file_capas[140];
  

	
  sprintf(s_file_capas,"./dat/techo_Lat%.2f_a%.1f_mes%d.csv",
	  Lat,a,month);
 
 
  printf(s_file_capas);
  printf("\n");


  // archivo datos temporales NOMBRE
  *f_in = fopen(s_file_capas,"w");
  fprintf(*f_in,"#Hora\tTint\tTa\tTsa\tHr\tHc\tT1\tT2\t<T12>\tKeq2\n");
  fprintf(*f_in,"#1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n");

}

void abrefile_muro(double a,double gamma,double beta,FILE **f_in,double Lat,int month) {
  char s_file_capas[140];
  
  char orientation[140];

  if (gamma == 0. ) sprintf(orientation,"sur");
  if (gamma == 180. ) sprintf(orientation,"norte");
  if (gamma == -90. ) sprintf(orientation,"este");
  if (gamma == 90. ) sprintf(orientation,"oeste");
  
  
  
  
  sprintf(s_file_capas,"./dat/muro_Lat%.2f_a%.1f_mes%d_%s.csv",
	  Lat,a,month,orientation);
	
 
 
  printf(s_file_capas);
  printf("\n");


  // archivo datos temporales NOMBRE
  *f_in = fopen(s_file_capas,"w");
  fprintf(*f_in,"#Hora\tTint\tTa\tTsa\tHr\tHc\tT1\tT2\t<T12>\tKeq2\n");
  fprintf(*f_in,"#1\t2\t3\t4\t5\t6\t7\t8\t9\t10\n");

}

abrefile_indice(double L1,double L2,double L3,double k1,double k2,double k3,double rhoc1,
		double rhoc2,double rhoc3,double La,char *user,int month,
		int cont_herramienta,FILE **f_in,
		double Tsamax,double Tsamin,double Tintmax,double Tintmin,double t_Tintmax,
		double t_Tsamax,double Qin,double Tintaverage,double dt,
		double TPIhot,double TPIcold,double DDHhot,double DDHcold) {
  char file_indices[180];
  double contador;

  contador = 86400/ dt;

  sprintf(file_indices,"./dat/%s_indices_%d_%d_Lext%.2f_Lmed%.2f_Lint%.2f_Kext%.3f_Kmed%.3f_Kint%.3f_rhoCext%.3e_rhoCmed%.3e_rhoCint%.3e_La%.1f.dat",user,cont_herramienta,month,L1,L2,L3,k1,k2,k3,rhoc1,rhoc2,rhoc3,La);
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
