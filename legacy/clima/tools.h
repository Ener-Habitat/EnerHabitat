

time_t doy2tm(const int year, const int dayofyear, struct tm *tm) {

  int leapyear, m, d;

  int d_in_m [] = {31,28,31,30,31,30,31,31,30,31,30,31};

  memset(tm, 0, sizeof(*tm));

  leapyear = (year % 4 == 0) && ((year % 100 != 0) || (year % 400 == 0));

  if( dayofyear < 1 || dayofyear > (leapyear ? 366 : 365) ) {
    return (time_t) -1;
  }

  if(leapyear) { d_in_m[1]=29; }

  for(m=0, d=dayofyear; 
      d>d_in_m[m]; 
      d-=d_in_m[m++]);

  tm->tm_mday=d;
  tm->tm_mon=m;
  tm->tm_year=year - 1900;

  return mktime(tm);
}


df_lt(double *DF,double *LT,double **To,double **Ti,int d1,int d2,int ndats,double dtat,
      double *Toutave,double *Tinave)
{
  int i,j;

  

  //calcula el valor de temperatura m'inima de cada d'ia y el momento en el que sucede
  double *Tmino,*Tmini;
  double *tTmino,*tTmini;
  double mini,mino;
  Tmino = one_d_double_array(365);
  tTmino = one_d_double_array(365);
  Tmini = one_d_double_array(365);
  tTmini = one_d_double_array(365);

  for (j = d1; j < d2+1; ++j) {
    mini = mino = 1000.;
    for (i = 0; i < ndats; ++i) {
      Tinave[j] += Ti[i][j];
      Toutave[j] += To[i][j];
      //printf("en el If To[%d][%d] = %f\n",i,j,To[i][j]);getchar();
      if (To[i][j]<=mino) {
	Tmino[j] = To[i][j];
	mino = Tmino[j];
	//printf("Tmino[%d] = %f\n",j,Tmino[j]);getchar();
	tTmino[j] = (double)i;
      }
      if (Ti[i][j]<=mini) {
	Tmini[j] = Ti[i][j];
	mini = Tmini[j];
	tTmini[j] =  (double) i;
      }
    }
  }
  //calcula el valor de temperatura m'axima de cada d'ia y el momento en el que sucede
  double *Tmaxo,*Tmaxi;
  double *tTmaxo,*tTmaxi;
  double maxi,maxo;
  
  Tmaxo = one_d_double_array(365);
  tTmaxo = one_d_double_array(365);
  Tmaxi = one_d_double_array(365);
  tTmaxi = one_d_double_array(365);
  
  for (j = d1; j < d2+1; ++j) {
    maxi = maxo = -1000.;
    for (i = 0; i < ndats; ++i) {
      if (To[i][j]>=maxo) {
	Tmaxo[j] = To[i][j];
	maxo = Tmaxo[j];
	//tTmaxo[j] = (double) i;
      }
      if (Ti[i][j]>=maxi) {
	Tmaxi[j] = Ti[i][j];
	maxi = Tmaxi[j];
	//tTmaxi[j] = (double) i;
      }
    }
  }
  //Calculo del tiempo de ocurrencia las temperaturas m'aximas filtrando el medio dia solar
  for (j = d1; j < d2+1; ++j) {
    maxi = maxo = -1000.;
    for (i = ndats/2; i < ndats; ++i) {
      if (To[i][j]>=maxo) {
	maxo = To[i][j];
	tTmaxo[j] = (double) i;
      }
      if (Ti[i][j]>=maxi) {
	maxi = Ti[i][j];
	tTmaxi[j] = (double) i;
      }
    }
  }
  //printf("#d1 = %d\td2 = %d\n",d1,d2);
  //printf("#day\tDF\tLT\t<Tin>\t<Tout>\n");

  for (i = d1; i < d2+1; ++i) {
    //printf("day = %d, Tmaxi = %f\t Tmini = %f\n",i,Tmaxi[i],Tmini[i]); 
    //printf("day = %d, Tmaxo = %f\t Tmino = %f\n",i,Tmaxo[i],Tmino[i]); getchar();
    Tinave[i] = Tinave[i]/ndats;
    DF[i] = ( Tmaxi[i] - Tmini[i]) / (Tmaxo[i] - Tmino[i] );
    LT[i] = (tTmaxi[i] - tTmaxo[i])*dtat/60.;
    //printf("%d\t%.2f\t%.2f\t%.2f\t%.2f\n",i,DF[i],LT[i],Tinave[i],Toutave[i]/ndats);
  }
  
}



double varianza  (double *DF, double *LT,double *Tinave,int d1,int d2) 
{
  
  int i,contador;
  DF_ave = LT_ave = Tave = 0.;
  contador = 0;
  for (i = d1; i < d2+1; ++i) {
    DF_ave += DF[i];
    LT_ave += LT[i];
    Tave += Tinave[i];
    ++contador;
  }
  
  DF_ave = DF_ave / contador;
  LT_ave = LT_ave / contador;
  Tave   = Tave/contador;

  DF_sigma = LT_sigma = Tin_sigma = 0.;
  
  for (i = d1; i < d2 + 1 ; ++i) {
    DF_sigma  += pow(DF[i] - DF_ave,2.);
    LT_sigma  += pow(LT[i] - LT_ave,2.);
    Tin_sigma += pow(Tinave[i] - Tave,2.);
  }
  
  DF_sigma  = pow(DF_sigma/contador,0.5);
  LT_sigma  = pow(LT_sigma/contador,0.5);
  Tin_sigma = pow(Tin_sigma/contador,0.5);
  
  //printf("DF_sigma  = %f\n",DF_sigma);
  //printf("LT_sigma  = %f\n",LT_sigma);

}
