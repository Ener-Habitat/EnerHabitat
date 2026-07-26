/**************** arrays.h  *********/


int *one_d_int_array(int n)
{
  int *ptr;

  ptr = (int *)malloc(n*sizeof(int));
  return ptr;
}

double *one_d_double_array(int n)
{
  double *ptr;

  ptr = (double *)malloc(n*sizeof(double));
  return ptr;
}

/*	s[i][j], 0<=i<lx, 0<=j<ly	*/
int **two_d_int_array(int lx, int ly) 
{
  int **ptr, i;
  
  ptr = (int **) malloc(lx*sizeof(int*));
  for(i=0;i<lx;++i) 
    ptr[i] = (int *) malloc(ly*sizeof(int));
  return ptr;
}

/*	s[i][j], 0<=i<lx, 0<=j<ly	*/
double **two_d_double_array(int lx, int ly) 
{
  int i;
  double  **ptr;

  ptr = (double **) malloc(lx*sizeof(double*));
  for(i=0;i<lx;++i) 
    ptr[i] = (double*) malloc(ly*sizeof(double));
  return ptr;
}

free_two_d_int_array(int n,int **s)
{
   int i;
   
   for (i=0;i<n;++i)
      free(s[i]);
}

free_two_d_double_array(int n,double **s)
{
   int i;
   
   for (i=0;i<n;++i)
      free(s[i]);
}

/*	s[i][j][k], 0<=i<lx, 0<=j<ly, 0<=k<lz  */
int ***three_d_int_array(int lx, int ly, int lz) 
{
  int ***ptr, i,j;
  
  ptr = (int ***) malloc(lx*sizeof(int**));
  for(i=0;i<lx;++i) {
     ptr[i] = (int **) malloc(ly*sizeof(int*));
    for(j=0;j<ly;++j) {
      ptr[i][j] = (int *) malloc(lz*sizeof(int));
    }
  }
  return ptr;
}

double ***three_d_double_array(int lx, int ly, int lz) 
{
  int i,j;
  double ***ptr;
  
  ptr = (double ***) malloc(lx*sizeof(double**));
  for(i=0;i<lx;++i) {
    ptr[i] = (double **) malloc(ly*sizeof(double*));
    for(j=0;j<ly;++j) {
      ptr[i][j] = (double *) malloc(lz*sizeof(double));
    }
  }
  return ptr;
}

free_three_d_int_array(int lx, int ly, int ***s)
{
  int i,j;
   
  for (i=0;i<lx;++i) {
    for (j=0;j<ly;++j) {
      free(s[i][j]);
    }
    free(s[i]);
  }
}

free_three_d_double_array(int lx, int ly, double ***s)
{
  int i,j;
   
  for (i=0;i<lx;++i) {
    for (j=0;j<ly;++j) {
      free(s[i][j]);
    }
    free(s[i]);
  }
}


/*	s[i][j][k][l], 0<=i<lx, 0<=j<ly, 0<=k<lz, 0<=l<t */
int ****four_d_int_array(int lx, int ly, int lz, int lt) 
{
  int ****ptr, i,j,k;
  
  ptr = (int ****) malloc(lx*sizeof(int***));
  for(i=0;i<lx;++i) {
    ptr[i] = (int ***) malloc(ly*sizeof(int**));
    for(j=0;j<ly;++j) {
      ptr[i][j] = (int **) malloc(lz*sizeof(int*));
      for(k=0;k<lz;++k) {
	ptr[i][j][k] = (int *) malloc(lt*sizeof(int));
      }
    }
  }
  return ptr;
}

free_four_d_int_array(int lx, int ly, int lz, int ****s)
{
  int i,j,k;
   
  for (i=0;i<lx;++i) {
    for (j=0;j<ly;++j) {
      for (k=0;k<lz;++k)
	 free(s[i][j][k]);
      free(s[i][j]);
    }
    free(s[i]);
  }
}

/*	s[i][j][k][l], 0<=i<lx, 0<=j<ly, 0<=k<lz, 0<=l<t */
double ****four_d_double_array(int lx, int ly, int lz, int lt) 
{
  double ****ptr;
  int i,j,k;
  
  ptr = (double ****) malloc(lx*sizeof(double***));
  for(i=0;i<lx;++i) {
    ptr[i] = (double ***) malloc(ly*sizeof(double**));
    for(j=0;j<ly;++j) {
      ptr[i][j] = (double **) malloc(lz*sizeof(double*));
      for(k=0;k<lz;++k) {
	ptr[i][j][k] = (double *) malloc(lt*sizeof(double));
      }
    }
  }
  return ptr;
}

free_four_d_double_array(int lx, int ly, int lz, double ****s)
{
  int i,j,k;
   
  for (i=0;i<lx;++i) {
    for (j=0;j<ly;++j) {
      for (k=0;k<lz;++k)
	 free(s[i][j][k]);
      free(s[i][j]);
    }
    free(s[i]);
  }
}

/////////////////////////////////////////////////////////////////////////////////
