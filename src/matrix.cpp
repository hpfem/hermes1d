#include "matrix.h"

#define TINY 1e-20


void ludcmp(double** a, int n, int* indx, double* d)
{
  int i, imax = 0, j, k;
  double big, dum, sum, temp;
  double* vv = new double[n];

  *d = 1.0;
  for (i = 0; i < n; i++)
  {
    big=0.0;
    for (j = 0; j < n; j++) 
      if ((temp = fabs(a[i][j])) > big) 
        big = temp;
    if (big == 0.0) error("Singular matrix!");
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a[i][j];
      for (k = 0; k < i; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)
    {
      sum = a[i][j];
      for (k = 0; k < j; k++) sum -= a[i][k]*a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i]*fabs(sum)) >= big)
      {
        big = dum;
        imax = i;
      }
    }
    if (j != imax)
    {
      for (k = 0; k < n; k++)
      {
        dum = a[imax][k];
        a[imax][k] = a[j][k];
        a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0) a[j][j] = TINY;
    if (j != n-1) 
    {
      dum = 1.0 / (a[j][j]);
      for (i = j+1; i < n; i++) a[i][j] *= dum;
    }
  }
  delete [] vv;
}


void choldc(double **a, int n, double p[])
{
  int i, j, k;
  for (i = 0; i < n; i++)
  {
    for (j = i; j < n; j++)
    {
      double sum = a[i][j];
      k = i;
      while (--k >= 0)
        sum -= a[i][k] * a[j][k];
      
      if (i == j)
      {
        if (sum <= 0.0)
          error("CHOLDC failed!");
        else
          p[i] = sqrt(sum);
      }
      else
        a[j][i] = sum / p[i];
    }
  }
}

void solve_linear_system_dense(DenseMatrix *mat, double *res)
{
    int n = mat->get_size();
    int *indx = new int[n];
    double **_mat = mat->get_mat();
    double d;
    ludcmp(_mat, n, indx, &d);
    lubksb(_mat, n, indx, res);
}
