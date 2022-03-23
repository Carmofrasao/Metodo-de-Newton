#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"
#include "Metodo_de_Newton_Modificado.h"

void pivotLU(SistLinear_t *SL, int i) {
  double max = fabs(SL->A[i][i]);
  int max_i = i;
  for (int j = i+1; j < SL->n; ++j) {
    double v = fabs(SL->A[j][i]);
    if (v > max) {
      max = v;
      max_i = j;
    }
  }
  if (max_i != i) {
    double *tmp = SL->A[i];
    SL->A[i] = SL->A[max_i];
    SL->A[max_i] = tmp;

    double aux = SL->b[i];
    SL->b[i] = SL->b[max_i];
    SL->b[max_i] = aux;
  }
} 

void pivotz(SistLinear_t *SL, int i) {
  double max = fabs(SL->L[i][i]);
  int max_i = i;
  for (int j = i+1; j < SL->n; ++j) {
    double v = fabs(SL->L[j][i]);
    if (v > max) {
      max = v;
      max_i = j;
    }
  }
  if (max_i != i) {
    double *tmp = SL->L[i];
    SL->L[i] = SL->L[max_i];
    SL->L[max_i] = tmp;

    double aux = SL->z[i];
    SL->z[i] = SL->z[max_i];
    SL->z[max_i] = aux;
  }
}

void retrossubpz(SistLinear_t *SL) {
  for (int i = 0; i < SL->n; ++i) {
    SL->z[i] = SL->b[i];
    for (int j = 0; j < i; j++)
      SL->z[i] -= SL->L[i][j] * SL->z[j];
  }
}

void retrossubs2(SistLinear_t *SL, double *X) {
  for (int i = SL->n-1; i >=0; --i) {
    X[i] = SL->z[i];
    for (int j = i+1; j < SL->n; j++)
      X[i] -= SL->A[i][j] * X[j];
    X[i] /= SL->A[i][i];
  }
}

void triangLU(SistLinear_t *SL) {
  for (int i = 0; i < SL->n; ++i) {
    pivotLU(SL, i);
    pivotz(SL, i);
    if (i == SL->n-1)
        SL->L[i][i] = 1;
    for (int k = i+1; k < SL->n; ++k) {
      double m = SL->A[k][i] / SL->A[i][i];
      if (isnan(m))
        printf("ERRO: %g\n", SL->A[i][i]);
      SL->U[k][i] = 0.0;
      SL->L[k][i] *= m;

      if (k-1 == i)
        SL->L[k-1][i] = 1;

      for (int j = i+1; j < SL->n; ++j)
        SL->A[k][j] -= SL->A[i][j] * m;
      SL->b[k] -= SL->b[i] * m;
    }
  }
}

void FatLU(SistLinear_t *SL, double *X) {
    triangLU(SL);
    retrossubpz(SL);
    retrossubs2(SL, X);
}
