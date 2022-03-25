#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Metodo_de_Newton_Padrao.h"
#include "SistLinear.h"

void pivot(SistLinear_t *SL, int i) {
  double max = fabs(SL->A[i][i]);
  int max_i = i;
  for (int j = i+1; j < SL->num_v; ++j) {
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

void retrossubs(SistLinear_t *SL) {
  for (int i = SL->num_v-1; i >=0; --i) {
    SL->X[i] = SL->b[i];
    for (int j = i+1; j < SL->num_v; j++)
      SL->X[i] -= SL->A[i][j] * SL->X[j];
    SL->X[i] /= SL->A[i][i];
  }
}

void triang(SistLinear_t *SL) {
  for (int i = 0; i < SL->num_v; ++i) {
    pivot(SL, i);
    for (int k = i+1; k < SL->num_v; ++k) {
      double m = SL->A[k][i] / SL->A[i][i];
      if (isnan(m))
        printf("ERRO: %g\n", SL->A[i][i]);
      SL->A[k][i] = 0.0;

      for (int j = i+1; j < SL->num_v; ++j)
        SL->A[k][j] -= SL->A[i][j] * m;
      SL->b[k] -= SL->b[i] * m;
    }
  }
}

  void eliminacaoGauss(SistLinear_t *SL) {
  triang(SL);
  retrossubs(SL);
}