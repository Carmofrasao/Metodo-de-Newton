#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"
#include "Metodo_de_Newton_Modificado.h"
#include "Metodo_de_Newton_Padrao.h"

void pivotLU(SistLinear_t *SL, int i, double**hes, double * grad) {
  double max = fabs(hes[i][i]);
  int max_i = i;
  for (int j = i+1; j < SL->num_v; ++j) {
    double v = fabs(hes[j][i]);
    if (v > max) {
      max = v;
      max_i = j;
    }
  }
  if (max_i != i) {
    double *tmp = hes[i];
    hes[i] = hes[max_i];
    hes[max_i] = tmp;

    double aux = grad[i];
    grad[i] = grad[max_i];
    grad[max_i] = aux;
  }
} 

void pivotz(SistLinear_t *SL, int i) {
  double max = fabs(SL->L[i][i]);
  int max_i = i;
  for (int j = i+1; j < SL->num_v; ++j) {
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
  for (int i = 0; i < SL->num_v; ++i) {
    SL->z[i] = SL->b[i];
    for (int j = 0; j < i; j++)
      SL->z[i] -= SL->L[i][j] * SL->z[j];
  }
}

void retrossubs2(SistLinear_t *SL, double**hes, double *delta) {
  for (int i = SL->num_v-1; i >=0; --i) {
    delta[i] = SL->z[i];
    for (int j = i+1; j < SL->num_v; j++)
      delta[i] -= hes[i][j] * delta[j];
    delta[i] /= hes[i][i];
  }
}

void triangLU(SistLinear_t *SL, double**hes, double * grad) {
  for (int i = 0; i < SL->num_v; ++i) {
    pivotLU(SL, i, hes, grad);
    pivotz(SL, i);
    if (i == SL->num_v-1)
        SL->L[i][i] = 1;
    for (int k = i+1; k < SL->num_v; ++k) {
      double m = hes[k][i] / hes[i][i];
      if (isnan(m))
        printf("ERRO: %g\n", hes[i][i]);
      SL->U[k][i] = 0.0;
      SL->L[k][i] *= m;

      if (k-1 == i)
        SL->L[k-1][i] = 1;

      for (int j = i+1; j < SL->num_v; ++j)
        hes[k][j] -= hes[i][j] * m;
      grad[k] -= grad[i] * m;
    }
  }
}

double * FatLU(SistLinear_t *SL, double *delta, double**hes, double * grad) {
  triangLU(SL, hes, grad);
  retrossubpz(SL);
  retrossubs2(SL, hes, delta);
  return delta;
}

double ** Newton_Modificado(SistLinear_t *SL)
{
  double ** m_res = (double**) calloc(SL->max_iter+1, sizeof(double*));
  for(int i = 0; i < SL->max_iter+1; i++)
  { 
    m_res[i] = (double*) calloc(SL->num_v, sizeof(double));
  }

  for (int i = 0; i < SL->max_iter; i++)
  {
    double aux = 0.0;

    double * grad = calc_grad(SL, SL->Xlu);
    
    double ** m_aux = calc_hes(SL, SL->Xlu);
  
    for (int i = 0; i < SL->num_v; i++)
    {
      if(aux <= fabs(grad[i]))
      {
        aux = grad[i];
      }
    }

    m_res[i] = SL->Xlu;

    if(fabs(aux) < SL->epsilon)
      return m_res;
    double * delta = (double*) calloc(SL->num_v, sizeof(double));

    delta = FatLU(SL, delta, m_aux, grad);

    for (int l = 0; l < SL->num_v; l++)
    {
      SL->Xeg[l] += delta[l];
    }

    aux = 0.0;

    for (int i = 0; i < SL->num_v-1; i++)
    {
      if(aux <= delta[i])
      {
        aux = delta[i];
      }
    }

    m_res[i+1] = SL->Xeg;

    if(fabs(aux) < SL->epsilon)
      return m_res;
  }
  return NULL;
}