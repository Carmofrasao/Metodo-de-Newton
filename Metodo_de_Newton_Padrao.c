#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "Metodo_de_Newton_Padrao.h"
#include "SistLinear.h"

void pivot(SistLinear_t *SL, double**hes, double * grad, int i) {
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

void retrossubs(SistLinear_t *SL, double *delta, double**hes, double * grad) {
  for (int i = SL->num_v-1; i >=0; --i) {
    delta[i] = grad[i];
    for (int j = i+1; j < SL->num_v; j++)
      delta[i] -= hes[i][j] * delta[j];
    delta[i] /= hes[i][i];
  }
}

void triang(SistLinear_t *SL, double**hes, double * grad) {
  for (int i = 0; i < SL->num_v; ++i) {

    pivot(SL, hes, grad, i);

    for (int k = i+1; k < SL->num_v; ++k) {
      double m = hes[k][i] / hes[i][i];
      if (isnan(m))
        printf("ERRO: %g\n", hes[i][i]);
      hes[k][i] = 0.0;

      for (int j = i+1; j < SL->num_v; ++j)
        hes[k][j] -= hes[i][j] * m;
      grad[k] -= grad[i] * m;
    }
  }
}

double* eliminacaoGauss(SistLinear_t *SL, double *delta, double**hes, double * grad) {
    triang(SL, hes, grad);
    retrossubs(SL, delta, hes, grad);
    return delta;
}

double ** Newton_Padrao(SistLinear_t *SL, double *TderivadasEG, double *TslEG)
{
  double ** m_res = (double**) calloc(SL->max_iter+1, sizeof(double*));
  if (!(m_res)){
    free(SL);
    printf("ERRO");
    return NULL;
  }
  for(int i = 0; i < SL->max_iter+1; i++)
  { 
    m_res[i] = (double*) calloc(SL->num_v, sizeof(double));
    if (!(m_res[i])){
      free(SL);
      printf("ERRO");
      return NULL;
    }
  }

  for(int z = 0; z < SL->num_v; z++)
    m_res[0][z] = SL->Xeg[z];

  for (int i = 0; i < SL->max_iter; i++)
  {
    double aux = 0.0;
    double * grad = calc_grad(SL, SL->Xeg, TderivadasEG);
    double ** m_aux = calc_hes(SL, SL->Xeg, TderivadasEG);
  
    for (int i = 0; i < SL->num_v; i++)
      aux += grad[i]*grad[i];
    aux = sqrt(aux);

    if(aux < SL->epsilon)
    {
      for (int l = i+1; l < SL->max_iter+1; l++)
        for(int z = 0; z < SL->num_v; z++)
          m_res[l][z] = NAN;
      return m_res;
    }
    
    double * delta = (double*) calloc(SL->num_v, sizeof(double));
    if (!(delta)){
      free(SL);
      printf("ERRO");
      return NULL;
    }

    double tTotal = timestamp();
    delta = eliminacaoGauss(SL, delta, m_aux, grad);
    *TslEG += timestamp() - tTotal;

    for (int l = 0; l < SL->num_v; l++)
      SL->Xeg[l] += delta[l];

    aux = 0.0;

    for (int i = 0; i < SL->num_v; i++)
      aux += delta[i]*delta[i];
    aux = sqrt(aux);

    for(int z = 0; z < SL->num_v; z++)
      m_res[i+1][z] = SL->Xeg[z];

    if(aux < SL->epsilon)
    {
      for (int l = i+1; l < SL->max_iter; l++)
        for(int z = 0; z < SL->num_v; z++)
          m_res[l+1][z] = NAN;
      return m_res;
    }
  }
  return m_res;
}