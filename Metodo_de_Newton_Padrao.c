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

double eliminacaoGauss(SistLinear_t *SL, double *delta, double**hes, double * grad) {
    triang(SL, delta, m_aux, res);
    retrossubs(SL, delta, m_aux, res);
}

double * calc_grad(SistLinear_t *SL)
{
  double * res = (double*)malloc(SL->num_v*sizeof(double));

  char aux[4];
  char Xn[4];

  for(int i = 0; i < SL->max_iter; i++)
  {
    memset(Xn, 0, sizeof(Xn));
    memset(aux, 0, sizeof(aux));

    for(int l = 0; l < SL->num_v; l++)
    {
      sprintf(aux, "%d", l);
      strcat(strcpy(Xn, "x"), aux);
      res[l] = evaluator_evaluate(SL->GRADIENTE[l], 1, Xn, SL->X[l]);
    }
  }
  return res;
}

double ** calc_hes(SistLinear_t *SL)
{
  double ** m_aux = (double**)malloc(SL->num_v*sizeof(double*));
  for(int i = 0; i < SL->num_v; i++)
  {
    m_aux[i] = (double*)malloc(SL->num_v*sizeof(double));
  }
  char aux[4];
  char Xn[4];

  for (int i = 0; i < SL->num_v; i++)
  {
    for(int l = 0; l < SL->num_v; l++)
    {
      memset(Xn, 0, sizeof(Xn));
      memset(aux, 0, sizeof(aux));
      sprintf(aux, "%d", l);
      strcat(strcpy(Xn, "x"), aux);
      m_aux[i][l] = evaluator_evaluate(SL->HESSIANA[i][l], 1, Xn, SL->X[l]);
    }
  }
  return m_aux;
}

double * Newton_Padrao(SistLinear_t *SL)
{
  for (int i = 0; i < SL->max_iter; i++)
  {
  
    double aux = 0.0;

    double * res = calc_grad(SL);

    double ** m_aux = calc_hes(SL);

    for (int i = 0; i < SL->num_v-1; i++)
    {
      if(aux <= SL->X[i])
      {
        aux = SL->X[i];
      }
    }
    
    if(fabs(aux) < SL->epsilon)
      return SL->X;
    double * delta = (double*) malloc(SL->num_v*sizeof(double));

    delta = eliminacaoGauss(SL, delta, m_aux, res);

    for (int l = 0; l < SL->num_v; l++)
    {
      SL->X[l] = SL->X[l] + delta[l];
    }

    aux = 0.0;

    for (int i = 0; i < SL->num_v-1; i++)
    {
      if(aux <= SL->X[i])
      {
        aux = delta[i];
      }
    }

    if(fabs(aux) < SL->epsilon)
      return SL->X;
  }
}