#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <matheval.h>

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

double * calc_grad(SistLinear_t *SL)
{
  double * res = (double*)malloc(SL->num_v*sizeof(double));

  char aux[4];
  char Xn[4];
  char **X = (char**) malloc(SL->num_v*sizeof(char*));

  for(int i = 0; i < SL->max_iter; i++)
  {
    memset(Xn, 0, sizeof(Xn));
    memset(aux, 0, sizeof(aux));

    for(int l = 0; l < SL->num_v; l++)
    {
      sprintf(aux, "%d", l+1);
      strcat(strcpy(Xn, "x"), aux);
      X[l] = Xn;
    }
    for(int l = 0; l < SL->num_v; l++)
    {
      res[l] = evaluator_evaluate(SL->GRADIENTE[l], SL->num_v, X, SL->Xeg);
    }
  }
  return res;
}

double ** calc_hes(SistLinear_t *SL)
{
  double ** m_aux = (double**)calloc(SL->num_v, sizeof(double*));
  for(int i = 0; i < SL->num_v; i++)
  {
    m_aux[i] = (double*)calloc(SL->num_v, sizeof(double));
  }
  char aux[4];
  char Xn[4];
  char **X = (char**) calloc(SL->num_v, sizeof(char*));

  for(int l = 0; l < SL->num_v; l++)
  {
    memset(Xn, 0, sizeof(Xn));
    memset(aux, 0, sizeof(aux));
    sprintf(aux, "%d", l+1);
    strcat(strcpy(Xn, "x"), aux);
    X[l] = Xn;
  }

  for (int i = 0; i < SL->num_v; i++)
  {
    for(int l = 0; l < SL->num_v; l++)
    {
      m_aux[i][l] = evaluator_evaluate(SL->HESSIANA[i][l], SL->num_v, X, SL->Xeg);
    }
  }

  return m_aux;
}

double ** Newton_Padrao(SistLinear_t *SL)
{
  double ** m_res = (double**) calloc(SL->max_iter, sizeof(double*));
  for(int i = 0; i < SL->max_iter; i++)
  { 
    m_res[i] = (double*) calloc(SL->max_iter, sizeof(double));
  }

  for (int i = 0; i < SL->max_iter; i++)
  {
    double aux = 0.0;

    double * grad = calc_grad(SL);
    
    double ** m_aux = calc_hes(SL);

    for (int i = 0; i < SL->num_v; i++)
    {
      if(aux <= fabs(grad[i]))
      {
        aux = grad[i];
      }
    }

    m_res[i] = SL->Xeg;

    if(fabs(aux) < SL->epsilon)
      return m_res;
    double * delta = (double*) calloc(SL->num_v, sizeof(double));

    delta = eliminacaoGauss(SL, delta, m_aux, grad);

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