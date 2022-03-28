#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Metodo_de_Newton_Inexato.h"
#include "SistLinear.h"

// Gauss-Seidel 
double * gaussSeidel(SistLinear_t *SL, double **m_aux, double *grad)
{
  double *res = (double*) calloc(SL->num_v, sizeof(double));
  if (!(res)){
    free(SL);
    printf("ERRO");
    return NULL;
  }
  double *prox_res = (double*) calloc(SL->num_v, sizeof(double));
  if (!(prox_res)){
    free(SL);
    printf("ERRO");
    return NULL;
  }

  int i = 0;
  double num_maior = 0.0;
  double sub = 0.0;

  do{
    prox_res = calcula_independentes(SL, m_aux, grad, res);
    for(int j = 0; j < SL->num_v; j++){
      sub = prox_res[j] - res[j];
      if(sub >= num_maior){
        num_maior = sub;
      }
    }
    res = prox_res;
    i++;
  }while((i <= 50) && (num_maior < SL->epsilon));

  return res;
}

double * calcula_independentes(SistLinear_t *SL, double **m_aux, double *grad, double * res)
{
  double result = 0.0;
  double *indep = res;

  for(int i = 0; i < SL->num_v; i++){
    for(int j = 0; j < SL->num_v; j++)
      result += m_aux[i][j] * indep[j];
    
    result -= m_aux[i][i] * indep[i];
    result = grad[i] - result;
    result = result / m_aux[i][i];
    if (isnan(m_aux[i][i]))
      printf("ERRO: %g\n", m_aux[i][i]);
    indep[i] = result;
    result = 0.0;
  }
  return indep;
}

double ** Newton_Inexato(SistLinear_t *SL, double *TderivadasGS, double * TlsGS)
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
    m_res[0][z] = SL->Xgs[z];

  for (int i = 0; i < SL->max_iter; i++)
  {
    double aux = 0.0;
    double * grad = calc_grad(SL, SL->Xgs, TderivadasGS);
    double ** m_aux = calc_hes(SL, SL->Xgs, TderivadasGS);
  
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
    delta = gaussSeidel(SL, m_aux, grad);
    *TlsGS += timestamp() - tTotal;

    for (int l = 0; l < SL->num_v; l++)
      SL->Xgs[l] += delta[l];

    aux = 0.0;

    for (int i = 0; i < SL->num_v; i++)
      aux += delta[i]*delta[i];
    aux = sqrt(aux);

    for(int z = 0; z < SL->num_v; z++)
      m_res[i+1][z] = SL->Xgs[z];

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