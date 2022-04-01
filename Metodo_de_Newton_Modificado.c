/* Código escrito por
    Anderson Aparecido do Carmo Frasão 
    GRR 20204069
    Erick Eckermann Cardoso 
    GRR 20186075
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"
#include "Metodo_de_Newton_Modificado.h"

//função para aplicar o pivoteamento na matriz hessiana
void pivotLU(SistLinear_t *SL, int i, double**hes, double * grad) {
  double max = fabs(hes[i][i]);
  int max_i = i;
   //verifica qual posição tem o maior valor na coluna i, linha j (i+1)
  for (int j = i+1; j < SL->num_v; ++j) {
    double v = fabs(hes[j][i]);
    if (v > max) {
      max = v;
      max_i = j;
    }
  }
  //pivoteamento
  if (max_i != i) {
    double *tmp = hes[i];
    hes[i] = hes[max_i];
    hes[max_i] = tmp;

    double aux = grad[i];
    grad[i] = grad[max_i];
    grad[max_i] = aux;
  }
} 

//função para aplicar o pivoteamento na matriz LESS com o vetor z (fatoração LU)
void pivotz(SistLinear_t *SL, int i) {
  double max = fabs(SL->L[i][i]);
  int max_i = i;
  //verifica qual posição tem o maior valor na coluna i, linha j (i+1)
  for (int j = i+1; j < SL->num_v; ++j) {
    double v = fabs(SL->L[j][i]);
    if (v > max) {
      max = v;
      max_i = j;
    }
  }
  //pivoteamento
  if (max_i != i) {
    double *tmp = SL->L[i];
    SL->L[i] = SL->L[max_i];
    SL->L[max_i] = tmp;

    double aux = SL->z[i];
    SL->z[i] = SL->z[max_i];
    SL->z[max_i] = aux;
  }
}

//função para calcular o vetor z
void retrossubpz(SistLinear_t *SL, double * grad) {
  for (int i = 0; i < SL->num_v; ++i) {
    SL->z[i] = grad[i];
    for (int j = 0; j < i; j++)
      SL->z[i] -= SL->L[i][j] * SL->z[j];
  }
}

//função para calcular o delta
void retrossubs2(SistLinear_t *SL, double**hes, double *delta) {
  for (int i = SL->num_v-1; i >=0; --i) {
    delta[i] = SL->z[i];
    for (int j = i+1; j < SL->num_v; j++)
      delta[i] -= hes[i][j] * delta[j];
    delta[i] /= hes[i][i];
    if (isnan(hes[i][i]))
      printf("ERRO: %g\n", hes[i][i]);
  }
}

//função para triangularizar a matriz hessiana
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

//função com a chamada para o procedimento fatoração LU
double * FatLU(SistLinear_t *SL, double *delta, double**hes, double * grad) {
  triangLU(SL, hes, grad);
  retrossubpz(SL, grad);
  retrossubs2(SL, hes, delta);
  return delta;
}

//função principal para o metodo newton Modificado
double ** Newton_Modificado(SistLinear_t *SL, double *TderivadasLU, double *TslLU, double ** m_aux)
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

  //iteração zero
  for(int z = 0; z < SL->num_v; z++)
    m_res[0][z] = SL->Xlu[z];

  int i;

  //for principal
  for (i = 0; i < SL->max_iter; i++)
  {
    double aux = 0.0;

    //atualiza o gradiente
    double * grad = calc_grad(SL, SL->Xlu, TderivadasLU);

    //atualiza a matriz hessiana a cada has steps (numero de variaveis)
    if(i % SL->num_v == 0)
      calc_hes(SL, SL->Xlu, TderivadasLU, m_aux);
  
    for (int i = 0; i < SL->num_v; i++)
      aux += grad[i]*grad[i];
    aux = sqrt(aux);

    //verifica se o gradiente é menor que epsolon
    if(aux < SL->epsilon)
    {
      free(grad);
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

    //calcula o vetor delta e o tempo de execução do calculo
    double tTotal = timestamp();
    delta = FatLU(SL, delta, m_aux, grad);
    *TslLU += timestamp() - tTotal;

    for (int l = 0; l < SL->num_v; l++)
      SL->Xlu[l] += delta[l];

    aux = 0.0;

    for (int i = 0; i < SL->num_v; i++)
      aux += delta[i]*delta[i];
    aux = sqrt(aux);

    //iteração i+1
    for(int z = 0; z < SL->num_v; z++)
      m_res[i+1][z] = SL->Xlu[z];

    free(delta);

    //verifica se delta é menor que epsolon
    if(aux < SL->epsilon)
    {
      free(grad);
      for (int l = i+1; l < SL->max_iter; l++)
        for(int z = 0; z < SL->num_v; z++)
          m_res[l+1][z] = NAN;
      return m_res;
    }
    free(grad);
  }
  return m_res;
}