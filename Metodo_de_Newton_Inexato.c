#include "utils.h"
#include "Metodo_de_Newton_Inexato.h"
#include <string.h>
#include <math.h>
#include <matheval.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>
#include <sys/time.h>
#include <stdlib.h>
#include "SistLinear.h"


//função para "limpar" string
void clean_fgets(char *pos) { 
  //strtok(pos, "\n");
}

// Gauss-Seidel 
double * gaussSeidel(SistLinear_t *SL, double *delta, double **m_aux, double *grad){
/*  int i,j,q;
  unsigned int n = e->num_v;
  double *r =  malloc((e->num_v) * sizeof (double)) ;
  double temp,sum,erroMaximo,erroCalculado;
  double **A = matriz;
  double *b = e->b;
  double *x = malloc((e->num_v) * sizeof(double));

  // A[n][n] = Matriz principal (e->f)
  // b[n] = vetor_independente (e->termos_independentes)

  for(i=0;i<n;i++){
      r[i] = 0;
  }
  
  q = 1;
  do{
      erroCalculado = 0;
      q++;
      for(i=0;i<n;i++){
          sum = 0;
          for(j=0;j<n;j++){
              if(i != j){
                  sum = sum + (A[i][j] * r[j]);
              }
          }
          temp = (-1.0 / A[i][i]) * sum + b[i] / A[i][i];
          erroMaximo = fabs(temp - r[i]);
          r[i] = temp;
          if(erroMaximo > erroCalculado){
              erroCalculado = erroMaximo;
          }
      }
  }while(erroCalculado >= e->epsilon && q<=e->max_iter);

  for(i=0;i<n;i++){ //copiar dados calculados para *x
    x[i] = r[i];
  }

  e->Xgs = x; //copiar dados para estrutura

  free(r);*/
}

double * calcula_independentes(SistLinear_t *SL, double **m_aux, double *grad)
{/*
  double *indep = (double *) malloc((SL->num_v) * sizeof(double));
}

double ** Newton_Inexato(SistLinear_t *SL)
{
  double ** m_res = (double**) calloc(2*SL->max_iter+1, sizeof(double*));
  for(int i = 0; i < 2*SL->max_iter+1; i++)
  { 
    m_res[i] = (double*) calloc(SL->num_v, sizeof(double));
  }

  m_res[0] = SL->Xgs;

  for (int i = 0; i < SL->max_iter; i++)
  {
    double aux = 0.0;

    double * grad = calc_grad(SL, SL->Xgs);
    
    double ** m_aux = calc_hes(SL, SL->Xgs);
  
    for (int i = 0; i < SL->num_v; i++)
    {
      if(aux <= fabs(grad[i]))
      {
        aux = grad[i];
      }
    }
  
    if(fabs(aux) < SL->epsilon)
      return m_res;
    double * delta = (double*) calloc(SL->num_v, sizeof(double));

    //calcula_independentes(e, matriz); //Calcula termos independentes
    //gaussSeidel(e, matriz);

    delta = gaussSeidel(SL, delta, m_aux, grad);

    for (int l = 0; l < SL->num_v; l++)
    {
      SL->Xgs[l] += delta[l];
    }

    aux = 0.0;

    for (int i = 0; i < SL->num_v-1; i++)
    {
      if(aux <= delta[i])
      {
        aux = delta[i];
      }
    }

    m_res[i+1] = SL->Xgs;

    if(fabs(aux) < SL->epsilon)
      return m_res;
  }
  return m_res;*/
}