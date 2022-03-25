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
  strtok(pos, "\n");
}

// Gauss-Seidel 
void gaussSeidel(SistLinear_t *e, double **matriz) {
  int i,j,q;
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

  free(r);
}

void calcula_independentes(SistLinear_t *e, double **matrix_diag){

  double *indep = (double *) malloc((e->num_v) * sizeof(double));
  for(int i = 0; i < (e->num_v); i++){
    indep[i] = matrix_diag[e->num_v - 1][i];
  }
  e->b = indep;

}

void Newton_Inexato(SistLinear_t *e){

  // usar libmatheval para gerar vetores com os valor de 0 até n para cada equação
  int i, j, k;
  double **matriz, linha[e->num_v], func;

  for(int i = 0; i < e->num_v; i++)
  {
    linha[i] = 0.0;
  }

  matriz = malloc ((e->num_v) * sizeof (double*)); //aloca espaco para matrix_diag com valores da diagonal 
  for (i=0; i < (e->num_v); i++){
    matriz[i] = malloc ((e->num_v) * sizeof (double));
  }

  for(i=0; i < e->num_v ; i++){
    for(j=0; j<e->num_v; j++){
      func = evaluator_evaluate_x(e->HESSIANA, j);
      linha[j] = func;
    }

    for(k=0; k<e->num_v; k++){
      matriz[i][k] = linha[k];
    }
  }

  calcula_independentes(e, matriz); //Calcula termos independentes
  gaussSeidel(e, matriz);
  
  for (i=0; i < (e->num_v); i++){
    free(matriz[i]);
  }
  free(matriz);

  free(e->Xgs);

  free(e->b);
} 