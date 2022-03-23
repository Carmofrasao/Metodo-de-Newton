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
  int i,j,q,d;
  unsigned int n = e->n;
  double *r =  malloc((e->n) * sizeof (double)) ;
  double temp,sum,erroMaximo,erroCalculado;
  double **A = matriz;
  double *b = e->b;
  double *x = malloc((e->n) * sizeof(double));
  e->r = malloc((e->n) * sizeof(double));

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

  e->r = x; //copiar dados para estrutura

  free(r);
}

void calcula_independentes(SistLinear_t *e, double **matrix_diag){

  double *indep = (double *) malloc((e->n) * sizeof(double));
  for(int i = 0; i < (e->n); i++){
    indep[i] = matrix_diag[e->n - 1][i];
  }
  e->b = indep;

}

void calcula_tempo(SistLinear_t *e, double **matriz){

  //calculamos tempo antes de executar o gaussSeidel e depois pra fazermos a diferença de ambos
  e->tempo=0;
  e->tempo = timestamp();
  gaussSeidel(e, matriz);
  e->tempo = timestamp() - e->tempo;


}

void Newton_Inexato(SistLinear_t *e){

  // usar libmatheval para gerar vetores com os valor de 0 até n para cada equação
  int i, j, k, l, row, columns;
  double **matriz, linha[e->n], func;

  for(int i = 0; i < e->n; i++)
  {
    linha[i] = 0.0;
  }

  matriz = malloc ((e->n) * sizeof (double*)); //aloca espaco para matrix_diag com valores da diagonal 
  for (i=0; i < (e->n); i++){
    matriz[i] = malloc ((e->n) * sizeof (double));
  }

  for(i=0; i < e->n ; i++){
    for(j=0; j<e->n; j++){
      func = evaluator_evaluate_x(e->f, j);
      linha[j] = func;
    }

    for(k=0; k<e->n; k++){
      matriz[i][k] = linha[k];
    }
  }

  calcula_independentes(e, matriz); //Calcula termos independentes
  calcula_tempo(e, matriz); //calcula tempo de execucao
  
  for (i=0; i < (e->n); i++){
    free(matriz[i]);
  }
  free(matriz);

  free(e->r);

  free(e->b);
} 