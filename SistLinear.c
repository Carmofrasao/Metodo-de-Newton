#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"

SistLinear_t *alocaSistLinear(unsigned int n) {

  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));

  if (SL) {
    SL->n = n;

    SL->A = (double **) malloc(sizeof(double *)*n);
    if (!(SL->A)) {
      free(SL);
      return NULL;
    }

    SL->M = (double *) malloc(sizeof(double)*n*n);
    if (!(SL->M)) {
      free(SL->A);
      free(SL);
      return NULL;
    }

    for (int i=0; i < n; ++i)
      SL->A[i] = SL->M + i*n;

    SL->b = (double *) malloc(sizeof(double)*n);
    if (!(SL->b)) {
      free(SL->M);
      free(SL->A);
      free(SL);
      return NULL;
    }
  }

  return SL;
}

void liberaSistLinear(SistLinear_t *SL) {
  free(SL->b);
  free(SL->M);
  free(SL->A);
  free(SL);
}

SistLinear_t *dupSL(SistLinear_t *SL) {
  SistLinear_t *dup = alocaSistLinear(SL->n);

  if (dup) {
    for(int i = 0; i < SL->n; ++i) {
      for(int j = 0; j < SL->n; ++j)
        dup->A[i][j] = SL->A[i][j];
      dup->b[i] = SL->b[i];
    }
  }
  
  return dup;
}

SistLinear_t *lerSistLinear() {

  unsigned int n;

  SistLinear_t *SL = NULL;
  
  if (scanf("%d",&n) != EOF) {  
    
    SL = alocaSistLinear(n);
    if (!SL) return NULL;
    
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
	      scanf("%lf", &(SL->A[i][j]));
    
    for(int i = 0; i < n; ++i)
      scanf("%lf", &(SL->b[i]));
  }
  
  return SL;
}