#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "SistLinear.h"

SistLinear_t *alocaSistLinear(unsigned int n) {

  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));

  if (SL) {
    SL->num_v = n;

    SL->GRADIENTE = (void**) calloc(n, sizeof(void*));
    if (!(SL->GRADIENTE)) {
      free(SL);
      return NULL;
    }

    SL->HESSIANA = (void***) calloc(n, sizeof(void**));
    if (!(SL->HESSIANA)) {
      free(SL);
      return NULL;
    }
    for (int i = 0; i < n; i++)
    {
      SL->HESSIANA[i] = (void**) calloc(n, sizeof(void*));
    }

    SL->eq_aux = (char*) calloc(1024, sizeof(char));
    if (!(SL->eq_aux)) {
      free(SL);
      return NULL;
    }

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

    SL->L = (double**) malloc(n*sizeof(double*));
    if (!(SL->L)) {
      free(SL->M);
      free(SL->A);
      free(SL);
      return NULL;
    }
    for (int i = 0; i < n; i++)
    {
      SL->L[i] = (double*) malloc(n*sizeof(double));
    }
    
    SL->U = (double**) malloc(n*sizeof(double*));
    if (!(SL->U)) {
      free(SL->M);
      free(SL->A);
      free(SL);
      return NULL;
    }
    for (int i = 0; i < n; i++)
    {
      SL->U[i] = (double*) malloc(n*sizeof(double));
    }

    SL->z = (double*) calloc(n, sizeof(double));
    if (!(SL->z)) {
      free(SL);
      return NULL;
    }
  
    SL->Xeg = (double*) calloc(n, sizeof(double));
    if (!(SL->Xeg)) {
      free(SL);
      return NULL;
    }

    SL->Xlu = (double*) calloc(n, sizeof(double));
    if (!(SL->Xlu)) {
      free(SL);
      return NULL;
    }

    SL->Xgs = (double*) calloc(n, sizeof(double));
    if (!(SL->Xgs)) {
      free(SL);
      return NULL;
    }
  }

  return SL;
}

void liberaSistLinear(SistLinear_t *SL) {
  for(int i = 0; i < SL->num_v; i++)
  {
    free(SL->L[i]);
  }
  free(SL->L);
  for(int i = 0; i < SL->num_v; i++)
  {
    free(SL->U[i]);
  }
  free(SL->U);
  free(SL->z);
  free(SL->Xeg);
  free(SL->Xlu);
  free(SL->Xgs);
  free(SL->b);
  free(SL->M);
  free(SL->A);
  for(int i = 0; i < SL->num_v; i++)
  {
    free(SL->A[i]);
  }
  free(SL);
}

SistLinear_t *dupSL(SistLinear_t *SL) {
  SistLinear_t *dup = alocaSistLinear(SL->num_v);

  if (dup) {
    for(int i = 0; i < SL->num_v; ++i) {
      for(int j = 0; j < SL->num_v; ++j)
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
    scanf("%s", SL->eq_aux);
    double aux = 0.0;
    for (int l = 0; l < SL->num_v; l++)
    {
      scanf("%le", &(aux));
      SL->Xeg[l] = aux;
      SL->Xlu[l] = aux;
      SL->Xgs[l] = aux;
    }
    
    scanf("%le", &(SL->epsilon));
    scanf("%i", &(SL->max_iter));
  }

  return SL;
}