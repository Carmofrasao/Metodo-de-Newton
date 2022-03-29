#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <matheval.h>
#include <string.h>
#include <assert.h>

#include "utils.h"
#include "SistLinear.h"

//criando matriz hessiana
void cria_hes(SistLinear_t *SL)
{
  char aux[4];
  char Xn[4];

  clean_fgets(SL->eq_aux);
  for(int n = 0; n < SL->num_v; n++)
  {
    for(int l = 0; l < SL->num_v; l++)
    {
      void * f_aux1;
      void * f_aux2;
      memset(Xn, 0, sizeof(Xn));
      memset(aux, 0, sizeof(aux));
      f_aux1 = evaluator_create(SL->eq_aux);
      assert(f_aux1);
      sprintf(aux, "%d", n+1);
      strcat(strcpy(Xn, "x"), aux); 
      f_aux2 = evaluator_derivative (f_aux1, Xn);
      assert(f_aux2);
      memset(Xn, 0, sizeof(Xn));
      memset(aux, 0, sizeof(aux));
      sprintf(aux, "%d", l+1);
      strcat(strcpy(Xn, "x"), aux); 
      SL->HESSIANA[n][l] = evaluator_derivative (f_aux2, Xn);
      assert(SL->HESSIANA[n][l]);
      evaluator_destroy(f_aux1);
      evaluator_destroy(f_aux2);
    }
  }
}

//criando vetor gradiente
void cria_grad(SistLinear_t *SL)
{
  char aux[4];
  char Xn[4];
    
  for(int l = 0; l < SL->num_v; l++)
  {
    void * f_aux;
    memset(Xn, 0, sizeof(Xn));
    memset(aux, 0, sizeof(aux));
    f_aux = evaluator_create(SL->eq_aux);
    assert(f_aux);
    sprintf(aux, "%d", l+1);
    strcat(strcpy(Xn, "x"), aux); 
    SL->GRADIENTE[l] = evaluator_derivative (f_aux, Xn);
    assert(SL->GRADIENTE[l]);
    evaluator_destroy(f_aux);
  }
}

//função para "limpar" string
void clean_fgets(char *pos) { 
  strtok(pos, "\n");
}

double * calc_grad(SistLinear_t *SL, double * X, double *tempo)
{
  double *res = (double*)malloc(SL->num_v*sizeof(double));
  if (!(res)) {
    free(SL);
    printf("ERRO");
    return 0;
  }
  char **Xs = (char**) malloc(SL->num_v*sizeof(char*));
  if (!(Xs)) {
    free(SL);
    printf("ERRO");
    return 0;
  }
  for (int i = 0; i < SL->num_v; i++)
  {
    Xs[i] = (char*) malloc(5*sizeof(char));
    if (!(Xs[i])) {
      free(SL);
      printf("ERRO");
      return 0;
    }
  }
  
  char aux[4];

  for(int i = 0; i < SL->max_iter; i++)
  {
    for(int l = 0; l < SL->num_v; l++)
    {
      memset(aux, 0, sizeof(aux));
      sprintf(aux, "%d", l+1);
      strcat(strcpy(Xs[l], "x"), aux);
    }
    for(int l = 0; l < SL->num_v; l++)
    {
      double tTotal = timestamp();
      res[l] = -evaluator_evaluate(SL->GRADIENTE[l], SL->num_v, Xs, X);
      *tempo += timestamp() - tTotal;
    }
  }
  for (int i = 0; i < SL->num_v; i++)
  {
    free(Xs[i]);
  }
  free(Xs);
  
  return res;
}

double ** calc_hes(SistLinear_t *SL, double * X, double *tempo)
{
  double ** m_aux = (double**)calloc(SL->num_v, sizeof(double*));
  if (!(m_aux)) {
    free(SL);
    printf("ERRO");
    return 0;
  }
  for(int i = 0; i < SL->num_v; i++)
  {
    m_aux[i] = (double*)calloc(SL->num_v, sizeof(double));
    if (!(m_aux[i])) {
      free(SL);
      printf("ERRO");
      return 0;
    }
  }
  char aux[4];
  
  char **Xs = (char**) calloc(SL->num_v, sizeof(char*));
  if (!(Xs)) {
    free(SL);
    printf("ERRO");
    return 0;
  }
  for(int l = 0; l < SL->num_v; l++)
  {
    Xs[l] = (char*) malloc(5*sizeof(char));
    if (!(Xs[l])) {
      free(SL);
      printf("ERRO");
      return 0;
    }
  }

  for(int l = 0; l < SL->num_v; l++)
  {
    memset(aux, 0, sizeof(aux));
    sprintf(aux, "%d", l+1);
    strcat(strcpy(Xs[l], "x"), aux);
  }

  for (int i = 0; i < SL->num_v; i++)
  {
    for(int l = 0; l < SL->num_v; l++)
    {
      double tTotal = timestamp();
      m_aux[i][l] = evaluator_evaluate(SL->HESSIANA[i][l], SL->num_v, Xs, X);
      *tempo += timestamp() - tTotal;
    }
  }
  for (int i = 0; i < SL->num_v; i++)
  {
    free(Xs[i]);
  }
  free(Xs);
  return m_aux;
}

SistLinear_t *alocaSistLinear(unsigned int n) {

  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));

  if (SL) {
    SL->num_v = n;

    SL->GRADIENTE = (void**) malloc(n * sizeof(void*));
    if (!(SL->GRADIENTE)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }

    SL->HESSIANA = (void***) malloc(n * sizeof(void**));
    if (!(SL->HESSIANA)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }
    for (int i = 0; i < n; i++)
    {
      SL->HESSIANA[i] = (void**) malloc(n * sizeof(void*));
      if (!(SL->HESSIANA[i])) {
        free(SL);
        printf("ERRO");
        return NULL;
      }
    }

    SL->eq_aux = (char*) calloc(1024, sizeof(char));
    if (!(SL->eq_aux)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }

    SL->L = (double**) calloc(n, sizeof(double*));
    if (!(SL->L)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }
    for (int i = 0; i < n; i++)
    {
      SL->L[i] = (double*) calloc(n, sizeof(double));
      if (!(SL->L[i])) {
        free(SL);
        printf("ERRO");
        return NULL;
      }
    }
    
    SL->U = (double**) calloc(n, sizeof(double*));
    if (!(SL->U)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }
    for (int i = 0; i < n; i++)
    {
      SL->U[i] = (double*) calloc(n, sizeof(double));
      if (!(SL->U[i])) {
        free(SL);
        printf("ERRO");
        return NULL;
      }
    }

    SL->z = (double*) calloc(n, sizeof(double));
    if (!(SL->z)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }
  
    SL->Xeg = (double*) calloc(n, sizeof(double));
    if (!(SL->Xeg)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }

    SL->Xlu = (double*) calloc(n, sizeof(double));
    if (!(SL->Xlu)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }

    SL->Xgs = (double*) calloc(n, sizeof(double));
    if (!(SL->Xgs)) {
      free(SL);
      printf("ERRO");
      return NULL;
    }
  }

  return SL;
}

void liberaSistLinear(SistLinear_t *SL) {
  for(int i = 0; i < SL->num_v; i++)
    free(SL->L[i]);
  free(SL->L);
  for(int i = 0; i < SL->num_v; i++)
    free(SL->U[i]);
  free(SL->U);
  for (int i = 0; i < SL->num_v; i++)
    for (int z = 0; z < SL->num_v; z++)
      evaluator_destroy(SL->HESSIANA[i][z]);
  for (int i = 0; i < SL->num_v; i++)
    evaluator_destroy(SL->GRADIENTE[i]);
  free(SL->eq_aux);
  free(SL->z);
  free(SL->Xeg);
  free(SL->Xlu);
  free(SL->Xgs);
  free(SL);
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