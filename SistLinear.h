#ifndef __SIST_LINEAR__
#define __SIST_LINEAR__

#include "utils.h"
#include <assert.h>

typedef struct {
  unsigned int num_v;   // numero de variaveis

  double **L;           // matriz dos multiplicadores (LESS)
  double **U;           // matriz dos multiplicadores (UPER)
  double *z;            // vetor z para Fatoração LU

  char *eq_aux;         //strings com equacao auxiliar pra montagem da matriz
  void ***HESSIANA;     //matriz de derivadas
  void **GRADIENTE;     //vetor de gradientes
  double epsilon;       //tolerancia no metodo de gauss seidel
  int max_iter;         //numero maximo de iterações
  double *Xeg;          //vetor de resultados finais eliminacao
  double *Xlu;          //vetor de resultados finais fat LU
  double *Xgs;          //vetor de resultados finais Gauss

} SistLinear_t;

void cria_hes(SistLinear_t *SL);

void cria_grad(SistLinear_t *SL);

void clean_fgets(char *pos);

double * calc_grad(SistLinear_t *SL, double * X, double *tempo);

void calc_hes(SistLinear_t *SL, double * X, double *tempo, double ** m_aux);

SistLinear_t *alocaSistLinear(unsigned int n);

void liberaSistLinear(SistLinear_t *SL);

SistLinear_t *lerSistLinear();

#endif