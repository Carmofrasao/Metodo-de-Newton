#ifndef __SIST_LINEAR__
#define __SIST_LINEAR__

#include "utils.h"

typedef struct {
  unsigned int num_v;   // numero de variaveis
  double *M;            // vetor nxn de posições da matriz
  double **A;           // matriz dos coeficientes do SL (vetor de ponteiros para posições de M)
  double *b;            // termos independentes do SL

  double **L;           // matriz dos multiplicadores (LESS)
  double **U;           // matriz dos multiplicadores (UPER)
  double *z;            // vetor z para Fatoração LU

  char *eq_aux;         //strings com equacao auxiliar pra montagem da matriz
  void ***HESSIANA;     //matriz de strings com as expressoes
  void **GRADIENTE;     //vetor de gradientes
  double epsilon;       //tolerancia no metodo de gauss seidel
  int max_iter;         //numero maximo de iterações
  double *Xeg;          //vetor de resultados finais eliminacao
  double *Xlu;          //vetor de resultados finais fat LU
  double *Xgs;          //vetor de resultados finais Gauss

} SistLinear_t;

SistLinear_t *alocaSistLinear(unsigned int n);

void liberaSistLinear(SistLinear_t *SL);

// devolve um outro Sistema Linear que é uma cópia do SL
SistLinear_t *dupSL(SistLinear_t *SL);

SistLinear_t *lerSistLinear();

#endif