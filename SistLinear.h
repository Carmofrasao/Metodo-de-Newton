#ifndef __SIST_LINEAR__
#define __SIST_LINEAR__

#include "utils.h"

/***************************
 ESTRUTURA SISTEMA LINEAR
****************************/
/*
  Método 3: vetor de ponteiros de linhas contíguas
  http://wiki.inf.ufpr.br/maziero/doku.php?id=prog2:alocacao_dinamica_de_matrizes&s[]=aloca%C3%A7%C3%A3o&s[]=de&s[]=matrizes
*/
typedef struct {
  unsigned int n;   // numero de variaveis
  double *M;        // vetor nxn de posições da matriz
  double **A;       // matriz dos coeficientes do SL (vetor de ponteiros para posições de M)
  double *b;        // termos independentes do SL

  double **L;       // matriz dos multiplicadores (LESS)
  double **U;       // matriz dos multiplicadores (UPER)
  double *z;        // vetor z para Fatoração LU

  //char *eq;        //strings com equacao auxiliar pra montagem da matriz
  void **f;       //matriz de strings com as expressoes
  double epsilon;   //tolerancia no metodo de gauss seidel
  int max_iter;     //numero maximo de iterações
  double *r;        //vetor de resultados finais
  double tempo;     //tempo de execucao
  double *ap;        //aproximação inicial

} SistLinear_t;

SistLinear_t *alocaSistLinear(unsigned int n);

void liberaSistLinear(SistLinear_t *SL);

// devolve um outro Sistema Linear que é uma cópia do SL
SistLinear_t *dupSL(SistLinear_t *SL);

// ordem de leitura: n A b
SistLinear_t *lerSistLinear();

#endif