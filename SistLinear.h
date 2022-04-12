/* Código escrito por
    Anderson Aparecido do Carmo Frasão 
    GRR 20204069
    Erick Eckermann Cardoso 
    GRR 20186075
*/

#ifndef __SIST_LINEAR__
#define __SIST_LINEAR__

#include "utils.h"
#include <assert.h>

typedef struct {
  unsigned int num_v;   // numero de variaveis

  double **L;           // matriz dos multiplicadores (LESS)
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

//criando matriz hessiana a partir da função original
void cria_hes(SistLinear_t *SL);

//criando vetor gradiente a partir da função original
void cria_grad(SistLinear_t *SL);

//calcula o vetor gradiente
double * calc_grad(SistLinear_t *SL, double * X, double *tempo);

//calculando a matriz hessiana
void calc_hes(SistLinear_t *SL, double * X, double *tempo, double ** m_aux);

//alocando espaço para o sistema linear
SistLinear_t *alocaSistLinear(unsigned int n);

//liberando a alocação de memoria
void liberaSistLinear(SistLinear_t *SL);

//lendo o sistema linear
SistLinear_t *lerSistLinear();

#endif