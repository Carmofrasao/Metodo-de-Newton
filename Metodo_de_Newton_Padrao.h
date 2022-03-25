#ifndef __METODO_DE_NEWTON_PADRAO__
#define __METODO_DE_NEWTON_PADRAO__

#include "SistLinear.h"

void pivot(SistLinear_t *SL, double**hes, double * grad, int i);

void retrossubs(SistLinear_t *SL, double *delta, double**hes, double * grad);

void triang(SistLinear_t *SL, double**hes, double * grad);

double* eliminacaoGauss(SistLinear_t *SL, double *delta, double**hes, double * grad) ;

double * calc_grad(SistLinear_t *SL);

double ** calc_hes(SistLinear_t *SL);

double ** Newton_Padrao(SistLinear_t *SL);

#endif