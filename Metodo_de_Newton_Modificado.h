#ifndef __METODO_DE_NEWTON_MODIFICADO__
#define __METODO_DE_NEWTON_MODIFICADO__

#include "SistLinear.h"

void pivotLU(SistLinear_t *SL, int i, double**hes, double * grad);

void pivotz(SistLinear_t *SL, int i);

void retrossubpz(SistLinear_t *SL);

void retrossubs2(SistLinear_t *SL, double**hes, double *delta);

void triangLU(SistLinear_t *SL, double**hes, double * grad);

double * FatLU(SistLinear_t *SL, double *delta, double**hes, double * grad);

double ** Newton_Modificado(SistLinear_t *SL);

#endif