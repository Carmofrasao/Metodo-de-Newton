#ifndef __METODO_DE_NEWTON_MODIFICADO__
#define __METODO_DE_NEWTON_MODIFICADO__

#include "SistLinear.h"

void pivotLU(SistLinear_t *SL, int i);

void pivotz(SistLinear_t *SL, int i);

void retrossubpz(SistLinear_t *SL);

void retrossubs2(SistLinear_t *SL, double *X);

void triangLU(SistLinear_t *SL);

void FatLU(SistLinear_t *SL, double *X);

#endif