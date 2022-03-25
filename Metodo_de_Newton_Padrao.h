#ifndef __METODO_DE_NEWTON_PADRAO__
#define __METODO_DE_NEWTON_PADRAO__

#include "SistLinear.h"

void pivot(SistLinear_t *SL, int i);

void retrossubs(SistLinear_t *SL);

void triang(SistLinear_t *SL);

void eliminacaoGauss(SistLinear_t *SL);

#endif