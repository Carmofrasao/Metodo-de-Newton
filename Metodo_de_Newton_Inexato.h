#ifndef __METODO_DE_NEWTON_INEXATO__
#define __METODO_DE_NEWTON_INEXATO__

#include "SistLinear.h"

double * calcula_independentes(SistLinear_t *SL, double **m_aux, double *grad, double * res);

double * gaussSeidel(SistLinear_t *SL, double **m_aux, double *grad);

//função principal para o metodo newton Inexato
double ** Newton_Inexato(SistLinear_t *SL, double *TderivadasGS, double * TlsGS, double ** m_aux);

#endif