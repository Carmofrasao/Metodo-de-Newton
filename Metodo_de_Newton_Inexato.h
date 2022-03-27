#ifndef __METODO_DE_NEWTON_INEXATO__
#define __METODO_DE_NEWTON_INEXATO__

#include <stdlib.h>
#include <inttypes.h>
#include "SistLinear.h"

void clean_fgets(char *pos);

double ** Newton_Inexato(SistLinear_t *SL);

void calcula_independentes(SistLinear_t *e, double **matrix_diag);

void calcula_tempo(SistLinear_t *e, double **matriz);

double * gaussSeidel(SistLinear_t *SL, double *delta, double **m_aux, double *grad);

#endif