#ifndef __METODO_DE_NEWTON_INEXATO__
#define __METODO_DE_NEWTON_INEXATO__

#include <stdlib.h>
#include <inttypes.h>
#include "SistLinear.h"

void clean_fgets(char *pos);

void Newton_Inexato(SistLinear_t *e);

void calcula_independentes(SistLinear_t *e, double **matrix_diag);

void calcula_tempo(SistLinear_t *e, double **matriz);

void gaussSeidel(SistLinear_t *e, double **matriz);

#endif