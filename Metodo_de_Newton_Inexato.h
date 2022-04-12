/*Alunos:
Anderson Aparecido do Carmo Frasão
GRR 20204069

Erick Eckermann Cardoso
GRR20186075
*/

#ifndef __METODO_DE_NEWTON_INEXATO__
#define __METODO_DE_NEWTON_INEXATO__

#include "SistLinear.h"

//Função utilizada em Gauss-Seidel que recebe Xi e retorna Xi + 1
double * calcula_independentes(SistLinear_t *SL, double **m_aux, double *grad, double * res);

// Método Gauss-Seidel para resolução de sistemas lineares
double * gaussSeidel(SistLinear_t *SL, double **m_aux, double *grad);

//função principal para o metodo newton Inexato
double ** Newton_Inexato(SistLinear_t *SL, double *TderivadasGS, double * TlsGS, double ** m_aux);

#endif