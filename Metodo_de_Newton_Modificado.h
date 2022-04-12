/* Código escrito por
    Anderson Aparecido do Carmo Frasão 
    GRR 20204069
    Erick Eckermann Cardoso 
    GRR 20186075
*/

#ifndef __METODO_DE_NEWTON_MODIFICADO__
#define __METODO_DE_NEWTON_MODIFICADO__

#include "SistLinear.h"

//função para aplicar o pivoteamento na matriz hessiana
void pivotLU(SistLinear_t *SL, int i, double**hes, double * grad);

//função para calcular o vetor z
void retrossubpz(SistLinear_t *SL, double * grad);

//função para calcular o delta
void retrossubs2(SistLinear_t *SL, double**hes, double *delta);

//função para triangularizar a matriz hessiana
void triangLU(SistLinear_t *SL, double**hes, double * grad);

//função com a chamada para o procedimento fatoração LU
double * FatLU(SistLinear_t *SL, double *delta, double**hes, double * grad);

//função principal para o metodo newton Modificado
double ** Newton_Modificado(SistLinear_t *SL, double *TderivadasLU, double *TslLU, double ** m_aux);

#endif