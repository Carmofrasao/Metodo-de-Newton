/* Código escrito por
    Anderson Aparecido do Carmo Frasão 
    GRR 20204069
    Erick Eckermann Cardoso 
    GRR 20186075
*/

#ifndef __METODO_DE_NEWTON_PADRAO__
#define __METODO_DE_NEWTON_PADRAO__

#include "SistLinear.h"

//função para aplicar o pivoteamento na matriz hessiana
void pivot(SistLinear_t *SL, double**hes, double * grad, int i);

//função para calcular o delta
void retrossubs(SistLinear_t *SL, double *delta, double**hes, double * grad);

//função para triangularizar a matriz hessiana
void triang(SistLinear_t *SL, double**hes, double * grad);

//função com a chamada para o procedimento eliminação de gauss
double* eliminacaoGauss(SistLinear_t *SL, double *delta, double**hes, double * grad);

//função principal para o metodo newton padrão
double ** Newton_Padrao(SistLinear_t *SL, double *TderivadasEG, double *TslEG, double ** m_aux);

#endif