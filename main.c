/* Código escrito por
    Anderson Aparecido do Carmo Frasão 
    GRR 20204069
    Erick Eckermann Cardoso 
    GRR 20186075
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "SistLinear.h"
#include <matheval.h>
#include "utils.h"
#include <inttypes.h>
#include <assert.h>
#include "Metodo_de_Newton_Inexato.h"
#include "Metodo_de_Newton_Modificado.h"
#include "Metodo_de_Newton_Padrao.h"


int main (){
  SistLinear_t *SL;
  char aux[4];
  char Xn[4];

  double TtotalEG = 0;
  double TtotalLU= 0;
  double TtotalGS = 0;
  double TderivadasEG = 0; 
  double TderivadasLU = 0;
  double TderivadasGS = 0;
  double TslEG = 0;
  double TslLU = 0;
  double TslGS = 0;
  
  while (SL = lerSistLinear())
  {
    double ** m_aux = (double**)calloc(SL->num_v, sizeof(double*));
    if (!(m_aux)) {
      free(SL);
      printf("ERRO");
      return 0;
    }
    for(int i = 0; i < SL->num_v; i++)
    {
      m_aux[i] = (double*)calloc(SL->num_v, sizeof(double));
      if (!(m_aux[i])) {
        free(SL);
        printf("ERRO");
        return 0;
      }
    }

    void * f_aux;
    double ** m_reseg;
    double ** m_reslu;
    double ** m_resgs;

    cria_hes(SL);
    cria_grad(SL);
    
    //Newton Padrão

    double tTotal = timestamp();
    m_reseg = Newton_Padrao(SL, &TderivadasEG, &TslEG, m_aux);
    TtotalEG = timestamp() - tTotal;    //calculando o tempo total do newton padrao

    //Metodo de Newton Modificado

    tTotal = timestamp();
    m_reslu = Newton_Modificado(SL, &TderivadasLU, &TslLU, m_aux);
    TtotalLU = timestamp() - tTotal;    //calculando o tempo total do newton modificado

    //Metodo de Newton Inexato

    tTotal = timestamp();
    m_resgs = Newton_Inexato(SL, &TderivadasGS, &TslGS, m_aux);
    TtotalGS = timestamp() - tTotal;    ////calculando o tempo total do newton inexato

    //inicio do processamento de impressão

    memset(Xn, 0, sizeof(Xn));
    memset(aux, 0, sizeof(aux));
    f_aux = evaluator_create(SL->eq_aux);
    assert(f_aux);

    char *X[SL->num_v];

    for (int i = 0; i < SL->num_v; i++)
    {
      char* Xi = (char*) malloc(5*sizeof(char));
      if (!(Xi)){
        free(SL);
        printf("ERRO");
        return 0;
      }
      sprintf(aux, "%d", i+1);
      strcat(strcpy(Xi, "x"), aux);
      X[i] = Xi;
    }

    // cabeçalho
    printf("%d\n", SL->num_v);
    printf("%s\n", SL->eq_aux);
    printf("#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");
    double final[3];
    // para cada iteração
    for (int i = 0; i <= SL->max_iter; i++) {
      final[0] = final[1] = final[2] = NAN;

      final[0] = evaluator_evaluate (f_aux, SL->num_v, X, m_reseg[i]);
      final[1] = evaluator_evaluate (f_aux, SL->num_v, X, m_reslu[i]);
      final[2] = evaluator_evaluate (f_aux, SL->num_v, X, m_resgs[i]);
      if (isnan(final[0]) && isnan(final[1]) && isnan(final[2]))
        break;
      
      printf("%d \t\t| ", i); // imprime iteração
      
      if (final[0] != NAN) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(final[0]) || isinf(final[0]))
          printf("\t\t\t| ");
        else
          printf("%1.14e\t| ", final[0]);
      }
      else
        printf("\t\t\t| ");

      if (final[1] != NAN) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(final[1]) || isinf(final[1]))
          printf("\t\t\t| ");
        else
          printf("%1.14e\t| ", final[1]);
      }
      else
        printf("\t\t\t| ");

      if (final[2] != NAN) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(final[2]) || isinf(final[2]))
          printf("\t\t\t ");
        else
          printf("%1.14e\t ", final[2]);
      }
      else
        printf("\t\t\t ");
      printf("\n");
    }

    // imprimir os tempos
    printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e\n", TtotalEG, TtotalLU, TtotalGS);
    printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e\n", TderivadasEG, TderivadasLU, TderivadasGS);
    printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e\n", TslEG, TslLU, TslGS);
    printf("#\n");
    printf("\n");
    
    for (int i = 0; i < SL->num_v; i++)
      free(X[i]);
    
    for(int i = 0; i < SL->max_iter+1; i++)
    { 
      free(m_reseg[i]);
      free(m_reslu[i]);
      free(m_resgs[i]);
    }
    free(m_reseg);
    free(m_reslu);
    free(m_resgs);
    for (int i = 0; i < SL->num_v; i++)
      free(m_aux[i]);
    free(m_aux);
    liberaSistLinear(SL);
    evaluator_destroy(f_aux);
  }
}