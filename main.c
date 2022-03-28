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
  
  int i = 1;
  char aux[4];
  char Xn[4];
  void * f_aux;

  double TtotalEG, TtotalLU, TtotalGS, TderivadasEG, TderivadasLU, TderivadasGS, TslEG, TslLU, TslGS;

  TtotalEG = 0.0;
  TtotalLU = 0.0;
  TtotalGS = 0.0;
  TderivadasEG = 0.0;
  TderivadasLU = 0.0;
  TderivadasGS = 0.0;
  TslEG = 0.0;
  TslLU = 0.0;
  TslGS = 0.0;
  
  while (SL = lerSistLinear())
  {
    double ** m_reseg = (double**) calloc(SL->max_iter+1, sizeof(double*));
    if (!(m_reseg)){
      free(SL);
      printf("ERRO");
      return 0;
    }
    for(int i = 0; i < SL->max_iter+1; i++)
    { 
      m_reseg[i] = (double*) calloc(SL->num_v, sizeof(double));
      if (!(m_reseg[i])){
        free(SL);
        printf("ERRO");
        return 0;
      }
    }

    double ** m_reslu = (double**) calloc(SL->max_iter+1, sizeof(double*));
    if (!(m_reslu)){
      free(SL);
      printf("ERRO");
      return 0;
    }
    for(int i = 0; i < SL->max_iter+1; i++)
    { 
      m_reslu[i] = (double*) calloc(SL->num_v, sizeof(double));
      if (!(m_reslu[i])){
        free(SL);
        printf("ERRO");
        return 0;
      }
    }

    double ** m_resgs = (double**) calloc(SL->max_iter+1, sizeof(double*));
    if (!(m_resgs)){
      free(SL);
      printf("ERRO");
      return 0;
    }
    for(int i = 0; i < SL->max_iter+1; i++)
    { 
      m_resgs[i] = (double*) calloc(SL->num_v, sizeof(double));
      if (!(m_resgs[i])){
        free(SL);
        printf("ERRO");
        return 0;
      }
    }

    cria_hes(SL);

    cria_grad(SL);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //Newton Padrão

    double tTotal = timestamp();
    m_reseg = Newton_Padrao(SL, &TderivadasEG, &TslEG);
    TtotalEG = timestamp() - tTotal;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Modificado

    tTotal = timestamp();
    m_reslu = Newton_Modificado(SL, &TderivadasLU, &TslLU);
    TtotalLU = timestamp() - tTotal;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Inexato

    tTotal = timestamp();
    m_resgs = Newton_Inexato(SL, &TderivadasGS, &TslGS);
    TtotalGS = timestamp() - tTotal;
  
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
      final[0] = NAN;
      final[1] = NAN;
      final[2] = NAN;
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
    
    ++i;
    printf("\n");
  }
  liberaSistLinear(SL);
  evaluator_destroy(f_aux);
}