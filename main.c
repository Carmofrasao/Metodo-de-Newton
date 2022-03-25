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
  
  int i = 0;
  void *f_aux;
  char aux[4];
  char Xn[4];

  double TtotalEG, TtotalLU, TtotalGS, TderivadasEG, TderivadasLU, TderivadasGS, TslEG, TslLU, TslGS;
  
  while (SL = lerSistLinear())
  {
    //criando matriz hessiana
    clean_fgets(SL->eq_aux);
    
    for(int n = 0; n < SL->num_v; n++)
    {
      for(int l = 0; l < SL->num_v; l++)
      {
        memset(Xn, 0, sizeof(Xn));
        memset(aux, 0, sizeof(aux));
        f_aux = evaluator_create(SL->eq_aux);
        assert(f_aux);
        printf("Funcao original\n");
        printf("%s\n", evaluator_get_string (f_aux));
        sprintf(aux, "%d", n+1);
        strcat(strcpy(Xn, "x"), aux); 
        f_aux = evaluator_derivative (f_aux, Xn);
        assert(f_aux);
        printf("Primeira derivada\n");
        printf("%s\n", evaluator_get_string (f_aux));
        memset(Xn, 0, sizeof(Xn));
        memset(aux, 0, sizeof(aux));
        sprintf(aux, "%d", l+1);
        strcat(strcpy(Xn, "x"), aux); 
        f_aux = evaluator_derivative (f_aux, Xn);
        assert(f_aux);
        printf("Segunda derivada\n");
        printf("%s\n", evaluator_get_string (f_aux));
        SL->HESSIANA[n][l] = f_aux;
      }
    }
    printf("\n");

    //criando vetor gradiente
    for(int l = 0; l < SL->num_v; l++)
    {
      memset(Xn, 0, sizeof(Xn));
      memset(aux, 0, sizeof(aux));
      f_aux = evaluator_create(SL->eq_aux);
      assert(f_aux);
      sprintf(aux, "%d", l+1);
      strcat(strcpy(Xn, "x"), aux); 
      f_aux = evaluator_derivative (f_aux, Xn);
      assert(f_aux);
      SL->GRADIENTE[l] = f_aux;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Metodo de Newton Padrao

    double tTotal = timestamp();
    SL->X = Newton_Padrao(SL);
    TtotalEG = timestamp() - tTotal;

    for (int i = 0; i < SL->num_v; i++)
    {
      printf("%f ", SL->X[i]);
    }
    printf("\n\n");
    
    
    /*
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Modificado

    memset(SL->X, 0, sizeof(SL->X));

    tTotal = timestamp();
    FatLU(SL);
    TtotalLU = timestamp() - tTotal;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Inexato

    memset(SL->X, 0, sizeof(SL->X));

    tTotal = timestamp();
    Newton_Inexato(SL);
    TtotalGS = timestamp() - tTotal;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // cabeçalho
    printf("Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");

    // para cada iteração
    for (...) {
      printf("%d \t\t| ", i); // imprime iteração

      if (...) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(fx) || isinf(fx))
          printf("%1.14e\t\t\t| ", fx);
        else
          printf("%1.14e\t| ", fx);
      }
      else
        printf("\t\t\t| ");

      // repete para as outras duas colunas...
    }

    // imprimir os tempos
    printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e\n", TtotalEG, TtotalLU, TtotalGS);
    printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e\n", TderivadasEG, TderivadasLU, TderivadasGS);
    printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e\n", TslEG, TslLU, TslGS);
    */
    ++i;
  }
  //liberaSistLinear(SL);
  evaluator_destroy(f_aux);
}