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
  void *f_aux;
  char aux[4];
  char Xn[4];

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
        sprintf(aux, "%d", n+1);
        strcat(strcpy(Xn, "x"), aux); 
        f_aux = evaluator_derivative (f_aux, Xn);
        assert(f_aux);
        memset(Xn, 0, sizeof(Xn));
        memset(aux, 0, sizeof(aux));
        sprintf(aux, "%d", l+1);
        strcat(strcpy(Xn, "x"), aux); 
        f_aux = evaluator_derivative (f_aux, Xn);
        assert(f_aux);
        SL->HESSIANA[n][l] = f_aux;
      }
    }

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
    
    //Newton Padrão

    double tTotal = timestamp();
    SL->Xeg = Newton_Padrao(SL);
    TtotalEG = timestamp() - tTotal;
    
    /*
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Modificado

    tTotal = timestamp();
    FatLU(SL);
    TtotalLU = timestamp() - tTotal;

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Inexato

    memset(SL->X, 0, sizeof(SL->X));

    tTotal = timestamp();
    Newton_Inexato(SL);
    TtotalGS = timestamp() - tTotal;
    */
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // cabeçalho
    printf("%d\n", SL->num_v);
    printf("%s\n", SL->eq_aux);
    printf("#Iteração \t| Newton Padrão \t| Newton Modificado \t| Newton Inexato\n");

    // para cada iteração
    for (int i = 0; i < SL->max_iter; i++) {
      printf("%d \t\t| ", i+1); // imprime iteração

      if (SL->Xeg[i] != 0.0) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(SL->Xeg[i]) || isinf(SL->Xeg[i]))
          printf("%1.14e\t\t\t| ", SL->Xeg[i]);
        else
          printf("%1.14e\t| ", SL->Xeg[i]);
      }
      else
        printf("\t\t\t| ");

      if (SL->Xlu[i] != 0.0) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(SL->Xlu[i]) || isinf(SL->Xlu[i]))
          printf("%1.14e\t\t\t| ", SL->Xlu[i]);
        else
          printf("%1.14e\t| ", SL->Xlu[i]);
      }
      else
        printf("\t\t\t| ");

      if (SL->Xgs[i] != 0.0) {  // se nesta iteração o valor da primeira coluna existe, imprime
        if (isnan(SL->Xgs[i]) || isinf(SL->Xgs[i]))
          printf("%1.14e\t\t\t ", SL->Xgs[i]);
        else
          printf("%1.14e\t ", SL->Xgs[i]);
      }
      else
        printf("\t\t\t ");
      printf("\n");
    }

    // imprimir os tempos
    printf("Tempo total \t| %1.14e\t| %1.14e\t| %1.14e\n", TtotalEG, TtotalLU, TtotalGS);
    printf("Tempo derivadas | %1.14e\t| %1.14e\t| %1.14e\n", TderivadasEG, TderivadasLU, TderivadasGS);
    printf("Tempo SL \t| %1.14e\t| %1.14e\t| %1.14e\n", TslEG, TslLU, TslGS);
    
    ++i;
    printf("\n");
  }
  //liberaSistLinear(SL);
  evaluator_destroy(f_aux);
}