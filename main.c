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
  
  while (SL = lerSistLinear())
  {

    //criando matriz hessiana
    clean_fgets(SL->eq_aux);
    
    char aux[4];
    char Xn[4];
    
    for(int n = 0; n < SL->num_v; n++)
    {
      for(int l = 0; l < SL->num_v; l++)
      {
        memset(Xn, 0, sizeof(Xn));
        memset(aux, 0, sizeof(aux));
        f_aux = evaluator_create(SL->eq_aux);
        assert(f_aux);
        sprintf(aux, "%d", n);
        strcat(strcpy(Xn, "x"), aux); 
        f_aux = evaluator_derivative (f_aux, Xn);
        assert(f_aux);
        memset(Xn, 0, sizeof(Xn));
        memset(aux, 0, sizeof(aux));
        sprintf(aux, "%d", l);
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
      sprintf(aux, "%d", l);
      strcat(strcpy(Xn, "x"), aux); 
      f_aux = evaluator_derivative (f_aux, Xn);
      assert(f_aux);
      SL->GRADIENTE[l] = f_aux;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*
    // Metodo de Newton Padrao
    double *X = (double *) malloc(sizeof(double)*SL->n);

    printf("\n***** Sistema %d --> n = %d\n", i+1, SL->n);
    double tTotal = timestamp();
    eliminacaoGauss(SL, X);
    tTotal = timestamp() - tTotal;
    printf("  --> X: ");
    prnVetorDouble(X, SL->n);
    printf("  --> Tempo: %lf ms\n", tTotal);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Modificado

    printf("\n***** Sistema %d --> n = %d\n", i+1, SL->n);

    tTotal = timestamp();
    FatLU(SL, X);
    tTotal = timestamp() - tTotal;
    printf("  --> X: ");
    prnVetorDouble(X, SL->n);
    printf("  --> Tempo: %lf ms\n", tTotal);
    

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //Metodo de Newton Inexato

    Newton_Inexato(SL);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    free(X);*/
    ++i;
  }
  //liberaSistLinear(SL);
  evaluator_destroy(f_aux);
}