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
  void *f;
  
  while (SL = lerSistLinear())
  {
    clean_fgets(SL->eq_aux);
    
    char aux[4];
    char buffer[4];
    
    for(int n = 0; n < SL->num_v; n++)
    {
      for(int l = 0; l < SL->num_v; l++)
      {
        memset(buffer, 0, sizeof(buffer));
        memset(aux, 0, sizeof(aux));
        f = evaluator_create(SL->eq_aux);
        assert(f);
        sprintf(aux, "%d", n);
        strcat(strcpy(buffer, "x"), aux); 
        f = evaluator_derivative (f, buffer);
        assert(f);
        memset(buffer, 0, sizeof(buffer));
        memset(aux, 0, sizeof(aux));
        sprintf(aux, "%d", l);
        strcat(strcpy(buffer, "x"), aux); 
        f = evaluator_derivative (f, buffer);
        assert(f);
        SL->HESSIANA[n][l] = f;
      }
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
  evaluator_destroy(f);
}