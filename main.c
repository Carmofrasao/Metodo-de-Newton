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
    clean_fgets(SL->eq);
    
    char aux[2];
    printf("%d %s\n", SL->n, SL->eq);
    for(int i = 0; i < SL->n; i++)
    {
      printf("%f ", SL->ap[i]);
    }
    printf("\n%e %d\n", SL->epsilon, SL->max_iter);
    for(int n = 0; n < SL->n; n++)
    {
      for(int l = 0; l < SL->n; l++)
      {
        f = evaluator_create(SL->eq);
        assert(f);
        char buffer[2];
        sprintf(aux, "%d", n);
        strcat(strcpy(buffer, "x"), aux); 
        f = evaluator_derivative (f, buffer);
        assert(f);
        printf("%s ", buffer);
        sprintf(aux, "%d", l);
        strcat(strcpy(buffer, "x"), aux); 
        f = evaluator_derivative (f, buffer);
        assert(f);
        printf("%s\n", buffer);
        SL->f[n][l] = f;
      }
    }
    printf("rodando %d\n\n", i);
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