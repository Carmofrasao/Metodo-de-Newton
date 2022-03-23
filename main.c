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

  SistLinear_t aux;

  SistLinear_t *SL = &aux;
  int i = 0;
  void *f;

  scanf("%d", &(SL->n));

  char *eq = (char*) calloc(SL->n, sizeof(char)); 
  
  scanf("%s", eq);

  while (eq != NULL)
  {
    //leitura das variáveis a partir de um arquivo 
    printf("chegou aqui");

    SL->ap = (double*) calloc(SL->n, sizeof(double));
    clean_fgets(eq);
    f = evaluator_create(eq);
    assert(f);
    SL->f[i] = f;

    for(int m = 0; m < SL->n; m++)
      scanf("%le", &SL->ap[m]);
    scanf("%le", &(SL->epsilon));
    scanf("%i", &(SL->max_iter));


    /////////////////////////////////////***************************************************************************************************

    // Metodo de Newton Padrao
    double *X = (double *) malloc(sizeof(double)*SL->n);

    printf("\n***** Sistema %d --> n = %d\n", i+1, SL->n);
    printf("chegou aqui");
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

    ++i;
    free(X);
    fscanf(stdin, "%d", &(SL->n)); 

    eq = (char*) calloc(SL->n, sizeof(char));
  
    fscanf(stdin, "%s", eq);
  }
  
  liberaSistLinear(SL);
  evaluator_destroy(f);
}