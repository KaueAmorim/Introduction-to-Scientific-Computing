#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "edo.h"
#include "gaussSeidel_EqDiff.h"

#define MAXIT 50

real_t pp (real_t x);
real_t qq (real_t x);
real_t rr (real_t x);

int main ()
{
  rtime_t tTotal1, tTotal2;
  int tamanhos[] = {5, 10, 100, 1000};
  int num_tamanhos = 4;
  
  for (int k = 0; k < num_tamanhos; k++) {
    int n = tamanhos[k];
    printf("\n========================================\n");
    printf("RESULTADOS PARA n = %d\n", n);
    printf("========================================\n");
    
    EDo edo = {n, 0, 1, -1, 0, pp, qq, rr}; // definição da EDO do exercício

    printf("\nSistema Linear Resultante (Diagonais e Termos Independentes):\n");
    prnEDOsl(&edo, 1);

    // VERSÃO 1: Usando vetores pré-calculados para diagonais e termos independentes
    printf("\n--- VERSÃO 1: Gauss-Seidel com vetores pré-calculados ---\n");
    real_t *Y1 = (real_t *) calloc(edo.n, sizeof(real_t)); // Resultado da EDO.
    
    Tridiag *SL = genTridiag(&edo);
    tTotal1 = gaussSeidel_3Diag(SL, Y1, MAXIT);

    printf("\nSolucao (Y1):\n");
    prnVetor(Y1, edo.n);

    real_t norma1 = normaL2_residuo_3Diag(SL, Y1);
    printf("\nNorma L2 do Residuo (Versao 1): %e\n", norma1);
    printf("Tempo gasto (Versao 1): %lf ms\n", tTotal1);

    // VERSÃO 2: Calculando diagonais e termos independentes diretamente
    printf("\n--- VERSÃO 2: Gauss-Seidel calculando valores diretamente ---\n");
    real_t *Y2 = (real_t *) calloc(edo.n, sizeof(real_t)); // Resultado da EDO.
    
    tTotal2 = gaussSeidel_EDO(&edo, Y2, MAXIT);

    printf("\nSolucao (Y2):\n");
    prnVetor(Y2, edo.n);

    real_t norma2 = normaL2_residuo_EDO(&edo, Y2);
    printf("\nNorma L2 do Residuo (Versao 2): %e\n", sqrt(norma2));
    printf("Tempo gasto (Versao 2): %lf ms\n", tTotal2);

    // COMPARAÇÃO DAS VERSÕES
    printf("\n--- COMPARAÇÃO ---\n");
    printf("Razao de tempo (Versao 2 / Versao 1): %.3f\n", tTotal2/tTotal1);
    
    // Calcular diferença entre as soluções
    real_t diff_max = 0.0, diff_norm = 0.0;
    for (int i = 0; i < n; i++) {
        real_t diff = fabs(Y1[i] - Y2[i]);
        if (diff > diff_max) diff_max = diff;
        diff_norm += diff * diff;
    }
    diff_norm = sqrt(diff_norm);
    
    printf("Diferenca maxima entre solucoes: %e\n", diff_max);
    printf("Norma L2 da diferenca entre solucoes: %e\n", diff_norm);
    
    // Liberação de memória
    free(Y1);
    free(Y2);
    free(SL->D);
    free(SL->Di);
    free(SL->Ds);
    free(SL->B);
    free(SL);
  }
  
  return 0;
}

real_t pp (real_t x)
{
  return x+1;
}

real_t qq (real_t x)
{
  return -2*x;
}

real_t rr (real_t x)
{
  // Escolher apenas um dos retornos abaixo. Qual o melhor?
  return (1-x*x)*exp(-x);
  //return (1+x)*(1-x)*exp(-x);
}


