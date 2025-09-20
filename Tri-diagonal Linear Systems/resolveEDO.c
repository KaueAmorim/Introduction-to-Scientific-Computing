/*
 * Nome: Kauê de Amorim Silva
 * GRR: 20244719
 * 
 * Programa principal para resolução de Equações Diferenciais Ordinárias (EDO).
 * 
 * Resolve EDOs da forma: y'' + py' + qy = r(x)
 * onde r(x) = r1*x + r2*x² + r3*cos(x) + r4*exp(x)
 * 
 * Usa discretização por diferenças finitas e método iterativo de Gauss-Seidel
 * para resolver o sistema linear tridiagonal resultante.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include <likwid.h>
#include "utils.h"
#include "edo.h"

#define MAXIT 100

/**
 * @brief Função principal do programa
 * 
 * Lê dados da entrada padrão no formato:
 * - Linha 1: número de pontos da malha (n)
 * - Linha 2: intervalo [a, b]
 * - Linha 3: condições de contorno y(a), y(b)  
 * - Linha 4: coeficientes p, q da EDO
 * - Linhas seguintes: coeficientes r1, r2, r3, r4 de r(x) (uma EDO por linha)
 * 
 * Para cada EDO:
 * 1. Gera e imprime o sistema linear tridiagonal
 * 2. Resolve usando Gauss-Seidel
 * 3. Imprime solução e estatísticas (iterações, resíduo, tempo)
 * 
 * @return 0 se execução bem-sucedida
 */
int main ()
{
  LIKWID_MARKER_INIT;

  // Configura arredondamento para baixo conforme especificação
  fesetround(FE_DOWNWARD);
  
  EDo edo;
  Tridiag *sl;
  real_t *sol;
  real_t norma_residuo, tempo;
  int iter;
  int edo_count = 0;

  // Leitura dos parâmetros gerais da EDO
  scanf("%d", &edo.n);
  scanf("%lf %lf", &edo.a, &edo.b);
  scanf("%lf %lf", &edo.ya, &edo.yb);
  scanf("%lf %lf", &edo.p, &edo.q);

  // Loop para processar múltiplas EDOs (diferentes r(x))
  while (scanf("%lf %lf %lf %lf", &edo.r1, &edo.r2, &edo.r3, &edo.r4) == 4) {
    edo_count++;
    string_t marker_name = markerName("Gauss-Seidel", edo_count);
    
    // Gera o sistema linear tridiagonal correspondente à EDO
    sl = genTridiag(&edo);

    // Imprime a ordem e matriz aumentada do sistema
    prnEDOsl(&edo);

    // Aloca vetor para armazenar a solução
    sol = (real_t *) malloc(edo.n * sizeof(real_t));
    
    // Inicializa chute inicial (vetor nulo)
    for (int i = 0; i < edo.n; ++i) {
        sol[i] = 0.0;
    }

    // Resolve o sistema usando Gauss-Seidel e mede o tempo
    tempo = timestamp();
    LIKWID_MARKER_START(marker_name);
    iter = gaussSeidel_3Diag(sl, sol, 100, &norma_residuo);
    LIKWID_MARKER_STOP(marker_name);
    tempo = timestamp() - tempo;

    // Imprime resultados: solução, iterações, resíduo, tempo
    prnSolucao(sol, edo.n);
    printf("%d\n", iter);
    printf(FORMAT "\n", norma_residuo);
    printf("  %.8e\n", tempo);

    // Libera memória alocada para esta EDO
    free(sl->D);
    free(sl->Di);
    free(sl->Ds);
    free(sl->B);
    free(sl);
    free(sol);
    free(marker_name);  // Libera string do nome do marcador
  }

  // Restaura modo de arredondamento padrão
  fesetround(FE_TONEAREST);

  LIKWID_MARKER_CLOSE;
  
  return 0;
}