#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include "utils.h"
#include "edo.h"

#define MAXIT 100
#define EPS 1.0e-5

void prnSolucao(real_t *sol, int n) {
    printf("\nSolucao:\n");
    for (int i = 0; i < n; ++i) {
        printf(FORMAT "\n", sol[i]);
    }
}

int gaussSeidel(Tridiag *sl, real_t *sol, real_t *norma_residuo) {
    int n = sl->n;
    real_t *residuo = (real_t *)malloc(n * sizeof(real_t));
    int iter;

    // Chute inicial: vetor nulo
    for (int i = 0; i < n; ++i) {
        sol[i] = 0.0;
    }

    for (iter = 0; iter < MAXIT; ++iter) {
        // Calcula nova solução
        for (int i = 0; i < n; ++i) {
            real_t soma = 0.0;
            if (i > 0) {
                soma += sl->Di[i] * sol[i - 1];
            }
            if (i < n - 1) {
                soma += sl->Ds[i] * sol[i + 1];
            }
            sol[i] = (sl->B[i] - soma) / sl->D[i];
        }

        // Calcula resíduo (r = b - Ax)
        *norma_residuo = 0.0;
        for (int i = 0; i < n; ++i) {
            real_t ax = sl->D[i] * sol[i];
            if (i > 0) {
                ax += sl->Di[i] * sol[i - 1];
            }
            if (i < n - 1) {
                ax += sl->Ds[i] * sol[i + 1];
            }
            residuo[i] = sl->B[i] - ax;
            *norma_residuo += residuo[i] * residuo[i];
        }
        *norma_residuo = sqrt(*norma_residuo);

        if (*norma_residuo <= EPS) {
            free(residuo);
            return iter + 1;
        }
    }

    free(residuo);
    return iter;
}


int main ()
{
  
  LIKWID_MARKER_INIT;

  // Define o modo de arredondamento para baixo para todos os cálculos de ponto flutuante
  fesetround(FE_DOWNWARD);
  
  EDo edo;
  Tridiag *sl;
  real_t *sol;
  real_t norma_residuo, tempo;
  int iter;

  // Leitura dos parâmetros da EDO
  scanf("%d", &edo.n);
  scanf("%lf %lf", &edo.a, &edo.b);
  scanf("%lf %lf", &edo.ya, &edo.yb);
  scanf("%lf %lf", &edo.p, &edo.q);

  // Loop para ler múltiplos coeficientes r(x)
  while (scanf("%lf %lf %lf %lf", &edo.r1, &edo.r2, &edo.r3, &edo.r4) == 4) {
    // Gera o sistema linear tridiagonal
    sl = genTridiag(&edo);

    // Imprime a matriz aumentada
    prnEDOsl(&edo);

    // Aloca vetor para a solução
    sol = (real_t *) malloc(edo.n * sizeof(real_t));

    // Resolve o sistema linear com Gauss-Seidel
    tempo = timestamp();
    iter = gaussSeidel(sl, sol, &norma_residuo);
    tempo = timestamp() - tempo;

    // Imprime a solução e as estatísticas
    prnSolucao(sol, edo.n);
    printf("\nIteracoes: %d\n", iter);
    printf("Norma L2 do residuo: " FORMAT "\n", norma_residuo);
    printf("Tempo: %.8e\n\n", tempo);


    // Libera a memória alocada para o sistema e a solução
    free(sl->D);
    free(sl->Di);
    free(sl->Ds);
    free(sl->B);
    free(sl);
    free(sol);
  }

  // Restaura o modo de arredondamento padrão
  fesetround(FE_TONEAREST);

  LIKWID_MARKER_CLOSE;
  
  return 0;
}