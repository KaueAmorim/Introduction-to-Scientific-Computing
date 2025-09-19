/*
 * Nome: Kaue (substituir pelo seu nome)
 * GRR: XXXXXXXXX (substituir pelo seu GRR)
 * 
 * Implementação das funções para resolução de Equações Diferenciais Ordinárias
 * usando diferenças finitas e sistemas lineares tridiagonais.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "edo.h"

#define MAXIT 100
#define EPS 1.0e-5

/**
 * @brief Gera um sistema linear tridiagonal a partir de uma EDO
 * 
 * Discretiza a EDO y'' + py' + qy = r(x) usando diferenças finitas centradas
 * em uma malha uniforme, gerando um sistema Ax = b onde A é tridiagonal.
 * 
 * A discretização usa:
 * - y''(xi) ≈ (y[i-1] - 2*y[i] + y[i+1])/h²
 * - y'(xi) ≈ (y[i+1] - y[i-1])/(2*h)
 * 
 * @param edo Ponteiro para a estrutura EDO com os parâmetros da equação
 * @return Ponteiro para o sistema tridiagonal, ou NULL se erro de alocação
 */

Tridiag *genTridiag (EDo *edo)
{
  Tridiag *sl;
  real_t x, rx;
  int n = edo->n;
  
  sl = (Tridiag *) malloc (sizeof(Tridiag));
  sl->n = edo->n;

  sl->D = (real_t *) calloc(n, sizeof(real_t));
  sl->Di = (real_t *) calloc(n, sizeof(real_t));
  sl->Ds = (real_t *) calloc(n, sizeof(real_t));
  sl->B = (real_t *) calloc(n, sizeof(real_t));

  real_t h = (edo->b - edo->a)/(n+1);

  for (int i=0; i < n; ++i) {
    x = edo->a + (i+1)*h;
    rx = edo->r1*x + edo->r2*x*x + edo->r3*cos(x) + edo->r4*exp(x);
    
    sl->B[i] = h*h * rx;
    sl->Di[i] = 1 - h * edo->p/2.0;
    sl->D[i] = -2 + h*h * edo->q;
    sl->Ds[i] = 1 + h * edo->p/2.0;
  }

  sl->B[0] -= edo->ya * (1 - h*edo->p/2.0);
  sl->B[n-1] -= edo->yb * (1 + h*edo->p/2.0);
  
  return sl;
}

/**
 * @brief Resolve um sistema tridiagonal usando o método iterativo de Gauss-Seidel
 * 
 * O método de Gauss-Seidel é um método iterativo que atualiza cada componente
 * da solução usando os valores mais recentes das outras componentes:
 * x[i]^(k+1) = (b[i] - sum(a[i,j]*x[j]^(k+1), j<i) - sum(a[i,j]*x[j]^(k), j>i)) / a[i,i]
 * 
 * @param sl Ponteiro para o sistema tridiagonal
 * @param sol Vetor onde será armazenada a solução (deve estar alocado com tamanho n)
 * @param norma_residuo Ponteiro onde será armazenada a norma L2 do resíduo final
 * @return Número de iterações realizadas até convergência ou limite máximo
 */
int gaussSeidel(Tridiag *sl, real_t *sol, real_t *norma_residuo) {
    int n = sl->n;
    real_t *residuo = (real_t *)malloc(n * sizeof(real_t));
    int iter;

    // Chute inicial: vetor nulo
    for (int i = 0; i < n; ++i) {
        sol[i] = 0.0;
    }
    for (iter = 0; iter < MAXIT; ++iter) {
        // Atualização Gauss-Seidel: usa valores já atualizados na mesma iteração
        for (int i = 0; i < n; ++i) {
            real_t soma = 0.0;
            if (i > 0) {
                soma += sl->Di[i] * sol[i - 1];  // valor já atualizado
            }
            if (i < n - 1) {
                soma += sl->Ds[i] * sol[i + 1];  // valor da iteração anterior
            }
            sol[i] = (sl->B[i] - soma) / sl->D[i];
        }

        // Calcula resíduo r = b - Ax para verificar convergência
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

        // Verifica critério de parada por tolerância
        if (*norma_residuo <= EPS) {
            break;
        }
    }

    free(residuo);
    return (iter < MAXIT && *norma_residuo <= EPS) ? iter + 1 : iter;
}

/**
 * @brief Imprime a solução do sistema em formato horizontal
 * 
 * Formata a saída da solução em uma única linha com valores separados por espaços,
 * usando o formato de precisão definido em FORMAT (15 casas decimais).
 * 
 * @param sol Vetor com a solução do sistema
 * @param n Tamanho do vetor solução
 */
void prnSolucao(real_t *sol, int n) {
    printf("\n");
    for (int i = 0; i < n; ++i) {
        printf(FORMAT, sol[i]);
    }
    printf("\n");
}


// Exibe SL na saída padrão
void prnEDOsl (EDo *edoeq)
{
  int n = edoeq->n, i, j;
  real_t x, b, d, di, ds,rx;
  real_t h = (edoeq->b - edoeq->a)/(n+1);

  printf ("%d\n", n);

  for (i=0; i < n; ++i) {
    x = edoeq->a + (i+1)*h;
    rx = edoeq->r1*x + edoeq->r2*x*x + edoeq->r3*cos(x) + edoeq->r4*exp(x);
    
    b = h*h * rx; 
    di = 1 - h * edoeq->p/2.0;
    d = -2 + h*h * edoeq->q;
    ds = 1 + h * edoeq->p/2.0;
      
    for (j=0; j < n; ++j) {
      if (i == j)
	printf (FORMAT,d);
      else if (j == i-1 && i)
	printf (FORMAT,di);
      else if (j == i+1 && i != n-1)
	printf (FORMAT,ds);
      else
	printf(FORMAT, 0.0);
    }
      
    if (i == 0)
      b -= edoeq->ya * (1 - h*edoeq->p/2.0);
    else if (i == n-1)
      b -= edoeq->yb * (1 + h*edoeq->p/2.0);

    printf (FORMAT, b);
      
    printf ("\n");
  }
}

/**
 * @brief Resolve sistema tridiagonal usando algoritmo de Thomas (eliminação direta)
 * 
 * Implementa o algoritmo de Thomas para resolução eficiente de sistemas
 * tridiagonais Ax = b em O(n) operações. O método consiste em duas fases:
 * 1. Eliminação para frente (forward elimination) 
 * 2. Substituição para trás (backward substitution)
 * 
 * Este método é específico para sistemas tridiagonais e muito mais eficiente
 * que a eliminação gaussiana completa O(n³).
 * 
 * @param sl Ponteiro para sistema tridiagonal (Di, D, Ds, B)
 * @param sol Vetor onde será armazenada a solução (deve estar alocado)
 * 
 * NOTA: Assume que o sistema é bem condicionado (sem elementos nulos na diagonal)
 */
void solveTridiag(Tridiag *sl, real_t *sol)
{
    int n = sl->n;
    real_t *c_prime = (real_t *) malloc(n * sizeof(real_t));
    real_t *d_prime = (real_t *) malloc(n * sizeof(real_t));

    // Fase 1: Eliminação para frente
    // Elimina a diagonal inferior, modificando diagonal principal e vetor b
    c_prime[0] = sl->Ds[0] / sl->D[0];
    d_prime[0] = sl->B[0] / sl->D[0];

    for (int i = 1; i < n; i++) {
        real_t m = 1.0 / (sl->D[i] - sl->Di[i] * c_prime[i-1]);
        c_prime[i] = sl->Ds[i] * m;
        d_prime[i] = (sl->B[i] - sl->Di[i] * d_prime[i-1]) * m;
    }

    // Fase 2: Substituição para trás  
    // Resolve o sistema triangular superior resultante
    sol[n-1] = d_prime[n-1];
    for (int i = n - 2; i >= 0; i--) {
        sol[i] = d_prime[i] - c_prime[i] * sol[i+1];
    }

    free(c_prime);
    free(d_prime);
}


