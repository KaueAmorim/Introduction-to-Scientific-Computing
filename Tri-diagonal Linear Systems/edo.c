/*
 * Nome: Kauê de Amorim Silva
 * GRR: 20244719
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
#define NORMA_STOP EPS

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

int gaussSeidel_3Diag(Tridiag *sl, real_t *Y, int maxiter, real_t *norma) {
    const int n = sl->n;
    int it = 0;

    do {
        // Primeira linha não possui nenhum elemento da diagonal inferior
        Y[0] = (sl->B[0] - sl->Ds[0] * Y[1]) / sl->D[0];

        for (int i = 1; i < n - 1; ++i)
            Y[i] = (sl->B[i] - sl->Di[i - 1] * Y[i - 1] - sl->Ds[i] * Y[i + 1]) / sl->D[i];

        // Última linha não possui nenhum elemento da diagonal superior
        Y[n - 1] = (sl->B[n - 1] - sl->Di[n - 2] * Y[n - 2] ) / sl->D[n - 1];

        *norma = normaL2_3Diag(sl, Y);
        it++;
    } while (*norma > NORMA_STOP && it < maxiter);

    return it;
}

real_t normaL2_3Diag (Tridiag *sl, real_t *Y) {
    int n = sl->n;
    real_t normaL2 = 0.0;
    real_t residuo;

    // Primeiro elemento
    residuo = sl->B[0] - (sl->D[0] * Y[0] + sl->Ds[0] * Y[1]);
    normaL2 += residuo * residuo;

    // Elementos do meio
    for (int i = 1; i < n - 1; ++i) {
        residuo = sl->B[i] - (sl->Di[i-1] * Y[i-1] + sl->D[i] * Y[i] + sl->Ds[i] * Y[i+1]);
        normaL2 += residuo * residuo;
    }

    // Último elemento
    residuo = sl->B[n-1] - (sl->Di[n-2] * Y[n-2] + sl->D[n-1] * Y[n-1]);
    normaL2 += residuo * residuo;

    return sqrt(normaL2);
}

void prnSolucao(real_t *sol, int n) {
    printf("\n");
    for (int i = 0; i < n; ++i) {
        printf(FORMAT, sol[i]);
    }
    printf("\n");
}

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