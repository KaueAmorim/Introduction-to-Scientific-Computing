#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fenv.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main () {

  // Configurar arredondamento para baixo para todo o programa
  fesetround(FE_DOWNWARD);

  real_t a, b;
  Polinomio pol;
  int it = 0; 
  real_t raiz = 0.0, erro = 0.0, tempo = 0.0;

  scanf("%d", &pol.grau);

  pol.p = (real_t*) malloc(sizeof(real_t) * (pol.grau + 1));

  for (int i=pol.grau; i >=0; --i)
    scanf("%lf", &pol.p[i]);

  scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.

  printf("RAPIDO\n\n");

  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = bisseccao(pol, a, b, i, &it, &raiz, calcPolinomio_rapido);
    tempo = timestamp() - tempo;

    printf("bissec  %.15e %.15e %4d  %.8e\n", raiz, erro, it, tempo);
  }
  
  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = newtonRaphson(pol, (a+b)/2, i, &it, &raiz, calcPolinomio_rapido);
    tempo = timestamp() - tempo;

    printf("newton  %.15e %.15e %4d  %.8e\n", raiz, erro, it, tempo);
  }

  printf("\nLENTO\n\n");

  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = bisseccao(pol, a, b, i, &it, &raiz, calcPolinomio_lento);
    tempo = timestamp() - tempo;

    printf("bissec  %.15e %.15e %4d  %.8e\n", raiz, erro, it, tempo);
  }

  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = newtonRaphson(pol, (a+b)/2, i, &it, &raiz, calcPolinomio_lento);
    tempo = timestamp() - tempo;

    printf("newton  %.15e %.15e %4d  %.8e\n", raiz, erro, it, tempo);
  }

  free(pol.p);
  fesetround(FE_TONEAREST);

  return 0;
}

