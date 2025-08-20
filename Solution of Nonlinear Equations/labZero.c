#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{

  real_t a, b;
  Polinomio pol;
  int it; 
  real_t raiz, erro, tempo;

  scanf("%d", &pol.grau);

  pol.p = (real_t*) malloc(sizeof(real_t) * (pol.grau + 1));

  for (int i=pol.grau; i >=0; --i)
    scanf("%lf", &pol.p[i]);

  scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.

  printf("RAPIDO\n");

  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = bisseccao(pol, a, b, i, true, &it, &raiz);
    tempo = timestamp() - tempo;

    printf("bissec %.15e %.15e %d %.8e", raiz, erro, it, tempo);
  }
  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = newtonRaphson(pol, a, i, true, &it, &raiz);
    tempo = timestamp() - tempo;

    printf("newton %.15e %.15e %d %.8e", raiz, erro, it, tempo);
  }

  printf("LENTO\n");

  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = bisseccao(pol, a, b, i, false, &it, &raiz);
    tempo = timestamp() - tempo;

    printf("bissec %.15e %.15e %d %.8e", raiz, erro, it, tempo);
  }
  for(int i = 1; i <= 3; i++){
    tempo = timestamp();
    erro = newtonRaphson(pol, a, i, false, &it, &raiz);
    tempo = timestamp() - tempo;

    printf("newton %.15e %.15e %d %.8e", raiz, erro, it, tempo);
  }


  // Restante do programa a partir daqui

  return 0;
}

