#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

#define DIFF 0.0

#define NRAND    ((real_t) random() / RAND_MAX)  // drand48() 
#define SRAND(a) srandom(a) // srand48(a)

real_t styblinski_tang_n(real_t *x, int n) 
{
  real_t soma = 0.0;
  for (int i = 0; i < n; i++) {
    soma += x[i]*x[i]*x[i]*x[i] - 16*x[i]*x[i] + 5*x[i];
  }
  return 0.5 * soma;
}

// Integral Monte Carlo da função Styblinski-Tang de 2 variáveis
real_t styblinskiTang(real_t a, real_t b, int n, int namostras)
{
  real_t soma = 0.0;
  real_t *x = malloc(n * sizeof(real_t));
  
  printf("Metodo de Monte Carlo (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), variaveis = 2\n", a, b, namostras);
  
  rtime_t t_inicial = timestamp();
  
  for (int k = 0; k < namostras; k++) {
    for (int i = 0; i < n; i++) {
      x[i] = a + (b - a) * NRAND;
    }
    soma += styblinski_tang_n(x, n);
  }

  real_t volume = pow(b - a, n);
  real_t resultado = volume * soma / namostras;
  
  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}

// Método dos retângulos (2D)
real_t retangulos_xy(real_t a, real_t b, int npontos) 
{
  real_t h = (b - a) / npontos;
  real_t soma = 0.0;
  
  printf("Metodo dos Retangulos (x, y).\n");
  printf("a = (%f), b = (%f), n = (%d), h = (%lg)\n", a, b, npontos, h);
  
  rtime_t t_inicial = timestamp();
  
  for (int i = 0; i < npontos; i++) {
    real_t x = a + (i + 0.5) * h;
    for (int j = 0; j < npontos; j++) {
      real_t y = a + (j + 0.5) * h;
      real_t ponto[2] = {x, y};
      soma += styblinski_tang_n(ponto, 2);
    }
  }

  real_t resultado = soma * h * h;
  
  rtime_t t_final = timestamp();
  printf("Tempo decorrido: %f seg.\n", t_final - t_inicial);
  
  return resultado;
}

int main(int argc, char **argv) {

  if (argc < 5) {
    printf("Utilização: %s inicial final n_amostras n_variaveis\n", argv[0]);
    return 1;
  }

  real_t a = atof(argv[1]);
  real_t b = atof(argv[2]);
  int n_amostras = atoi(argv[3]);
  int n_variaveis = atoi(argv[4]);

  SRAND(20252);

  printf("========================================\n");
  printf("Integracao da funcao Styblinski-Tang\n");
  printf("========================================\n");

  if (n_variaveis == 2) {
    real_t resultado_ret = retangulos_xy(a, b, (int) sqrt(n_amostras));
    real_t resultado_mc = styblinskiTang(a, b, 2, n_amostras);

    printf("\nComparação (2D)\n");
    printf("Retângulos = %.10lf\n", resultado_ret);
    printf("Monte Carlo = %.10lf\n", resultado_mc);
  }
  else if (n_variaveis == 4 || n_variaveis == 8) {
    printf("\nMetodo de Monte Carlo para %d variaveis.\n", n_variaveis);
    real_t resultado_mc = styblinskiTang(a, b, n_variaveis, n_amostras);
    printf("Resultado (Monte Carlo %dD): %.6f\n", n_variaveis, resultado_mc);
  } else {
    printf("Número de variáveis inválido. Use 2, 4 ou 8.\n");
    return 1;
  }
  
  return 0;
}