#include <stdio.h>
#include <stdlib.h> /* exit, malloc, calloc, etc. */
#include <string.h>
#include <getopt.h> /* getopt */
#include <time.h>
#include <likwid.h>

#include "matriz.h"
#include "utils.h"

/**
 * Exibe mensagem de erro indicando forma de uso do programa e termina
 * o programa.
 */

static void usage(char *progname)
{
  fprintf(stderr, "Forma de uso: %s [ <ordem> ] \n", progname);
  exit(1);
}

/**
 * Programa principal
 * Forma de uso: matmult [ -n <ordem> ]
 * -n <ordem>: ordem da matriz quadrada e dos vetores
 *
 */

int main(int argc, char *argv[])
{
  int n = DEF_SIZE;

  MatRow mRow_1, mRow_2, resMat, resMat_otim;
  Vetor vet, res, res_otim;
  rtime_t tempo;

  /* =============== TRATAMENTO DE LINHA DE COMANDO =============== */

  if (argc < 2)
    usage(argv[0]);

  n = atoi(argv[1]);

  /* ================ FIM DO TRATAMENTO DE LINHA DE COMANDO ========= */

  LIKWID_MARKER_INIT;

  srandom(20232);

  res = geraVetor(n, 0);
  res_otim = geraVetor(n, 0);
  resMat = geraMatRow(n, n, 1);
  resMat_otim = geraMatRow(n, n, 1);

  mRow_1 = geraMatRow(n, n, 0);
  mRow_2 = geraMatRow(n, n, 0);

  vet = geraVetor(n, 0);

  if (!res || !res_otim || !resMat || !resMat_otim || !mRow_1 || !mRow_2 || !vet)
  {
    fprintf(stderr, "Falha em alocação de memória !!\n");
    liberaVetor((void *)mRow_1);
    liberaVetor((void *)mRow_2);
    liberaVetor((void *)resMat);
    liberaVetor((void *)resMat_otim);
    liberaVetor((void *)vet);
    liberaVetor((void *)res);
    liberaVetor((void *)res_otim);
    exit(2);
  }

#ifdef _DEBUG_
  prnMat(mRow_1, n, n);
  prnMat(mRow_2, n, n);
  prnVetor(vet, n);
  printf("=================================\n\n");
#endif /* _DEBUG_ */

  // Multiplicação Matriz-Vetor
  LIKWID_MARKER_START("matVet");
  tempo = timestamp();
  multMatVet(mRow_1, vet, n, n, res);
  tempo = timestamp() - tempo;
  LIKWID_MARKER_STOP("matVet");
  printf("%d,%.10lg,", n, tempo);

  // Multiplicação Matriz-Vetor Otimizada
  LIKWID_MARKER_START("matVet_otim");
  tempo = timestamp();
  multMatVet_otim(mRow_1, vet, n, n, res_otim);
  tempo = timestamp() - tempo;
  LIKWID_MARKER_STOP("matVet_otim");
  printf("%.10lg,", tempo);

  // Multiplicação Matriz-Matriz
  LIKWID_MARKER_START("matMat");
  tempo = timestamp();
  multMatMat(mRow_1, mRow_2, n, resMat);
  tempo = timestamp() - tempo;
  LIKWID_MARKER_STOP("matMat");
  printf("%.10lg,", tempo);

  // Multiplicação Matriz-Matriz Otimizada
  LIKWID_MARKER_START("matMat_otim");
  tempo = timestamp();
  multMatMat_otim(mRow_1, mRow_2, n, resMat_otim);
  tempo = timestamp() - tempo;
  LIKWID_MARKER_STOP("matMat_otim");
  printf("%.10lg\n", tempo);

#ifdef _DEBUG_
  prnVetor(res, n);
  prnMat(resMat, n, n);
#endif /* _DEBUG_ */

  liberaVetor((void *)mRow_1);
  liberaVetor((void *)mRow_2);
  liberaVetor((void *)resMat);
  liberaVetor((void *)resMat_otim);
  liberaVetor((void *)vet);
  liberaVetor((void *)res);
  liberaVetor((void *)res_otim);

  LIKWID_MARKER_CLOSE;

  return 0;
}
