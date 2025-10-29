#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

void imprimirSistemaLinear(real_t ***A, real_t *B, int n, int k);

void criaKDiagonal(int n, int k, real_t ***A, real_t **B);

void genSimetricaPositiva(real_t ***A, real_t **b, int n, int k, int *new_k, real_t ***ASP, real_t **bsp, rtime_t *tempo);
void geraDLU (real_t ***A, int n, int k, real_t **D, real_t ***L, real_t ***U, rtime_t *tempo);
void geraPreCond(real_t **D, real_t ***L, real_t ***U, real_t w, int n, int k, real_t ***M, rtime_t *tempo);
void aplicaPreCondicionador(real_t ***M, real_t *r, real_t *v, int n, int k, real_t w);
real_t calcResiduoSL (real_t ***A, real_t **b, real_t **X, int n, int k, rtime_t *tempo);
int gradientesConjugados(real_t ***A, real_t *b, real_t *x, real_t ***M, int n, int k, 
                        real_t epsilon, int maxit, real_t *norma, rtime_t *tempo_iter, real_t w);

#endif // __SISLIN_H__

