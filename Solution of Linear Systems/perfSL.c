#include <stdio.h>
#include <fenv.h>
#include <stdlib.h>
#include "utils.h"
#include "sislin.h"
#include "eliminacaoGauss.h"
#include "gaussSeidel.h"

int main() {
    // Define o modo de arredondamento para baixo para todos os cálculos de ponto flutuante
    fesetround(FE_DOWNWARD);

    // Lê o sistema linear da entrada padrão
    SistLinear_t *s_orig = lerSisLin();
    if (!s_orig) {
        fprintf(stderr, "Erro ao ler o sistema linear.\n");
        return -1;
    }

    // Aloca vetores para as soluções e resíduos
    real_t *x_eg = (real_t *) calloc(s_orig->n, sizeof(real_t));
    real_t *r_eg = (real_t *) calloc(s_orig->n, sizeof(real_t));
    real_t *y_gs = (real_t *) calloc(s_orig->n, sizeof(real_t));
    real_t *r_gs = (real_t *) calloc(s_orig->n, sizeof(real_t));

    if (!x_eg || !r_eg || !y_gs || !r_gs) {
        fprintf(stderr, "Erro de alocação de memória.\n");
        liberaSisLin(s_orig);
        free(x_eg);
        free(r_eg);
        free(y_gs);
        free(r_gs);
        return -1;
    }

    // --- Eliminação de Gauss ---
    SistLinear_t *s_eg = dupSisLin(s_orig);
    if (s_eg) {
        rtime_t tempo_eg = timestamp();
        retrosubst(s_eg, x_eg);
        tempo_eg = timestamp() - tempo_eg;

        residuo(s_orig, x_eg, r_eg, s_orig->n);

        printf("EG:\n");
        printf("%.8f ms\n", tempo_eg);
        prnVetor(x_eg, s_orig->n);
        prnVetor(r_eg, s_orig->n);
        printf("\n");

        liberaSisLin(s_eg);
    }

    // --- Gauss-Seidel ---
    SistLinear_t *s_gs = dupSisLin(s_orig);
    if (s_gs) {
        real_t erro = TOL;
        int maxit = MAXIT;
        real_t norma;
        
        rtime_t tempo_gs = timestamp();
        int it = gaussSeidel(s_gs, y_gs, erro, maxit, &norma);
        tempo_gs = timestamp() - tempo_gs;

        residuo(s_orig, y_gs, r_gs, s_orig->n);

        printf("GS [ %d iterações ]:\n", it);
        printf("%.8f ms\n", tempo_gs);
        prnVetor(y_gs, s_orig->n);
        prnVetor(r_gs, s_orig->n);

        liberaSisLin(s_gs);
    }

    // Libera toda a memória alocada
    liberaSisLin(s_orig);
    free(x_eg);
    free(r_eg);
    free(y_gs);
    free(r_gs);

    // Restaura o modo de arredondamento padrão
    fesetround(FE_TONEAREST);

    return 0;
}