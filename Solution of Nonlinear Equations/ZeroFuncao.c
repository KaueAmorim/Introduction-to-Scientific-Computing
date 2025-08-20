#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, bool ehRapido, int *it, real_t *raiz) {

}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, bool ehRapido, int *it, real_t *raiz) {
    real_t xm_old, xm_new;

    xm_new = (a + b) / 2;

    real_t f1=0.0, df1=0.0, f2=0.0, df2=0.0;
    calcPolinomio_rapido(p, a, &f1, &df1);
    calcPolinomio_rapido(p, xm_new, &f2, &df2);

    if (f1 * f2 < 0) {
        b = xm_new;
    }
    else if (f1 * f2 > 0) {
        a = xm_new;
    }
    else {
        return xm_new;
    }

    do {
        xm_old = xm_new;
        xm_new = (a + b) / 2;

        calcPolinomio_rapido(p, a, &f1, &df1);
        calcPolinomio_rapido(p, xm_new, &f2, &df2);

        if (f1 * f2 < 0) {
            b = xm_new;
        }
        else if (f1 * f2 > 0) {
            a = xm_new;
        }
        else {
            return xm_new;
        }
    } while (fabs(xm_new - xm_old) > criterioParada);

    *raiz = xm_new;
    return 0.0;
}


void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    *px = 0;
    *dpx = 0;
    for (int i = p.grau; i > 0; --i) {
        *px = *px * x + p.p[i];
        *dpx = *dpx * x + *px;
    }
    *px = *px * x + p.p[0];
}


void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    *px = 0;
    *dpx = 0;
    for (int i = p.grau; i > 0; --i) {
        *px += p.p[i] * pow(x, i);
        *dpx += i * p.p[i] * pow(x, i - 1);
    }
    *px += p.p[0] * pow(x, 0);
}