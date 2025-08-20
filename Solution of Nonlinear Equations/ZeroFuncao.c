#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"
#include "DoubleType.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, void (*f)(Polinomio, real_t, real_t *, real_t *)){
    
    real_t xm_new = x0, xm_old, erro = 0.0;
    real_t fx=0.0, dfx=0.0;
    int parar = 0;
    *it = 0;

    do {
        xm_old = xm_new;
        f(p, xm_new, &fx, &dfx);
        if (dfx == 0.0) {
            // Derivada zero, não pode continuar
            break;
        }

        xm_new = xm_old - (fx / dfx);

        (*it)++;

        switch (criterioParada) {
            case 1:
                erro = fabs((xm_new - xm_old)/xm_new);
                parar = erro <= EPS;
                break;
            case 2:
                erro = fabs(fx);
                parar = erro <= ZERO;
                break;
            case 3:
                Double_t d1, d2;
                d1.f = xm_new;
                d2.f = xm_old;

                // Corrige para valores negativos
                int64_t ulps_diff = llabs((int64_t)d1.i - (int64_t)d2.i);
                if(ulps_diff > 0) {
                    ulps_diff--;
                }

                erro = (real_t)ulps_diff;
                parar = ulps_diff <= ULPS;
                break;
            default:
                return 0.0; // Critério de parada inválido
        }
    } while (!parar && *it < MAXIT);

    *raiz = xm_new;

    return erro;
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, void (*f)(Polinomio, real_t, real_t *, real_t *)) {

    real_t xm_old, xm_new, erro;
    real_t f1=0.0, df1=0.0, f2=0.0, df2=0.0;
    int parar = 0;
    *it = 2;

    xm_new = (a + b) / 2;
    f(p, a, &f1, &df1);
    f(p, xm_new, &f2, &df2);

    if (f1 * f2 < 0)
        b = xm_new;
    else if (f1 * f2 > 0)
        a = xm_new;
    else
        return xm_new;

    do {
        xm_old = xm_new;
        xm_new = (a + b) / 2;

        (*it)++;

        f(p, a, &f1, &df1);
        f(p, xm_new, &f2, &df2);

        if (f1 * f2 < 0)
            b = xm_new;
        else if (f1 * f2 > 0)
            a = xm_new;
        else
            return xm_new;

        switch (criterioParada) {
            case 1:
                erro = fabs((xm_new - xm_old)/xm_new);
                parar = erro <= EPS;
                break;
            case 2:
                erro = fabs(f2);
                parar = erro <= ZERO;
                break;
            case 3:
                Double_t d1, d2;
                d1.f = xm_new;
                d2.f = xm_old;

                // Corrige para valores negativos
                int64_t ulps_diff = llabs((int64_t)d1.i - (int64_t)d2.i);
                if(ulps_diff > 0) {
                    ulps_diff--;
                }

                erro = (real_t)ulps_diff;
                parar = ulps_diff <= ULPS;
                break;
            default:
                return 0.0; // Critério de parada inválido
        }
    } while (!parar && *it <= MAXIT);

    *raiz = xm_new;

    return erro;
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