#include <stdio.h>
#include <string.h>
#include <math.h>

#include "utils.h"

/*  Retorna tempo em milisegundos desde EPOCH

    Forma de uso:
 
    rtime_t tempo;
    tempo = timestamp();
    <trecho de programa do qual se deseja medir tempo>
    tempo = timestamp() - tempo;
*/

rtime_t timestamp (void)
{
  struct timespec tp;
  clock_gettime(CLOCK_MONOTONIC_RAW, &tp);
  return ( (rtime_t) tp.tv_sec*1.0e3 + (rtime_t) tp.tv_nsec*1.0e-6 );
}

/* Gera string '<baseName>_n'
 * Por exemplo, se baseName = "ABC" e n = 10,
 *  Função retorna a string "ABC_10"
 * Útil para gerar marcadores para LIKWID
 */
string_t markerName(string_t baseName, int n)
{
    string_t mark = (string_t) malloc( (strlen(baseName)+1) + numDigits(n) + 1 );

  sprintf(mark, "%s_%u", baseName,n);

  // printf("*** %s\n", mark);

  return mark;

}

void solveTridiag(Tridiag *sl, real_t *sol)
{
    int n = sl->n;
    real_t *c_prime = (real_t *) malloc(n * sizeof(real_t));
    real_t *d_prime = (real_t *) malloc(n * sizeof(real_t));

    // Forward elimination
    c_prime[0] = sl->Ds[0] / sl->D[0];
    d_prime[0] = sl->B[0] / sl->D[0];

    for (int i = 1; i < n; i++) {
        real_t m = 1.0 / (sl->D[i] - sl->Di[i] * c_prime[i-1]);
        c_prime[i] = sl->Ds[i] * m;
        d_prime[i] = (sl->B[i] - sl->Di[i] * d_prime[i-1]) * m;
    }

    // Backward substitution
    sol[n-1] = d_prime[n-1];
    for (int i = n - 2; i >= 0; i--) {
        sol[i] = d_prime[i] - c_prime[i] * sol[i+1];
    }

    free(c_prime);
    free(d_prime);
}

