#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "utils.h"
#include "edo.h"

#define MAXIT 100

int main ()
{
  rtime_t tTotal;
  EDo edo = {5, 0, 1, -1, 0, pp, qq, rr}; // definição da EDO do exercício
  real_t *Y = (real_t *) calloc(edo.n, sizeof(real_t)); // Resultado da EDO.

  Tridiag *SL = genTridiag(&edo);
  tTotal = gaussSeidel_3Diag(SL, Y, MAXIT);

  prnTriDiagonal(SL);
  prnVetor(Y, edo.n);

  gaussSeidel_EDO(&edo, Y, MAXIT);
  
  return 0;
}