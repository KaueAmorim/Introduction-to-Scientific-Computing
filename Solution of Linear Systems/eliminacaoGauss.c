/* Matriz  'normal' (vetor  de  ponteiros (linhas  matriz) para  vetores
   (colunas da matriz), estilo 'Mazieiro/Prog 2'
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "eliminacaoGauss.h"

/* Encontra a linha do SL 'C' que tem o maior valor na coluna 'k' da matriz
   de coeficientes.
   RETORNO:  índice da  linha que  tem o  maior valor  na coluna  'k' da
             matriz de coeficientes
*/
static int encontraMax(SistLinear_t *C, int k)
{
   int linha = k;
   real_t max = fabs(C->A[k][k]);

   for(int i = k+1; i < C->n; i++){
      if(fabs(C->A[i][k]) > max){
         max = fabs(C->A[i][k]);
         linha = i;
      }
   }

   return linha;
}

/* Troca de lugar entre si as linhas 'k' e 'p' do SL 'C' */
static void trocaLinha (SistLinear_t *C, int k, int p)
{
   if(k == p){
      return;
   }

   real_t temp;

   for(int i = 0; i < C->n; i++){
      temp = C->A[k][i];
      C->A[k][i] = C->A[p][i];
      C->A[p][i] = temp;
   }

   temp = C->b[k];
   C->b[k] = C->b[p];
   C->b[p] = temp;
}


/* Seja um S.L. de ordem 'n'
   C = A|B em Ax=B
 */
void triangulariza( SistLinear_t *C )
{
   for(int i = 0; i < C->n - 1; i++) {

      int iPivo = encontraMax(C, i);
      if(iPivo != i){
         trocaLinha(C, i, iPivo);
      }

      for(int k = i+1; k < C->n; k++) {
         real_t m = C->A[k][i] / C->A[i][i];
         for(int j = i+1; j < C->n; j++) {
            C->A[k][j] -= C->A[i][j] * m;
         }
         C->b[k] -= C->b[i] * m;
      }
   }
}

void retrosubst( SistLinear_t *C, real_t *X )
{
   triangulariza(C);

   for(int i = (C->n - 1); i >= 0; i--) {
      X[i] = C->b[i];
      for(int j = i+1; j < C->n; j++) {
         X[i] -= C->A[i][j] * X[j];
      }
      X[i] /= C->A[i][i];
   }
}