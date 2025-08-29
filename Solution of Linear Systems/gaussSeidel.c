#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"
#include "gaussSeidel.h"

/**
 * Executa o método de Gauss-Seidel para resolver um sistema linear.
 *
 * @param C Ponteiro para a estrutura do sistema linear.
 * @param X Vetor de incógnitas (deve ser inicializado com zeros pelo chamador).
 * @param erro Critério de parada (tolerância do erro).
 * @param maxit Número máximo de iterações.
 * @param norma Ponteiro para armazenar a norma do erro da última iteração.
 * @return O número de iterações executadas em caso de convergência, ou -1 caso contrário.
 */
int gaussSeidel(SistLinear_t *C, real_t *X, real_t erro, int maxit, real_t *norma)
{
    int i, j, k;
    real_t soma;

    // Loop principal das iterações
    for (k = 0; k < maxit; k++) {
        // Inicializa a norma (erro máximo) para a iteração atual
        *norma = 0.0;

        // Loop para cada equação (linha) do sistema
        for (i = 0; i < C->n; i++) {
            // Guarda o valor de X[i] da iteração anterior para calcular o erro
            real_t x_old = X[i];

            soma = C->b[i]; // Inicia a soma com o termo independente b_i
            
            // Calcula o somatório da fórmula de Gauss-Seidel
            for (j = 0; j < C->n; j++) {
                if (i != j) {
                    // Utiliza os valores de X já atualizados nesta mesma iteração (k)
                    // para j < i, que é a característica principal do método.
                    soma -= C->A[i][j] * X[j];
                }
            }

            // Calcula o novo valor de X[i]
            X[i] = soma / C->A[i][i];

            // Calcula o erro absoluto para a incógnita atual
            real_t diff = fabs(X[i] - x_old);
            
            // Atualiza a norma se o erro atual for o maior encontrado nesta iteração
            if (diff > *norma) {
                *norma = diff;
            }
        }

        // Verifica o critério de parada
        if (*norma < erro) {
            return k + 1; // Sucesso! Retorna o número de iterações (k começa em 0)
        }
    }

    // Se o loop terminar, o método não convergiu no número máximo de iterações
    return -1; 
}

