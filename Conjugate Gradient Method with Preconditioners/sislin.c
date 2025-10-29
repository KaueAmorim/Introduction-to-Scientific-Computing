#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

/** Imprime sistema linear
 * @param ***A matriz de coeficientes k-diagonal
 * @param *B vetor de termos independentes
 * @param n ordem do sistema linear
 * @param k quantidade de diagonais
 */
void imprimirSistemaLinear(real_t ***A, real_t *B, int n, int k)
{
    int i, j;
    int m = k/2;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            real_t val = 0.0;
            int diag_offset = j - i;
            if (ABS(diag_offset) <= m) {
                int diag_idx = diag_offset + m;
                val = (*A)[diag_idx][i];
            }
            printf("%.16g ", val);
        }
        printf(" | %.16g\n", B[i]);
    }
    printf("\n");
}

static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k );
static inline real_t generateRandomB( unsigned int k );

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
    static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
    return ( (i==j) ? (real_t)(k<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline real_t generateRandomB( unsigned int k )
{
    static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
    return (real_t)(k<<2) * (real_t)random() * invRandMax;
}


/**
 * Cria matriz 'A' k-diagonal e vetor de termos independentes B 
 * @param n dimensao do sistema linear (n > 10)
 * @param k numero de diagonais da matriz A (k > 1 e impar)
 * @param **A matriz de coeficientes (k x n)
 * @param **B vetor de termos independentes
 */
void criaKDiagonal(int n, int k, real_t ***A, real_t **B)
{
    *A = malloc(sizeof(real_t*) * k);

    int i, j;
    int m = k/2;

    for (i = 0; i < k; i++) {
        (*A)[i] = malloc(sizeof(real_t) * n);

        if (i < m) { // diagonal inferior
            int diag_offset = m - i;
            for (j = 0; j < diag_offset; j++) 
                (*A)[i][j] = 0.0;
            for (j = diag_offset; j < n; j++) 
                (*A)[i][j] = generateRandomA(j, j - diag_offset, k);
        } else if (i > m) { // diagonal superior
            int diag_offset = i - m;
            for (j = 0; j < n - diag_offset; j++) 
                (*A)[i][j] = generateRandomA(j, j + diag_offset, k);
            for (j = n - diag_offset; j < n; j++) 
                (*A)[i][j] = 0.0;
        } else { // diagonal principal (i == m)
            for (j = 0; j < n; j++) 
                (*A)[i][j] = generateRandomA(j, j, k);
        }
    }

    *B = malloc(sizeof(real_t) * n);
    for (i = 0; i < n; i++) 
        (*B)[i] = generateRandomB(k);
}

/** Gera matriz simetrica positiva a partir de uma matriz k-diagonal
 * @param ***A matriz de coeficientes
 * @param **b vetor de termos independentes
 * @param n dimensão do sistema linear
 * @param k número de diagonais da matriz A
 * @param ***ASP matriz A simetrica positiva
 * @param **bsp vetor de termos independentes para a matriz simetrica positiva
*/
void genSimetricaPositiva(real_t ***A, real_t **b, int n, int k, int* new_k,
                          real_t ***ASP, real_t **bsp, rtime_t *tempo)
{
    *tempo = timestamp();

    // Quando A tem k diagonais (k = 2m+1), A^T*A terá 2k-1 diagonais
    // mas limitado ao máximo de 2n-1
    *new_k = 2 * k - 1;
    if (*new_k > 2*n - 1) {
        *new_k = 2*n - 1;
    }

    // Aloca memória para a matriz A^T * A (simétrica positiva definida)
    *ASP = malloc(sizeof(real_t*) * (*new_k));
    for (int i = 0; i < *new_k; i++) {
        (*ASP)[i] = calloc(n, sizeof(real_t));
    }

    // Aloca memória para o vetor A^T * b
    *bsp = calloc(n, sizeof(real_t));

    int m = k/2;

    // Calcula A^T * b
    // (A^T * b)[i] = soma_j A^T[i][j] * b[j] = soma_j A[j][i] * b[j]
    // A[j][i] está em A[diag_idx][j] onde diag_idx = (i - j) + m
    for (int i = 0; i < n; i++) {
        real_t sum = 0.0;
        for (int j = 0; j < n; j++) {
            int diag_offset = i - j;  // A[j][i] tem offset (i - j) relativo à linha j
            if (abs(diag_offset) <= m) {
                int diag_idx = diag_offset + m;
                sum += (*A)[diag_idx][j] * (*b)[j];
            }
        }
        (*bsp)[i] = sum;
    }

    // Calcula A^T * A (matriz com k+2 diagonais)
    // (A^T * A)[i][j] = soma_l A^T[i][l] * A[l][j] = soma_l A[l][i] * A[l][j]
    int m_new = *new_k / 2;
    
    for (int i = 0; i < n; i++) {
        for (int d = 0; d < *new_k; d++) {
            int diag_offset = d - m_new;
            int j = i + diag_offset;
            
            if (j < 0 || j >= n) continue;
            
            real_t sum = 0.0;
            
            // Soma sobre todas as linhas l onde tanto A[l][i] quanto A[l][j] existem
            // A[l][i] existe se |i - l| <= m
            // A[l][j] existe se |j - l| <= m
            int l_min = 0;
            int l_max = n - 1;
            
            // Limita l para que A[l][i] exista
            if (i - m > l_min) l_min = i - m;
            if (i + m < l_max) l_max = i + m;
            
            // Limita l para que A[l][j] exista
            if (j - m > l_min) l_min = j - m;
            if (j + m < l_max) l_max = j + m;
            
            for (int l = l_min; l <= l_max; l++) {
                // Acessa A[l][i] -> está em A[diag_i][l] onde diag_i = (i - l) + m
                int offset_i = i - l;
                int diag_i = offset_i + m;
                
                // Acessa A[l][j] -> está em A[diag_j][l] onde diag_j = (j - l) + m
                int offset_j = j - l;
                int diag_j = offset_j + m;
                
                sum += (*A)[diag_i][l] * (*A)[diag_j][l];
            }
            
            (*ASP)[d][i] = sum;
        }
    }

    *tempo = timestamp() - *tempo;
}

/** Decomposicao DLU 
 * @param *A matriz de coeficientes KxN
 * @param n ordem do sistema linear
 * @param k quantidade de diagonais
 * @param *D vetor com coeficientes da diagonal principal
 * @param **L matriz triangular inferior
 * @param **U matriz triangular superior
 * @param *tempo tempo utilizado para o calculo
 */
void geraDLU (real_t ***A, int n, int k,
              real_t **D, real_t ***L, real_t ***U, rtime_t *tempo)
{
    *tempo = timestamp();

    int m = k/2;

    // Aloca memória para D (diagonal principal)
    *D = malloc(sizeof(real_t) * n);

    // Aloca memória para L (triangular inferior)
    *L = malloc(sizeof(real_t*) * m);
    for (int i = 0; i < m; i++)
        (*L)[i] = calloc(n, sizeof(real_t));

    // Aloca memória para U (triangular superior)
    *U = malloc(sizeof(real_t*) * m);  
    for (int i = 0; i < m; i++)
        (*U)[i] = calloc(n, sizeof(real_t));

    // Extrai diagonal principal (D)
    for (int j = 0; j < n; j++)
        (*D)[j] = (*A)[m][j];  // Diagonal principal está na posição m

    // Extrai parte triangular inferior (L)
    // L[i] armazena a (i+1)-ésima diagonal inferior de A
    // Em A k-diagonal: A[m-i-1] contém a (i+1)-ésima diagonal inferior
    for (int i = 0; i < m; i++) {
        int diag_idx = m - i - 1;
        for (int j = 0; j < n; j++)
            (*L)[i][j] = (*A)[diag_idx][j];
    }

    // Extrai parte triangular superior (U)
    // U[i] armazena a (i+1)-ésima diagonal superior de A
    // Em A k-diagonal: A[m+i+1] contém a (i+1)-ésima diagonal superior
    for (int i = 0; i < m; i++) {
        int diag_idx = m + i + 1;
        for (int j = 0; j < n; j++)
            (*U)[i][j] = (*A)[diag_idx][j];
    }

    *tempo = timestamp() - *tempo;
}

/**
 * Gera matriz de pré-condicionamento M
 * @param *D vetor da diagonal principal
 * @param *L matriz triangular inferior
 * @param *U matriz triangular superior
 * @param w parâmetro do pré-condicionador 
 * @param n ordem do sistema linear
 * @param k quantidade de diagonais
 * @param **M matriz de pré-condicionamento (armazenada como k-diagonal)
 * @param *tempo tempo utilizado para o calculo
 */
void geraPreCond(real_t **D, real_t ***L, real_t ***U, real_t w, int n, int k,
                 real_t ***M, rtime_t *tempo)
{
    *tempo = timestamp();

    int m = k/2;

    if (w == -1.0) { 
        // Sem pré-condicionador: M = I (matriz identidade)
        *M = malloc(sizeof(real_t*) * k);
        for (int i = 0; i < k; i++)
            (*M)[i] = calloc(n, sizeof(real_t));
        
        for (int j = 0; j < n; j++)
            (*M)[m][j] = 1.0;

    } else if (w == 0.0) { 
        // Pré-condicionador de Jacobi: M = D
        *M = malloc(sizeof(real_t*) * k);
        for (int i = 0; i < k; i++)
            (*M)[i] = calloc(n, sizeof(real_t));

        // Copia diagonal principal (inversa de D)
        for (int j = 0; j < n; j++) {
            if ((*D)[j] != 0.0) {
                (*M)[m][j] = 1.0 / (*D)[j];
            } else {
                fprintf(stderr, "Erro: elemento diagonal zero na posição %d!\n", j);
                exit(-1);
            }
        }

    } else if (w >= 1.0 && w < 2.0) { 
        // Gauss-Seidel (w=1.0) ou SSOR (1.0 < w < 2.0)
        // M = (D + ωL)D^-1(D + ωU)
        // D, L e U são armazenados para resolver sistemas triangulares
        
        *M = malloc(sizeof(real_t*) * k);
        for (int i = 0; i < k; i++)
            (*M)[i] = calloc(n, sizeof(real_t));

        // Armazena diagonal principal D
        for (int j = 0; j < n; j++) {
            if ((*D)[j] != 0.0) {
                (*M)[m][j] = (*D)[j];
            } else {
                fprintf(stderr, "Erro: elemento diagonal zero na posição %d!\n", j);
                exit(-1);
            }
        }

        // Armazena parte triangular inferior L
        // L[0] contém a primeira diagonal inferior
        // M deve ter: M[m-1] = primeira diag inf, M[m-2] = segunda diag inf, etc
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                (*M)[m - 1 - i][j] = (*L)[i][j];
            }
        }

        // Armazena parte triangular superior U
        // U[0] contém a primeira diagonal superior
        // M deve ter: M[m+1] = primeira diag sup, M[m+2] = segunda diag sup, etc
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                (*M)[m + 1 + i][j] = (*U)[i][j];
            }
        }

    } else {
        fprintf(stderr, "Erro: valor de w inválido: %lf\n", w);
        exit(-1);
    }

    *tempo = timestamp() - *tempo;
}


/** Calcula o residuo do sistema linear 
 * @param **A matriz de coeficientes
 * @param *b vetor de termos independentes
 * @param *X vetor das incognitas
 * @param n ordem do sistema linear
 * @param k quantidade de diagonais
 * @param *tempo tempo utilizado para o calculo
 * @return residuo (norma euclidiana ||b - Ax||_2)
 */
real_t calcResiduoSL (real_t ***A, real_t **b, real_t **X,
                      int n, int k, rtime_t *tempo)
{
    *tempo = timestamp();

    real_t *r = malloc(n * sizeof(real_t));

    // Inicializa r = b
    for (int i = 0; i < n; i++) 
        r[i] = (*b)[i];

    int m = k/2;

    // Calcula r = b - Ax
    for (int diag = 0; diag < k; diag++) {
        int diag_offset = diag - m;

        for (int j = 0; j < n; j++) {
            if ((*A)[diag][j] != 0.0) {
                // Para matriz k-diagonal: A[diag][j] está na posição (j, j+diag_offset)
                int row_idx = j;
                int col_idx = j + diag_offset;

                if (col_idx >= 0 && col_idx < n) {
                    r[row_idx] -= (*A)[diag][j] * (*X)[col_idx];
                }
            }
        }
    }

    // Calcula norma euclidiana ||r||_2 = sqrt(sum(r[i]^2))
    real_t residuo = 0.0;
    for (int i = 0; i < n; i++) {
        residuo += r[i] * r[i];
    }
    residuo = sqrt(residuo);

    free(r);
    *tempo = timestamp() - *tempo;

    return residuo;
}

/**
 * Multiplica matriz k-diagonal por vetor
 * @param ***A matriz k-diagonal (k x n)
 * @param *x vetor de tamanho n
 * @param *result vetor com o resultado da multiplicacao
 * @param n quantidade de linhas da matriz A e tamanho do vetor X
 * @param k quantidade de colunas da matriz A
 */
void multiplicaMatrizVetor(real_t ***A, real_t *x, real_t *result, int n, int k)
{
    int m = k/2;

    // Inicializa resultado com zeros
    for (int i = 0; i < n; i++) 
        result[i] = 0.0;

    // Multiplica cada diagonal
    for (int diag = 0; diag < k; diag++) {
        int diag_offset = diag - m;

        for (int j = 0; j < n; j++) {
            if ((*A)[diag][j] != 0.0) {
                // Para matriz k-diagonal: A[diag][j] está na posição (j, j+diag_offset)
                int row_idx = j;
                int col_idx = j + diag_offset;

                if (col_idx >= 0 && col_idx < n) {
                    result[row_idx] += (*A)[diag][j] * x[col_idx];
                }
            }
        }
    }
}

/**
 * Aplica pré-condicionador: resolve M*v = r
 * @param ***M matriz pre-condicionada de tamanho kxn
 * @param *r vetor de termos independentes
 * @param *v vetor solução
 * @param n tamanho dos vetores e quantidade de linhas da matriz M
 * @param k quantidade de colunas da matriz M
 */
void aplicaPreCondicionador(real_t ***M, real_t *r, real_t *v, int n, int k, real_t w)
{
    int m = k/2;

    if (w == -1.0) {
        // M = I: v = r
        for (int i = 0; i < n; i++) {
            v[i] = r[i];
        }
        
    } else if (w == 0.0) {
        // Jacobi: M = D^-1 | v = M * r
        multiplicaMatrizVetor(M, r, v, n, k);
        
    } else if (w >= 1.0 && w < 2.0) {
        // Pré-condicionador SSOR: resolve M^-1 * r
        // M = (D + ωL) * D^-1 * (D + ωU)
        // Aplicação via substituições triangulares:
        // 1. Forward:  (D + ωL) * z = r
        // 2. Backward: (D + ωU) * v = D * z
        
        real_t *z = malloc(n * sizeof(real_t));
        
        // Resolve (D + ωL) * z = r
        // Equivalente a: D*z + ωL*z = r
        //                z = D^-1 * (r - ωL*z)
        for (int i = 0; i < n; i++) {
            real_t sum = r[i];
            
            // Subtrai ωL * z das diagonais inferiores
            for (int diag = 0; diag < m; diag++) {
                int diag_offset = diag - m;  // Negativo
                int j = i + diag_offset;     // j < i
                
                if (j >= 0 && (*M)[diag][i] != 0.0) {
                    sum -= w * (*M)[diag][i] * z[j];
                }
            }
            
            real_t d_ii = (*M)[m][i];
            if (ABS(d_ii) < 1e-14) {
                fprintf(stderr, "Erro: diagonal muito pequena em i=%d (d=%g)\n", i, d_ii);
                free(z);
                exit(-1);
            }
            z[i] = sum / d_ii;
        }
        
        // Resolve (D + ωU) * v = D * z
        // Equivalente a: D*v + ωU*v = D*z
        //                v = D^-1 * (D*z - ωU*v)
        for (int i = n - 1; i >= 0; i--) {
            real_t d_ii = (*M)[m][i];
            real_t sum = d_ii * z[i];
            
            // Subtrai ωU * v das diagonais superiores
            for (int diag = m + 1; diag < k; diag++) {
                int diag_offset = diag - m;  // Positivo
                int j = i + diag_offset;     // j > i
                
                if (j < n && (*M)[diag][i] != 0.0) {
                    sum -= w * (*M)[diag][i] * v[j];
                }
            }
            
            if (ABS(d_ii) < 1e-14) {
                fprintf(stderr, "Erro: diagonal muito pequena em i=%d (d=%g)\n", i, d_ii);
                free(z);
                exit(-1);
            }
            v[i] = sum / d_ii;
        }
        
        free(z);
        
    } else {
        fprintf(stderr, "Erro: valor de w=%lf não suportado em aplicaPreCondicionador\n", w);
        exit(-1);
    }
}

/**
 * Calcula produto interno entre dois vetores
 * @param *a vetor de tamanho n
 * @param *b segundo vetor de tamanho n
 * @param n tamanho dos vetores a e b
 * @return produto interno (a^T * b)
 */
real_t produtoInterno(real_t *a, real_t *b, int n)
{
    real_t result = 0.0;
    for (int i = 0; i < n; i++) {
        result += a[i] * b[i];
    }
    return result;
}

/**
 * Calcula norma máxima da diferença entre dois vetores
 * @param *x_old vetor antigo
 * @param *x_new vetor novo
 * @param n tamanho dos vetores
 * @return norma máxima (||x_new - x_old||)
 */
real_t normaMaxima(real_t *x_old, real_t *x_new, int n)
{
    real_t max_diff = 0.0;
    for (int i = 0; i < n; i++) {
        real_t diff = ABS(x_new[i] - x_old[i]);
        if (diff > max_diff) {
            max_diff = diff;
        }
    }
    return max_diff;
}

/**
 * Método dos Gradientes Conjugados com pré-condicionador
 * @param ***A matriz k-diagonal do sistema linear
 * @param *b vetor de termos independentes
 * @param *x vetor solução
 * @param ***M matriz de pré-condicionamento
 * @param n ordem do sistema linear
 * @param k quantidade de diagonais
 * @param epsilon critério de convergência
 * @param maxit número máximo de iterações
 * @param *norma norma máxima da diferença entre iterações consecutivas
 * @param *tempo_iter tempo médio por iteração
 * @return número de iterações realizadas
 */
int gradientesConjugados(real_t ***A, real_t *b, real_t *x, real_t ***M, int n, int k, 
                         real_t epsilon, int maxit, real_t *norma, rtime_t *tempo_iter, real_t w)
{
    rtime_t tempo_inicio = timestamp();

    // Aloca vetores auxiliares
    real_t *r = malloc(n * sizeof(real_t));      // resíduo
    real_t *v = malloc(n * sizeof(real_t));      // direção de busca
    real_t *y = malloc(n * sizeof(real_t));      // pré-condicionado com residuo
    real_t *z = malloc(n * sizeof(real_t));      // A*v
    real_t *x_old = malloc(n * sizeof(real_t));

    // Inicializa x = 0
    for (int i = 0; i < n; i++) {
        x[i] = 0.0;
        x_old[i] = 0.0;
    }

    // Calcula resíduo inicial: r0 = b - A*x0 = b (pois x0 = 0)
    for (int i = 0; i < n; i++) {
        r[i] = b[i];
    }

    aplicaPreCondicionador(M, r, y, n, k, w); // y = M^-1 * r
    
    // Inicializa v = y
    for (int i = 0; i < n; i++) {
        v[i] = y[i];
    }

    real_t aux = produtoInterno(y, r, n);  // y^T * r

    int iter;
    for (iter = 0; iter < maxit; iter++) {
        multiplicaMatrizVetor(A, v, z, n, k); // z = A*v

        real_t vtz = produtoInterno(v, z, n); // v^T * z
        if (ABS(vtz) < 1e-14) {
            // Convergência numérica - o denominador tornou-se muito pequeno
            // Isso indica que está muito próximo da solução ótima
            if (iter == 0) {
                fprintf(stderr, "Erro: problema na primeira iteração - matriz mal condicionada\n");
                break;
            }
            break;
        }
        real_t s = aux / vtz;

        // Atualiza solução: x_(k+1) = x_(k) + s * v
        for (int i = 0; i < n; i++) {
            x[i] += s * v[i];
        }

        // Atualiza resíduo: r = r - s * z
        for (int i = 0; i < n; i++) {
            r[i] -= s * z[i];
        }

        // Aplica pré-condicionador: M * y = r
        aplicaPreCondicionador(M, r, y, n, k, w);

        // Verifica convergência usando norma L2
        *norma = normaMaxima(x_old, x, n);
        if (*norma < epsilon) {
            iter++;
            break;
        }

        real_t aux1 = produtoInterno(y, r, n);
        real_t m = aux1 / aux;
        aux = aux1;

        // Atualiza direção: v = y + m * v
        for (int i = 0; i < n; i++) {
            v[i] = y[i] + m * v[i];
        }

        // Copia x para x_old
        for (int i = 0; i < n; i++)
            x_old[i] = x[i];
    }

    rtime_t tempo_total = timestamp() - tempo_inicio;
    *tempo_iter = (iter > 0) ? tempo_total / iter : 0.0;

    // Libera memória
    free(r);
    free(z);
    free(v);
    free(y);
    free(x_old);

    return iter;
}
