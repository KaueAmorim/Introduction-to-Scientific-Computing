#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h> // Para uso de função 'memset()'

#include "matriz.h"

/**
 * Função que gera valores para para ser usado em uma matriz
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
*  @return valor gerado para a posição i,j
  */
static inline real_t generateRandomA( unsigned int i, unsigned int j)
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return ( (i==j) ? (real_t)(BASE<<1) : 1.0 )  * (real_t)random() * invRandMax;
}

/**
 * Função que gera valores aleatórios para ser usado em um vetor
 * @return valor gerado
 *
 */
static inline real_t generateRandomB( )
{
  static real_t invRandMax = 1.0 / (real_t)RAND_MAX;
  return (real_t)(BASE<<2) * (real_t)random() * invRandMax;
}



/* ----------- FUNÇÕES ---------------- */

/**
 *  Funcao geraMatRow: gera matriz como vetor único, 'row-oriented'
 * 
 *  @param m     número de linhas da matriz
 *  @param n     número de colunas da matriz
 *  @param zerar se 0, matriz  tem valores aleatórios, caso contrário,
 *               matriz tem valores todos nulos
 *  @return  ponteiro para a matriz gerada
 *
 */

MatRow geraMatRow (int m, int n, int zerar)
{
  MatRow matriz = (real_t *) malloc(m*n*sizeof(real_t));

  if (matriz) {
    if (zerar)
      memset(matriz,0,m*n*sizeof(real_t));
    else
      for (int i=0; i < m; ++i)
	for (int j=0; j < n; ++j)
	  matriz[i*n + j] = generateRandomA(i, j);
  }
  
  return (matriz);
}


/**
 *  Funcao geraVetor: gera vetor de tamanho 'n'
 * 
 *  @param n  número de elementos do vetor
 *  @param zerar se 0, vetor  tem valores aleatórios, caso contrário,
 *               vetor tem valores todos nulos
 *  @return  ponteiro para vetor gerado
 *
 */

Vetor geraVetor (int n, int zerar)
{
  Vetor vetor = (real_t *) malloc(n*sizeof(real_t));

  if (vetor) {
    if (zerar)
      memset(vetor,0,n*sizeof(real_t));
    else
      for (int i=0; i < n; ++i)
	vetor[i] = generateRandomB();
  }
  
  return (vetor);
}

/**
 *  \brief: libera vetor
 * 
 *  @param  ponteiro para vetor
 *
 */
void liberaVetor (void *vet)
{
	free(vet);
}


/**
 *  Funcao multMatVet:  Efetua multiplicacao entre matriz 'mxn' por vetor
 *                       de 'n' elementos
 *  @param mat matriz 'mxn'
 *  @param m número de linhas da matriz
 *  @param n número de colunas da matriz
 *  @param res vetor que guarda o resultado. Deve estar previamente alocado e com
 *             seus elementos inicializados em 0.0 (zero)
 *  @return vetor de 'm' elementos
 *
 */

void multMatVet (MatRow mat, Vetor v, int m, int n, Vetor res)
{
    
  /* Efetua a multiplicação */
  if (res) {
    for (int i=0; i < m; ++i)
      for (int j=0; j < n; ++j)
        res[i] += mat[n*i + j] * v[j];
  }
}

/**
 *  Funcao multMatVet_otim:  Versão otimizada da multiplicacao matriz-vetor
 *  Utiliza loop unrolling e melhor localidade de cache
 */
void multMatVet_otim (MatRow mat, Vetor v, int m, int n, Vetor res)
{
  if (res) {
    for (int i=0; i < m; ++i) {
      register real_t sum = 0.0;
      register int pos = i * n;
      
      // Loop unrolling por 4
      int j;
      for (j = 0; j < n-3; j += 4) {
        sum += mat[pos + j] * v[j] + 
               mat[pos + j + 1] * v[j + 1] +
               mat[pos + j + 2] * v[j + 2] +
               mat[pos + j + 3] * v[j + 3];
      }
      
      // Processa elementos restantes
      for (; j < n; ++j) {
        sum += mat[pos + j] * v[j];
      }
      
      res[i] = sum;
    }
  }
}


/**
 *  Funcao multMatMat: Efetua multiplicacao de duas matrizes 'n x n' 
 *  @param A matriz 'n x n'
 *  @param B matriz 'n x n'
 *  @param n ordem da matriz quadrada
 *  @param C   matriz que guarda o resultado. Deve ser previamente gerada com 'geraMatPtr()'
 *             e com seus elementos inicializados em 0.0 (zero)
 *
 */

void multMatMat (MatRow A, MatRow B, int n, MatRow C)
{

  /* Efetua a multiplicação */
  for (int i=0; i < n; ++i)
    for (int j=0; j < n; ++j)
      for (int k=0; k < n; ++k)
	C[i*n+j] += A[i*n+k] * B[k*n+j];
}

/**
 *  Funcao multMatMat_otim: Versão otimizada da multiplicacao matriz-matriz
 *  Utiliza blocking/tiling para melhor uso de cache
 */
void multMatMat_otim (MatRow A, MatRow B, int n, MatRow C)
{
  const int BLOCK_SIZE = 64;  // Tamanho do bloco otimizado para cache
  
  for (int ii = 0; ii < n; ii += BLOCK_SIZE) {
    for (int jj = 0; jj < n; jj += BLOCK_SIZE) {
      for (int kk = 0; kk < n; kk += BLOCK_SIZE) {
        
        // Limites dos blocos
        int i_max = (ii + BLOCK_SIZE < n) ? ii + BLOCK_SIZE : n;
        int j_max = (jj + BLOCK_SIZE < n) ? jj + BLOCK_SIZE : n;
        int k_max = (kk + BLOCK_SIZE < n) ? kk + BLOCK_SIZE : n;
        
        // Multiplicação dentro do bloco
        for (int i = ii; i < i_max; ++i) {
          for (int j = jj; j < j_max; ++j) {
            register real_t sum = C[i*n + j];
            for (int k = kk; k < k_max; ++k) {
              sum += A[i*n + k] * B[k*n + j];
            }
            C[i*n + j] = sum;
          }
        }
      }
    }
  }
}


/**
 *  Funcao prnMat:  Imprime o conteudo de uma matriz em stdout
 *  @param mat matriz
 *  @param m número de linhas da matriz
 *  @param n número de colunas da matriz
 *
 */

void prnMat (MatRow mat, int m, int n)
{
  for (int i=0; i < m; ++i) {
    for (int j=0; j < n; ++j)
      printf(DBL_FIELD, mat[n*i + j]);
    printf("\n");
  }
  printf(SEP_RES);
}

/**
 *  Funcao prnVetor:  Imprime o conteudo de vetor em stdout
 *  @param vet vetor com 'n' elementos
 *  @param n número de elementos do vetor
 *
 */

void prnVetor (Vetor vet, int n)
{
  for (int i=0; i < n; ++i)
    printf(DBL_FIELD, vet[i]);
  printf(SEP_RES);
}

