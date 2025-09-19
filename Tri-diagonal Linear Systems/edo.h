/*
 * Nome: Kaue (substituir pelo seu nome)
 * GRR: XXXXXXXXX (substituir pelo seu GRR)
 * 
 * Definições e estruturas para resolução de Equações Diferenciais Ordinárias (EDO)
 * usando métodos de diferenças finitas e sistemas lineares tridiagonais.
 */

#ifndef __EQDIFF_H__
#define __EQDIFF_H__

typedef double real_t;

#define FORMAT "%23.15e"

/**
 * @struct Tridiag
 * @brief Estrutura para representar um sistema linear tridiagonal
 * 
 * @param D Vetor da diagonal principal (tamanho n)
 * @param Di Vetor da diagonal inferior (tamanho n) 
 * @param Ds Vetor da diagonal superior (tamanho n)
 * @param B Vetor dos termos independentes (tamanho n)
 * @param n Ordem do sistema linear
 */
typedef struct {
  real_t *D, *Di, *Ds, *B;
  int n;
} Tridiag;

/**
 * @struct EDo
 * @brief Estrutura para representar uma Equação Diferencial Ordinária
 * 
 * Representa EDOs da forma: y'' + py' + qy = r(x)
 * onde r(x) = r1*x + r2*x² + r3*cos(x) + r4*exp(x)
 * 
 * @param n Número de pontos internos na malha
 * @param a Limite inferior do intervalo
 * @param b Limite superior do intervalo  
 * @param ya Condição de contorno y(a)
 * @param yb Condição de contorno y(b)
 * @param p Coeficiente de y' na EDO
 * @param q Coeficiente de y na EDO
 * @param r1 Coeficiente do termo x em r(x)
 * @param r2 Coeficiente do termo x² em r(x)
 * @param r3 Coeficiente do termo cos(x) em r(x)
 * @param r4 Coeficiente do termo exp(x) em r(x)
 */
typedef struct {
  int n;
  real_t a, b;
  real_t ya, yb;
  real_t p, q, r1, r2, r3, r4;
} EDo;

// Declarações de funções

/**
 * @brief Gera um sistema linear tridiagonal a partir de uma EDO
 * @param edoeq Ponteiro para a estrutura EDO
 * @return Ponteiro para o sistema tridiagonal gerado, ou NULL em caso de erro
 * 
 * Códigos de erro:
 * - NULL: Falha na alocação de memória
 */
Tridiag *genTridiag (EDo *edoeq);

/**
 * @brief Imprime a matriz aumentada do sistema linear na saída padrão
 * @param edoeq Ponteiro para a estrutura EDO
 */
void prnEDOsl (EDo *edoeq);

/**
 * @brief Resolve um sistema tridiagonal usando o método de Gauss-Seidel
 * @param sl Ponteiro para o sistema tridiagonal
 * @param sol Vetor onde será armazenada a solução (deve estar alocado)
 * @param norma_residuo Ponteiro onde será armazenada a norma L2 do resíduo final
 * @return Número de iterações realizadas
 * 
 * Critérios de parada:
 * - Número máximo de iterações (100)
 * - Norma L2 do resíduo ≤ 10⁻⁵
 * 
 * Códigos de erro/retorno:
 * - > 0: Número de iterações até convergência
 * - 100: Máximo de iterações atingido (pode não ter convergido)
 */
int gaussSeidel(Tridiag *sl, real_t *sol, real_t *norma_residuo);

/**
 * @brief Imprime a solução do sistema em formato horizontal
 * @param sol Vetor solução
 * @param n Tamanho do vetor
 */
void prnSolucao(real_t *sol, int n);

/**
 * @brief Resolve sistema tridiagonal usando algoritmo de Thomas (eliminação direta)
 * @param sl Ponteiro para sistema tridiagonal
 * @param sol Vetor onde será armazenada a solução (deve estar alocado)
 * 
 * Códigos de erro:
 * - Função pode falhar se algum elemento da diagonal principal for zero
 * - Não há verificação explícita de erro, assume sistema bem condicionado
 */
void solveTridiag(Tridiag *sl, real_t *sol);

#endif // __EQDIFF_H__

