/*
 * Nome: Kauê de Amorim Silva
 * GRR: 20244719
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
 * 
 * Discretiza a EDO y'' + py' + qy = r(x) usando diferenças finitas centradas
 * em uma malha uniforme, gerando um sistema Ax = b onde A é tridiagonal.
 * 
 * A discretização usa:
 * - y''(xi) ≈ (y[i-1] - 2*y[i] + y[i+1])/h²
 * - y'(xi) ≈ (y[i+1] - y[i-1])/(2*h)
 * 
 * @param edo Ponteiro para a estrutura EDO com os parâmetros da equação
 * @return Ponteiro para o sistema tridiagonal, ou NULL se erro de alocação
 */
Tridiag *genTridiag (EDo *edoeq);

/**
 * @brief Resolve sistema tridiagonal usando método de Gauss-Seidel (versão otimizada)
 * 
 * Implementação otimizada que aproveita a estrutura tridiagonal específica,
 * evitando verificações condicionais desnecessárias dentro dos loops principais.
 * 
 * @param sl Ponteiro para o sistema tridiagonal
 * @param Y Vetor solução (entrada: chute inicial, saída: solução)
 * @param maxiter Número máximo de iterações permitidas
 * @param norma Ponteiro onde será armazenada a norma L2 do resíduo final
 * @return Número de iterações realizadas até convergência ou limite máximo
 */
int gaussSeidel_3Diag(Tridiag *sl, real_t *Y, int maxiter, real_t *norma);

/**
 * @brief Calcula a norma L2 do resíduo para sistemas tridiagonais
 * 
 * Função otimizada que calcula ||b - Ax||₂ aproveitando a estrutura
 * tridiagonal da matriz, tratando explicitamente os casos especiais
 * da primeira e última linha.
 * 
 * @param sl Ponteiro para o sistema tridiagonal  
 * @param Y Vetor solução atual
 * @return Norma L2 do resíduo (||b - Ax||₂)
 */
real_t normaL2_3Diag(Tridiag *sl, real_t *Y);

/**
 * @brief Imprime a solução do sistema em formato horizontal
 * 
 * Formata a saída da solução em uma única linha com valores separados por espaços,
 * usando o formato de precisão definido em FORMAT (15 casas decimais).
 * 
 * @param sol Vetor com a solução do sistema
 * @param n Tamanho do vetor solução
 */
void prnSolucao(real_t *sol, int n);

/**
 * @brief Imprime a matriz aumentada do sistema linear na saída padrão
 * @param edoeq Ponteiro para a estrutura EDO
 */
void prnEDOsl (EDo *edoeq);

#endif // __EQDIFF_H__

