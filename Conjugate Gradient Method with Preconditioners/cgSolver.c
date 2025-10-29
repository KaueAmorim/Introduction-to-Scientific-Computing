#include "sislin.h"

void imprimirSolucao(real_t **A, real_t *B, real_t *X, int n, int k) {
    printf("==================\n");
    imprimirSistemaLinear(&A, B, n, k);
    for (int i = 0; i < n; i++) {
        printf("%.16g ", X[i]);
    }
    printf("\n==================\n");
}

int main () {
    srandom(20252);

    int n, k, maxit;
    real_t w, e;

    if (scanf("%d %d %lf %d %lf", &n, &k, &w, &maxit, &e) != 5) {
        fprintf(stderr, "Erro na leitura dos parametros de entrada!\n");
        return -1;
    }
    
    // Validação dos parâmetros
    if (n <= 10) {
        fprintf(stderr, "Erro: n deve ser maior que 10!\n");
        return -1;
    }
    if (k <= 1 || (k % 2) == 0 || k > (2*n - 1)) {
        fprintf(stderr, "Erro: k deve ser maior que 1, menor ou igual a (2n - 1) e impar!\n");
        return -1;
    }
    if (w < -1.0 || w >= 2.0) {
        fprintf(stderr, "Erro: w deve estar no intervalo [-1, 2)!\n");
        return -1;
    }
    if (maxit <= 0) {
        fprintf(stderr, "Erro: maxit deve ser positivo!\n");
        return -1;
    }
    if (e <= 0.0) {
        fprintf(stderr, "Erro: epsilon deve ser positivo!\n");
        return -1;
    }

    // Gera matriz k-diagonal A e vetor b
    real_t **A, *b;
    criaKDiagonal(n, k, &A, &b);

    #ifdef __DEBUG__
    printf("Sistema linear original:\n");
    imprimirSistemaLinear(&A, b, n, k);
    #endif

    // Transforma em sistema simétrico positivo definido
    real_t **ASP, *bsp;
    int new_k = 0;
    rtime_t tempo_pc;
    genSimetricaPositiva(&A, &b, n, k, &new_k, &ASP, &bsp, &tempo_pc);

    #ifdef __DEBUG__
    printf("Sistema linear transformado:\n");
    imprimirSistemaLinear(&ASP, bsp, n, new_k);
    #endif

    // Decompõe matriz em D, L, U
    real_t *D, **L, **U;
    rtime_t tempo_dlu;
    geraDLU(&ASP, n, new_k, &D, &L, &U, &tempo_dlu);

    // Gera pré-condicionador
    real_t **M;
    rtime_t tempo_precond;
    geraPreCond(&D, &L, &U, w, n, new_k, &M, &tempo_precond);

    #ifdef __DEBUG__
    printf("Sistema linear com pre condicionador:\n");
    imprimirSistemaLinear(&M, bsp, n, new_k);
    #endif
    
    tempo_pc += tempo_dlu + tempo_precond;

    // Aloca vetor solução
    real_t *x = malloc(n * sizeof(real_t));

    // Resolve sistema usando Gradientes Conjugados
    real_t norma;
    rtime_t tempo_iter;
    int iteracoes = gradientesConjugados(&ASP, bsp, x, &M, n, new_k, e, maxit, &norma, &tempo_iter, w);

    // Calcula resíduo final do sistema original
    rtime_t tempo_residuo;
    real_t residuo = calcResiduoSL(&A, &b, &x, n, k, &tempo_residuo);
    
    #ifdef __DEBUG__
    rtime_t tempo_residuo_trans;
    real_t residuo_transformado = calcResiduoSL(&ASP, &bsp, &x, n, new_k, &tempo_residuo_trans);
    printf("Residuo do sistema transformado: %.16g\n", residuo_transformado);
    printf("Residuo do sistema original: %.16g\n", residuo);
    printf("Num iteracoes: %d\n", iteracoes);
    #endif

    // Saída dos resultados
    printf("%d\n", n);
    
    // Imprime vetor solução
    for (int i = 0; i < n; i++) {
        printf("%.16g", x[i]);
        if (i < n-1) printf(" ");
    }
    printf("\n");
    
    // Imprime norma, resíduo e tempos
    printf("%.8g\n", norma);
    printf("%.16g\n", residuo);
    printf("%.8g\n", tempo_pc);
    printf("%.8g\n", tempo_iter);
    printf("%.8g\n", tempo_residuo);

    // Verifica se convergiu
    if (iteracoes >= maxit) {
        fprintf(stderr, "Aviso: método não convergiu em %d iterações!\n", maxit);
    }

    #ifdef __DEBUG__
    printf("Sistema linear original (com solucao):\n");
    imprimirSolucao(A, b, x, n, k);
    #endif

    // Libera memória
    for (int i = 0; i < k; i++) {
        free(A[i]);
    }
    free(A);
    free(b);
    
    for (int i = 0; i < new_k; i++) {
        free(ASP[i]);
        free(M[i]);
    }
    free(ASP);
    free(M);
    free(bsp);
    free(D);
    
    int m = new_k/2;
    for (int i = 0; i < m; i++) {
        free(L[i]);
        free(U[i]);
    }
    free(L);
    free(U);
    free(x);

    return 0;
}
