# LEIAME - Resolução de Equações Diferenciais Ordinárias

**Aluno:** Kauê de Amorim Silva  
**GRR:** 20244719

## Descrição do Programa

Este programa resolve Equações Diferenciais Ordinárias (EDO) de segunda ordem da forma:

```
y'' + py' + qy = r(x)
```

onde r(x) = r1*x + r2*x² + r3*cos(x) + r4*exp(x)

O programa utiliza discretização por diferenças finitas e resolve o sistema linear tridiagonal resultante através do método iterativo de Gauss-Seidel.

## Compilação

Para compilar o programa, execute:
```bash
make
```

## Execução

Execute o programa fornecendo os dados através da entrada padrão:
```bash
./resolveEDO < arquivo_entrada.dat
```

### Formato da Entrada

1. **Linha 1:** número de pontos da malha (n)
2. **Linha 2:** intervalo [a, b] 
3. **Linha 3:** condições de contorno y(a), y(b)
4. **Linha 4:** coeficientes p, q da EDO
5. **Linhas seguintes:** coeficientes r1, r2, r3, r4 de r(x) (uma EDO por linha)

### Formato da Saída

Para cada EDO processada, o programa imprime:
- A ordem do sistema e matriz aumentada
- A solução numérica (valores de y nos pontos internos)
- Número de iterações do método de Gauss-Seidel
- Norma L2 do resíduo final
- Tempo de execução da resolução

## Limitações do Programa

### 1. Convergência do Método Iterativo
- O método de Gauss-Seidel pode não convergir para certas combinações de coeficientes p e q
- A convergência depende da diagonal dominância da matriz tridiagonal gerada
- Máximo de 100 iterações por resolução (definido na constante MAXIT)

### 2. Critério de Parada
- Tolerância fixa de 1.0e-5 para a norma L2 do resíduo
- Não há controle dinâmico da tolerância baseado no problema

### 3. Limitações Numéricas
- Usa arredondamento para baixo (FE_DOWNWARD) conforme especificação
- Pode haver perda de precisão em casos com números muito pequenos ou muito grandes
- Precisão limitada a double (15-17 dígitos significativos)

### 4. Restrições da EDO
- Funciona apenas para EDOs lineares de segunda ordem
- Formato específico para r(x): combinação linear de x, x², cos(x), exp(x)

### 5. Limitações de Entrada
- Não há validação robusta dos dados de entrada
- Programa pode falhar com entradas inválidas ou malformadas
- Não trata casos degenerados (n ≤ 0, a ≥ b, etc.)

### 6. Gestão de Memória
- Falhas de alocação podem causar terminação abrupta
- Não há recuperação em caso de erro de memória

## Como Utilizar o likwid.sh

O script `likwid.sh` é usado para medir o desempenho do programa utilizando a ferramenta LIKWID.

### Pré-requisitos
- LIKWID instalado no sistema
- Permissões para modificar configurações de CPU
- Variáveis de ambiente LIKWID_INCLUDE e LIKWID_LIB configuradas

### Execução
```bash
./likwid.sh
```

### O que o Script Faz

1. **Configuração de Performance:**
   - Define o governador da CPU para "performance" no núcleo 3
   - Garante frequência constante durante as medições

2. **Compilação:**
   - Remove executáveis antigos (`make purge`)
   - Recompila o programa com flags LIKWID

3. **Medição:**
   - Executa o programa com `likwid-perfctr`
   - Mede métricas FLOPS_DP (operações de ponto flutuante dupla precisão)
   - Salva resultados em arquivo CSV temporário

4. **Resultados:**
   - Extrai e exibe contadores de instruções de ponto flutuante
   - Remove arquivos temporários

5. **Restauração:**
   - Retorna governador da CPU para "powersave"

### Personalização

Para modificar as medições, edite as variáveis no início do script:
- `METRICA`: tipo de métrica a medir (padrão: FLOPS_DP)
- `NUCLEO`: núcleo da CPU a usar (padrão: 3)
- `CSV`: nome do arquivo de saída temporário

### Exemplo de Uso
```bash
# Torna o script executável (se necessário)
chmod +x likwid.sh

# Executa medição de performance
./likwid.sh < edos.dat
```

O script é especialmente útil para:
- Análise de performance do método de Gauss-Seidel
- Comparação de eficiência entre diferentes tamanhos de problema
- Validação de otimizações de código