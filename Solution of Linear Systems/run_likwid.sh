#!/bin/bash

# Remove todos os objetos e o executável
make purge

# Compila com LIKWID e -O0, forçando recompilação
make clean
make CFLAGS="-O0 -DLIKWID_PERFMON" LFLAGS="-llikwid -lm" --always-make

# Executa com LIKWID e filtra a saída
likwid-perfctr -C 0 -g FLOPS_DP ./perfSL 2>&1 | \
grep -E 'FP_ARITH_INST_RETIRED_SCALAR_DOUBLE|DP \\(MFLOP/s\\)'