#!/usr/bin/gnuplot

# Script para gerar gráficos de análise de desempenho

set terminal pngcairo size 1200,800
set key outside right
set grid
set logscale x

# Gráfico de Tempo
set output 'Dados/grafico_tempo.png'
set title 'Tempo de Execução vs Tamanho da Matriz'
set xlabel 'Tamanho da Matriz (N)'
set ylabel 'Tempo (ms)'
plot 'Dados/Tempos.csv' using 1:2 with linespoints title 'matVet', \
     'Dados/Tempos.csv' using 1:3 with linespoints title 'matVet_otim', \
     'Dados/Tempos.csv' using 1:4 with linespoints title 'matMat', \
     'Dados/Tempos.csv' using 1:5 with linespoints title 'matMat_otim'

# Gráfico de Cache Miss L3
set output 'Dados/grafico_l3cache.png'
set title 'Taxa de Cache Miss L3 vs Tamanho da Matriz'
set xlabel 'Tamanho da Matriz (N)'
set ylabel 'L3 Miss Ratio'
plot 'Dados/L3CACHE.csv' using 1:2 with linespoints title 'matVet', \
     'Dados/L3CACHE.csv' using 1:3 with linespoints title 'matVet_otim', \
     'Dados/L3CACHE.csv' using 1:4 with linespoints title 'matMat', \
     'Dados/L3CACHE.csv' using 1:5 with linespoints title 'matMat_otim'

# Gráfico de Energia
set output 'Dados/grafico_energia.png'
set title 'Consumo de Energia vs Tamanho da Matriz'
set xlabel 'Tamanho da Matriz (N)'
set ylabel 'Energia (J)'
plot 'Dados/ENERGY.csv' using 1:2 with linespoints title 'matVet', \
     'Dados/ENERGY.csv' using 1:3 with linespoints title 'matVet_otim', \
     'Dados/ENERGY.csv' using 1:4 with linespoints title 'matMat', \
     'Dados/ENERGY.csv' using 1:5 with linespoints title 'matMat_otim'

# Gráfico de FLOPS
set output 'Dados/grafico_flops.png'
set title 'Desempenho FLOPS vs Tamanho da Matriz'
set xlabel 'Tamanho da Matriz (N)'
set ylabel 'MFLOP/s'
plot 'Dados/FLOPS_DP.csv' using 1:2 with linespoints title 'matVet', \
     'Dados/FLOPS_DP.csv' using 1:3 with linespoints title 'matVet_otim', \
     'Dados/FLOPS_DP.csv' using 1:4 with linespoints title 'matMat', \
     'Dados/FLOPS_DP.csv' using 1:5 with linespoints title 'matMat_otim'

print "Gráficos gerados em Dados/"