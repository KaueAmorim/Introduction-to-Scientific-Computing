#!/usr/bin/gnuplot -c
set grid
set style data point
set style function line
set style line 1 lc 3 pt 7 ps 0.3
set boxwidth 1
set xtics
set xrange [64:4500]
set xlabel  "N (bytes)"

set logscale x
# set logscale y

#
# TEMPO
#
ARQ=ARG1."/Tempos.csv"
set key left top
set logscale y
set ylabel  "Tempo (ms)"
set title   "Tempo"
set terminal qt 0 title "Tempos"
plot ARQ using 1:2 title "MatVet" lc rgb "green" with linespoints, \
     '' using 1:3 title "MatVet-uj" lc rgb "red" with linespoints, \
     '' using 1:4 title "MatMat" lc rgb "magenta" with linespoints, \
     '' using 1:5 title "MatMat-uj" lc rgb "cyan" with linespoints

pause -1

## set datafile separator {whitespace | tab | comma | "<chars>"}
set datafile separator comma

#
# FLOPS_DP
#
ARQ=ARG1."/FLOPS_DP.csv"
set key right top
unset logscale y
set ylabel  "FLOPS DP [MFlops/s]"
set title   "FLOPS DP"
set terminal qt 1 title "FLOPS DP"
plot ARQ using 1:2 title "MatVet" lc rgb "green" with linespoints, \
     '' using 1:3 title "MatVet-uj" lc rgb "red" with linespoints, \
     '' using 1:4 title "MatMat" lc rgb "magenta" with linespoints, \
     '' using 1:5 title "MatMat-uj" lc rgb "cyan" with linespoints

pause -1

#
# ENERGY
#
ARQ=ARG1."/ENERGY.csv"
set key center top
unset logscale y
set ylabel  "Energia [J]"
set title   "Energia"
set terminal qt 2 title "Energia"
plot ARQ using 1:2 title "MatVet" lc rgb "green" with linespoints, \
     '' using 1:3 title "MatVet-uj" lc rgb "red" with linespoints, \
     '' using 1:4 title "MatMat" lc rgb "magenta" with linespoints, \
     '' using 1:5 title "MatMat-uj" lc rgb "cyan" with linespoints

pause -1

#
# L3CACHE
#
ARQ=ARG1."/L3CACHE.csv"
set key left top
unset logscale y
set ylabel  "L3 miss ratio"
set title   "L3 miss ratio"
set terminal qt 3 title "L3 miss ratio"
plot ARQ using 1:2 title "MatVet" lc rgb "green" with linespoints, \
     '' using 1:3 title "MatVet-uj" lc rgb "red" with linespoints, \
     '' using 1:4 title "MatMat" lc rgb "magenta" with linespoints, \
     '' using 1:5 title "MatMat-uj" lc rgb "cyan" with linespoints

pause -1


# Gerando figura PNG
#set terminal png
#set output "plot_FLOPS_DP.png"
#plot 'plot_FLOPS_DP.csv' using 1:2 title "MatPtrVet" with linespoints, \
#     '' using 1:3 title "MatVet" with linespoints, \
#     '' using 1:4 title "MatMatPtr" with linespoints, \
#     '' using 1:5 title "MatMat" with linespoints
#replot
#unset output

