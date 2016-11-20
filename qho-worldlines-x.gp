# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$\ii t$'
unset key
L = 5
set xrange [0 : L]
plot 'qho-ensemble.data' using 2 : 1 with linespoints linetype 1 pointtype 7
