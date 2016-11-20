# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$i$'
set ylabel '$d$'
unset key
plot 'qho-drift.data' using 1 : 2 with linespoints linetype 1 pointtype 7
