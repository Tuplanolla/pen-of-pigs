# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$i$'
set ylabel '$\delta$'
set log y
plot 'qho-drift.data' using 1 : 2 title 'SS' \
     with linespoints linetype 1 pointtype 7, \
     'qho-drift.data' using 1 : 3 title 'CoMD' \
     with linespoints linetype 2 pointtype 7
