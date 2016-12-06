# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$i_T + i_P$'
set ylabel '$\delta$'
set log y
plot 'qho-drift.data' using ($1 + $2) : 5 with lines title 'SSM', \
     'qho-drift.data' using ($1 + $2) : 8 with lines title 'CMD'
