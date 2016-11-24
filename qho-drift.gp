# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$i_T + i_P$'
set ylabel '$\delta$'
set log y
plot 'qho-drift.data' using 1 : 2 with lines title 'SS', \
     'qho-drift.data' using 1 : 3 with lines title 'CoMD'
