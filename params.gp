set terminal epslatex
set output 'params.tex'
set xlabel '$i$'
set ylabel '$h$'
set log y
set autoscale xfix
plot 'run-latest/params.data' using ($1 + $2) : 5 with lines title 'SSM', \
  '' using ($1 + $2) : 8 with lines title 'CMD'
