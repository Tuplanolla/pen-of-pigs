set terminal epslatex
set output 'params.tex'
set xlabel '$i_T + i_P$'
set ylabel '$h$'
set log y
plot 'run-latest/params.data' using ($1 + $2) : 5 \
  with lines title 'SSM', \
  '' using ($1 + $2) : 8 \
  with lines title 'CMD'
