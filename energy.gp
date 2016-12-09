set terminal epslatex
set output 'energy.tex'
C = `cat 'run-latest/energy-corrtime.data'`
set xlabel '$i_P$'
set ylabel '$E$'
plot 'run-latest/energy.data' using 1 : 2 \
  with lines linetype 1 title '$E$', \
  'run-latest/energy.data' using 1 : 3 \
  with lines linetype 2 title '$\langle E\rangle$', \
  'run-latest/energy.data' using \
  1 : ($3 - sqrt(C) * $4) \
  with lines linetype 2 notitle, \
  'run-latest/energy.data' using \
  1 : ($3 + sqrt(C) * $4) \
  with lines linetype 2 notitle
