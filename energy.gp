set terminal epslatex
set output 'energy.tex'
C = `cat 'run-latest/energy-corrtime.data'`
set xlabel '$i_P$'
set ylabel '$E$'
plot 'run-latest/energy.data' using 1 : 2 \
  with lines linetype 1 title '$E_T$', \
  'run-latest/energy.data' using 1 : 5 \
  with lines linetype 2 title '$E_L$', \
  'run-latest/energy.data' using 1 : 8 \
  with lines linetype 3 title '$E_M$', \
  for [i = -1 : 1 : 1] \
  'run-latest/energy.data' using \
  1 : ($3 + i * sqrt(C) * $4) \
  with lines linetype 1 notitle, \
  for [i = -1 : 1 : 1] \
  'run-latest/energy.data' using \
  1 : ($6 + i * sqrt(C) * $7) \
  with lines linetype 2 notitle, \
  for [i = -1 : 1 : 1] \
  'run-latest/energy.data' using \
  1 : ($9 + i * sqrt(C) * $10) \
  with lines linetype 3 notitle
