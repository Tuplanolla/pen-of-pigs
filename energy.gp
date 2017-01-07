set terminal epslatex
set output 'energy.tex'
C_C = `cat 'run-latest/energy-corrtime.data' | cut -d ' ' -f 1`
C_V = `cat 'run-latest/energy-corrtime.data' | cut -d ' ' -f 2`
C_M = `cat 'run-latest/energy-corrtime.data' | cut -d ' ' -f 3`
set xlabel '$i_P$'
set ylabel '$E$'
set autoscale xfix
plot 'run-latest/energy.data' using 1 : 2 \
  with points linetype 1 pointtype 1 title '$E_C$', \
  for [i = -1 : 1 : 1] '' using 1 : ($3 + i * sqrt(C_C) * $4) \
  with lines linetype 1 notitle, \
  '' using 1 : 5 \
  with points linetype 2 pointtype 2 title '$E_V$', \
  for [i = -1 : 1 : 1] '' using 1 : ($6 + i * sqrt(C_V) * $7) \
  with lines linetype 2 notitle, \
  '' using 1 : 8 \
  with points linetype 3 pointtype 3 title '$E_M$', \
  for [i = -1 : 1 : 1] '' using 1 : ($9 + i * sqrt(C_M) * $10) \
  with lines linetype 3 notitle
