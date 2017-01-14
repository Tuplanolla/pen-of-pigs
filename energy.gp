set terminal epslatex
set output 'energy.tex'
C_T = `cat 'run-latest/corrtime.data' | cut -d ' ' -f 1`
C_M = `cat 'run-latest/corrtime.data' | cut -d ' ' -f 2`
C_V = `cat 'run-latest/corrtime.data' | cut -d ' ' -f 3`
set xlabel '$i$'
set ylabel '$E / N$'
set autoscale xfix
set yrange [-50 < * : * < 50]
plot 'run-latest/energy.data' using 1 : 2 \
  with points linetype 1 pointtype 1 title 'T', \
  for [i = -1 : 1 : 1] '' using 1 : ($3 + i * sqrt(C_T) * $4) \
  with lines linetype 1 notitle, \
  '' using 1 : 5 \
  with points linetype 2 pointtype 2 title 'M', \
  for [i = -1 : 1 : 1] '' using 1 : ($6 + i * sqrt(C_M) * $7) \
  with lines linetype 2 notitle, \
  '' using 1 : 8 \
  with points linetype 3 pointtype 3 title 'V', \
  for [i = -1 : 1 : 1] '' using 1 : ($9 + i * sqrt(C_V) * $10) \
  with lines linetype 3 notitle
