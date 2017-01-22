set terminal epslatex
set output 'energy.tex'
stats 'run-latest/corrtime.data' using 1 name 'stats_C_T'
C_T = stats_C_T_mean
stats 'run-latest/corrtime.data' using 2 name 'stats_C_M'
C_M = stats_C_M_mean
stats 'run-latest/corrtime.data' using 3 name 'stats_C_V'
C_V = stats_C_V_mean
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
