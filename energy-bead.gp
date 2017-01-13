set terminal epslatex
set output 'energy-bead.tex'
C_T = `cat 'run-latest/energy-corrtime.data' | cut -d ' ' -f 2`
C_M = `cat 'run-latest/energy-corrtime.data' | cut -d ' ' -f 3`
set xlabel '$\tau$'
set ylabel '$E / N$'
set autoscale xfix
plot 'run-latest/energy-bead.data' using 1 : 2 \
  with lines linetype 1 title 'T', \
  for [i = -1 : 1 : 2] '' using 1 : ($2 + i * sqrt(C_T) * $3) \
  with lines linetype 1 notitle, \
  '' using 1 : 4 \
  with lines linetype 2 title 'M', \
  for [i = -1 : 1 : 2] '' using 1 : ($4 + i * sqrt(C_M) * $5) \
  with lines linetype 2 notitle
