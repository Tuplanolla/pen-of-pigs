set terminal epslatex
set output 'energy-bead.tex'
C = `cat 'run-latest/energy-corrtime.data'`
set xlabel '$\tau$'
set ylabel '$E$'
plot 'run-latest/energy-bead.data' using 1 : 2 \
  with lines linetype 1 title '$E_L$', \
  for [i = -1 : 1 : 2] '' using 1 : ($2 + i * sqrt(C) * $3) \
  with lines linetype 1 notitle, \
  # '' using 1 : 5 \
  with lines linetype 2 title '$E_M$', \
  for [i = -1 : 1 : 2] '' using 1 : ($4 + i * sqrt(C) * $5) \
  with lines linetype 2 notitle
