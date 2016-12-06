set terminal epslatex
set output 'pots-1.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 1) {exit}
L = `cat 'run-latest/length.data'`
set xlabel '$x$'
set ylabel '$V$'
unset key
set xrange [0 - L / 2 : L + L / 2]
plot for [i = 2 : 3] \
  for [dx = -L : L : L] \
  'run-latest/pots.data' using \
  ($1 + dx) : i \
  with lines linetype 1
