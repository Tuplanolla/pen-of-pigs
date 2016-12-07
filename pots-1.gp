set terminal epslatex
set output 'pots-1.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 1) {exit}
L = `cat 'run-latest/length.data'`
set xlabel '$x$'
set ylabel '$V$'
unset key
set xrange [0 - L / 2 : L + L / 2]
set yrange [-50 < * : * < 50]
plot for [i = 0 : 1] \
  for [dx = -L : L : L] \
  'run-latest/pots.data' using \
  ($1 + dx) : i + 2 \
  with lines linetype i + 1
