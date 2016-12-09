set terminal epslatex
set output 'posdist-1.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 1) {exit}
L = `cat 'run-latest/length.data'`
set xlabel '$x$'
set ylabel '$p$'
unset key
set xrange [0.0 - L / 2.0 : L + L / 2.0]
plot for [dx = -L : L : L] \
  'run-latest/posdist.data' using \
  ($1 + dx) : 2 \
  with lines linetype 1
