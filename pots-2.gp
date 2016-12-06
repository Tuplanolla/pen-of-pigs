set terminal epslatex
set output 'pots-2.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 2) {exit}
L = `cat 'run-latest/length.data'`
K = `cat 'run-latest/nsubdiv.data'` + 1
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$V$'
unset key
set dgrid3d K, K
set hidden3d
set xyplane relative 0.1
set xrange [0 - L / 2 : L + L / 2]
set yrange [0 - L / 2 : L + L / 2]
splot for [i = 0 : 1] \
  for [dx = -L : L : L] for [dy = -L : L : L] \
  'run-latest/pots.data' using \
  ($1 + dx) : ($2 + dy) : i + 3 \
  with lines linetype 2 * i + 1
