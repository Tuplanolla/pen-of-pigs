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
set xrange [0 - L / 2 : L + L / 2]
set yrange [0 - L / 2 : L + L / 2]
set dgrid3d K, K
set hidden3d
splot for [i = 3 : 4] \
  for [dx = -L : L : L] for [dy = -L : L : L] \
  'run-latest/pots.data' using \
  ($1 + dx) : ($2 + dy) : i \
  with lines linetype 1
