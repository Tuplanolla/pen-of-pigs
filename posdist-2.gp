set terminal epslatex
set output 'posdist-2.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 2) {exit}
L = `cat 'run-latest/length.data'`
K = `cat 'run-latest/nsubdiv.data'`
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$p$'
unset key
set dgrid3d K, K
set hidden3d
set xyplane relative 0.1
set xrange [0 : L]
set yrange [0 : L]
splot for [dx = -L : L : L] for [dy = -L : L : L] \
  'run-latest/posdist.data' using \
  ($1 + dx) : ($1 + dy) : 3 \
  with lines linetype 1
