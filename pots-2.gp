set terminal epslatex
set output 'pots-2.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 2) {exit}
periodic = `cat 'run-latest/periodic.data'`
L = `cat 'run-latest/length.data'`
K = `cat 'run-latest/nsubdiv.data'`
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$V$'
unset key
set hidden3d
set xyplane relative 0.1
set dgrid3d K, K
set zrange [-50 < * : * < 50]
if (periodic) {
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  splot for [i = 0 : 1 : 1] \
    for [dx = -L : L : L] for [dy = -L : L : L] \
    'run-latest/pots.data' using \
    ($1 + dx) : ($2 + dy) : i + 3 \
    with lines linetype 2 * i + 1
} else {
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  splot for [i = 0 : 1 : 1] \
    'run-latest/pots.data' using 1 : 2 : i + 3 \
    with lines linetype 2 * i + 1
}
