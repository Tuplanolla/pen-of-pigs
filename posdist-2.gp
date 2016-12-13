set terminal epslatex
set output 'posdist-2.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 2) {exit}
periodic = `cat 'run-latest/periodic.data'`
L = `cat 'run-latest/length.data'`
K = `cat 'run-latest/nsubdiv.data'`
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$p$'
unset key
set hidden3d
set xyplane relative 0.1
set dgrid3d K, K
if (periodic) {
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  splot for [dx = -L : L : L] for [dy = -L : L : L] \
    'run-latest/posdist.data' using \
    ($1 + dx) : ($1 + dy) : 3 \
    with lines linetype 1
} else {
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  splot 'run-latest/posdist.data' using 1 : 2 : 3 \
    with lines linetype 1
}
