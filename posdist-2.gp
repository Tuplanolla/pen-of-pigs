set terminal epslatex
set output 'posdist-2.tex'
stats 'run-latest/ndim.data' using 1 name 'stats_d'
d = stats_d_mean
if (d != 2) {exit}
stats 'run-latest/periodic.data' using 1 name 'stats_periodic'
periodic = stats_periodic_mean
stats 'run-latest/length.data' using 1 name 'stats_L'
L = stats_L_mean
stats 'run-latest/nsubdiv.data' using 1 name 'stats_K'
K = stats_K_mean
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$p$'
unset key
set xyplane relative 0.1
set hidden3d
set dgrid3d K, K
if (periodic) {
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  splot for [dx = -L : L : L] for [dy = -L : L : L] \
    'run-latest/posdist.data' using ($1 + dx) : ($1 + dy) : 3 \
    with lines linetype 1
} else {
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  splot 'run-latest/posdist.data' using 1 : 2 : 3 \
    with lines linetype 1
}
