set terminal epslatex
set output 'pots-2.tex'
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
set zlabel '$V$'
unset key
set xyplane relative 0.1
set hidden3d
set dgrid3d K, K
set zrange [-50 < * : * < 50]
if (periodic) {
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  splot for [i = 0 : 1 : 1] \
    for [dx = -L : L : L] for [dy = -L : L : L] \
    'run-latest/pots.data' using ($1 + dx) : ($2 + dy) : i + 3 \
    with lines linetype 2 * i + 1
} else {
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  splot for [i = 0 : 1 : 1] \
    'run-latest/pots.data' using 1 : 2 : i + 3 \
    with lines linetype 2 * i + 1
}
