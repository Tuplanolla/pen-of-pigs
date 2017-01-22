set terminal epslatex
set output 'pots-1.tex'
stats 'run-latest/ndim.data' using 1 name 'stats_d'
d = stats_d_mean
if (d != 1) {exit}
stats 'run-latest/periodic.data' using 1 name 'stats_periodic'
periodic = stats_periodic_mean
stats 'run-latest/length.data' using 1 name 'stats_L'
L = stats_L_mean
set xlabel '$x$'
set ylabel '$V$'
unset key
set yrange [-50 < * : * < 50]
if (periodic) {
  set xrange [0.0 - L / 2.0 : L + L / 2.0]
  plot for [i = 0 : 1 : 1] \
    for [dx = -L : L : L] \
    'run-latest/pots.data' using ($1 + dx) : i + 2 \
    with lines linetype i + 1
} else {
  set xrange [0.0 - L / 2.0 : 0.0 + L / 2.0]
  plot for [i = 0 : 1 : 1] \
    'run-latest/pots.data' using 1 : i + 2 \
    with lines linetype i + 1
}
