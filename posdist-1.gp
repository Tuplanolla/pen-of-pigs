set terminal epslatex
set output 'posdist-1.tex'
stats 'run-latest/ndim.data' using 1 name 'stats_d'
d = stats_d_mean
if (d != 1) {exit}
stats 'run-latest/periodic.data' using 1 name 'stats_periodic'
periodic = stats_periodic_mean
stats 'run-latest/length.data' using 1 name 'stats_L'
L = stats_L_mean
set xlabel '$x$'
set ylabel '$p$'
unset key
if (periodic) {
  set xrange [-L / 2.0 : L + L / 2.0]
  plot for [dx = -L : L : L] \
    'run-latest/posdist.data' using ($1 + dx) : 2 \
    with lines linetype 1
} else {
  set xrange [-L / 2.0 : L / 2.0]
  plot 'run-latest/posdist.data' using 1 : 2 \
    with lines linetype 1
}
