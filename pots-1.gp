set terminal epslatex
set output 'pots-1.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 1) {exit}
periodic = `cat 'run-latest/periodic.data'`
L = `cat 'run-latest/length.data'`
set xlabel '$x$'
set ylabel '$V$'
unset key
set yrange [-50 < * : * < 50]
if (periodic) {
  set xrange [0.0 - L / 2.0 : L + L / 2.0]
  plot for [i = 0 : 1 : 1] \
    for [dx = -L : L : L] \
    'run-latest/pots.data' using \
    ($1 + dx) : i + 2 \
    with lines linetype i + 1
} else {
  set xrange [0.0 - L / 2.0 : 0.0 + L / 2.0]
  plot for [i = 0 : 1 : 1] \
    'run-latest/pots.data' using 1 : i + 2 \
    with lines linetype i + 1
}
