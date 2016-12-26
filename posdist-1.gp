set terminal epslatex
set output 'posdist-1.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 1) {exit}
periodic = `cat 'run-latest/periodic.data'`
L = `cat 'run-latest/length.data'`
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
