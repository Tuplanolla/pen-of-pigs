set terminal epslatex
set output 'raddist.tex'
stats 'run-latest/npoly.data' using 1 name 'stats_N'
N = stats_N_mean
if (N < 2) {exit}
stats 'run-latest/length.data' using 1 name 'stats_L'
L = stats_L_mean
set xlabel '$r$'
set ylabel '$g$'
unset key
set xrange [0.0 : L / 2.0]
plot 'run-latest/raddist.data' using 1 : 2 \
  with lines linetype 1
