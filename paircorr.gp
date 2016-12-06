set terminal epslatex
set output 'paircorr.tex'
N = `cat 'run-latest/npoly.data'`
if (N < 2) {exit}
L = `cat 'run-latest/length.data'`
set xlabel '$r$'
set ylabel '$g$'
unset key
set xrange [0 : L]
plot 'run-latest/paircorr.data' using 1 : 2 \
  with lines linetype 1
