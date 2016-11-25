# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$V$'
unset key
L = `cat 'qho-length.data'`
set xrange [0 - L / 2 : L + L / 2]
plot for [dx = -L : L : L] \
  'qho-potential.data' using \
  ($1 + dx) : 2 \
  with lines linetype 1
