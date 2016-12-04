# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$V$'
unset key
L = `cat 'qho-length-1d.data'`
set xrange [0 - L / 2 : L + L / 2]
plot for [i = 2 : 3] \
  for [dx = -L : L : L] \
  'qho-potential-1d.data' using \
  ($1 + dx) : i \
  with lines linetype 1
