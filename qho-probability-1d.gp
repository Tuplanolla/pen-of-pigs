# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$P$'
unset key
L = `cat 'qho-length-1d.data'`
set xrange [0 : L]
plot 'qho-probability-1d.data' using \
  1 : 2 \
  with lines linetype 1
