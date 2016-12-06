# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$r$'
set ylabel '$g$'
unset key
L = `cat 'qho-length.data'`
set xrange [0 : L]
plot 'qho-paircorrelation.data' using \
  1 : 2 \
  with lines linetype 1
