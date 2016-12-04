# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$P$'
unset key
N = `cat 'qho-subdivisions-2d.data'`
L = `cat 'qho-length-2d.data'`
set xrange [0 : L]
set yrange [0 : L]
set dgrid3d N, N
set hidden3d
splot 'qho-probability-2d.data' using \
  1 : 2 : 3 \
  with lines linetype 1
