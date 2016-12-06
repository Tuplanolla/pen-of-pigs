# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$p$'
unset key
N = `cat 'qho-subdivisions-2d.data'`
L = `cat 'qho-length-2d.data'`
set xrange [0 : L]
set yrange [0 : L]
set dgrid3d N, N
set hidden3d
splot for [dx = -L : L : L] for [dy = -L : L : L] \
  'qho-probability-2d.data' using \
  ($1 + dx) : ($1 + dy) : 3 \
  with lines linetype 1
