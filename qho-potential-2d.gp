# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$V$'
unset key
N = `cat 'qho-subdivisions.data'` + 1
L = `cat 'qho-length.data'`
set xrange [0 - L / 2 : L + L / 2]
set yrange [0 - L / 2 : L + L / 2]
set dgrid3d N, N
set hidden3d
splot for [dx = -L : L : L] for [dy = -L : L : L] \
  'qho-potential-2d.data' using \
  ($1 + dx) : ($2 + dy) : 3 \
  with lines linetype 1
