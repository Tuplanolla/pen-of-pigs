# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$V$'
unset key
N = `cat 'qho-subdivisions-2d.data'` + 1
L = `cat 'qho-length-2d.data'`
set xrange [0 - L / 2 : L + L / 2]
set yrange [0 - L / 2 : L + L / 2]
set dgrid3d N, N
set hidden3d
splot for [i = 3 : 4] \
  for [dx = -L : L : L] for [dy = -L : L : L] \
  'qho-potential-2d.data' using \
  ($1 + dx) : ($2 + dy) : i \
  with lines linetype 1
