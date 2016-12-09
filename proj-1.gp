set terminal epslatex
set output 'proj-1.tex'
d = `cat 'run-latest/ndim.data'`
if (d != 1) {exit}
L = `cat 'run-latest/length.data'`
M = `cat 'run-latest/nbead.data'`
set xlabel '$x$'
set ylabel '$\tau$'
unset key
p(q) = q > L / 2.0 ? q - L : q < -L / 2.0 ? q + L : q
fx(i, x, y) = i == 0 ? \
  (x2 = NaN, y2 = NaN, x1 = x, y1 = y, x2) : \
  (x2 = x1, y2 = y1, x1 = x, y1 = y, x2)
fy(i, x, y) = y2
dfx(i, x, y) = p(x1 - x2)
dfy(i, x, y) = p(y1 - y2)
set xrange [0.0 : L]
set yrange [0 : M]
plot for [dx = -L : L : L] \
  'run-latest/polys.data' using \
  (fx($1, $2 + dx, $1)) : (fy($1, $2 + dx, $1)) : \
  (dfx($1, $2 + dx, $1)) : (dfy($1, $2 + dx, $1)) \
  with vectors linetype 1 filled
