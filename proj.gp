set terminal epslatex
set output 'proj.tex'
stats 'run-latest/periodic.data' using 1 name 'stats_periodic'
periodic = stats_periodic_mean
stats 'run-latest/length.data' using 1 name 'stats_L'
L = stats_L_mean
set xlabel '$x$'
set ylabel '$\tau$'
unset key
p(q) = q > L / 2.0 ? q - L : q < -L / 2.0 ? q + L : q
fx(i, x, y) = i == 0 ? \
  (x2 = NaN, y2 = NaN, x1 = x, y1 = y, x2) : \
  (x2 = x1, y2 = y1, x1 = x, y1 = y, x2)
fy(i, x, y) = y2
dfy(i, x, y) = y1 - y2
set autoscale yfix
if (periodic) {
  dfx(i, x, y) = p(x1 - x2)
  set xrange [0.0 : L]
  plot for [dx = -L : L : L] \
    'run-latest/polys.data' using \
    (fx($1, $2 + dx, $1)) : (fy($1, $2 + dx, $1)) : \
    (dfx($1, $2 + dx, $1)) : (dfy($1, $2 + dx, $1)) \
    with vectors linetype 1 filled
} else {
  dfx(i, x, y) = x1 - x2
  set xrange [-L / 2.0 : L / 2.0]
  plot 'run-latest/polys.data' using \
    (fx($1, $2, $1)) : (fy($1, $2, $1)) : \
    (dfx($1, $2, $1)) : (dfy($1, $2, $1)) \
    with vectors linetype 1 filled
}
