set terminal epslatex
set output 'polys-2.tex'
stats 'run-latest/ndim.data' using 1 name 'stats_d'
d = stats_d_mean
if (d < 2) {exit}
stats 'run-latest/periodic.data' using 1 name 'stats_periodic'
periodic = stats_periodic_mean
stats 'run-latest/length.data' using 1 name 'stats_L'
L = stats_L_mean
set xlabel '$x$'
set ylabel '$y$'
unset key
p(q) = q > L / 2.0 ? q - L : q < -L / 2.0 ? q + L : q
fx(i, x, y) = i == 0 ? \
  (x2 = NaN, y2 = NaN, x1 = x, y1 = y, x2) : \
  (x2 = x1, y2 = y1, x1 = x, y1 = y, x2)
fy(i, x, y) = y2
if (periodic) {
  dfx(i, x, y) = p(x1 - x2)
  dfy(i, x, y) = p(y1 - y2)
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  set object 1 rectangle from 0, 0 to L, L fillstyle empty
  plot for [dx = -L : L : L] for [dy = -L : L : L] \
    'run-latest/polys.data' using \
    ($2 + dx) : ($3 + dy) every ::0::0 \
    with points linetype 1 pointtype 7, \
    for [dx = -L : L : L] for [dy = -L : L : L] \
    '' using \
    (fx($1, $2 + dx, $3 + dy)) : (fy($1, $2 + dx, $3 + dy)) : \
    (dfx($1, $2 + dx, $3 + dy)) : (dfy($1, $2 + dx, $3 + dy)) \
    with vectors linetype 1 filled
} else {
  dfx(i, x, y) = x1 - x2
  dfy(i, x, y) = y1 - y2
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  unset object 1
  plot 'run-latest/polys.data' using 2 : 3 every ::0::0 \
    with points linetype 1 pointtype 7, \
    '' using \
    (fx($1, $2, $3)) : (fy($1, $2, $3)) : \
    (dfx($1, $2, $3)) : (dfy($1, $2, $3)) \
    with vectors linetype 1 filled
}
