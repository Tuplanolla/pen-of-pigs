# set terminal epslatex
# set output 'epslatex.tex'
set terminal wxt noraise
set xlabel '$x$'
set ylabel '$y$'
unset key
L = `cat 'qho-length.data'`
p(q) = q > L / 2 ? q - L : q < -L / 2 ? q + L : q
fx(i, x, y) = i == 0 ? \
  (x2 = NaN, y2 = NaN, x1 = x, y1 = y, x2) : \
  (x2 = x1, y2 = y1, x1 = x, y1 = y, x2)
fy(i, x, y) = y2
dfx(i, x, y) = p(x1 - x2)
dfy(i, x, y) = p(y1 - y2)
set xrange [0 - L / 2 : L + L / 2]
set yrange [0 - L / 2 : L + L / 2]
set object 1 rectangle from 0, 0 to L, L fillstyle empty
plot for [dx = -L : L : L] for [dy = -L : L : L] \
  'qho-ensemble.data' using \
  ($2 + dx) : ($3 + dy) every ::0::0 \
  with points linetype 1 pointtype 7, \
  for [dx = -L : L : L] for [dy = -L : L : L] \
  'qho-ensemble.data' using \
  (fx($1, $2 + dx, $3 + dy)) : (fy($1, $2 + dx, $3 + dy)) : \
  (dfx($1, $2 + dx, $3 + dy)) : (dfy($1, $2 + dx, $3 + dy)) \
  with vectors linetype 1
pause 10
reread
