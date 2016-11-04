# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$y$'
unset key
fx(i, x, y) = i == 0 ? \
  (x2 = NaN, y2 = NaN, x1 = x, y1 = y, x2) : \
  (x2 = x1, y2 = y1, x1 = x, y1 = y, x2)
fy(i, x, y) = y2
dfx(i, x, y) = x1 - x2
dfy(i, x, y) = y1 - y2
plot 'qho-ensemble.data' using \
  2 : 3 every ::0::0 \
  with points linetype 1 pointtype 7, \
  'qho-ensemble.data' using \
  (fx($1, $2, $3)) : (fy($1, $2, $3)) : (dfx($1, $2, $3)) : (dfy($1, $2, $3)) \
  with vectors linetype 1
