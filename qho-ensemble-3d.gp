# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$z$'
unset key
fx(i, x, y, z) = i == 0 ? \
  (x2 = NaN, y2 = NaN, z2 = NaN, x1 = x, y1 = y, z1 = z, x2) : \
  (x2 = x1, y2 = y1, z2 = z1, x1 = x, y1 = y, z1 = z, x2)
fy(i, x, y, z) = y2
fz(i, x, y, z) = z2
dfx(i, x, y, z) = x1 - x2
dfy(i, x, y, z) = y1 - y2
dfz(i, x, y, z) = z1 - z2
splot 'qho-ensemble.data' using \
  2 : 3 : 4 every ::0::0 \
  with points linetype 1 pointtype 7, \
  'qho-ensemble.data' using \
  (fx($1, $2, $3, $4)) : (fy($1, $2, $3, $4)) : (fz($1, $2, $3, $4)) : \
  (dfx($1, $2, $3, $4)) : (dfy($1, $2, $3, $4)) : (dfz($1, $2, $3, $4)) \
  with vectors linetype 1
