set terminal epslatex
set output 'polys-3.tex'
d = `cat 'run-latest/ndim.data'`
if (d < 3) {exit}
periodic = `cat 'run-latest/periodic.data'`
L = `cat 'run-latest/length.data'`
set xlabel '$x$'
set ylabel '$y$'
set zlabel '$z$'
unset key
set xyplane relative 0.1
p(q) = q > L / 2.0 ? q - L : q < -L / 2.0 ? q + L : q
fx(i, x, y, z) = i == 0 ? \
  (x2 = NaN, y2 = NaN, z2 = NaN, x1 = x, y1 = y, z1 = z, x2) : \
  (x2 = x1, y2 = y1, z2 = z1, x1 = x, y1 = y, z1 = z, x2)
fy(i, x, y, z) = y2
fz(i, x, y, z) = z2
if (periodic) {
  dfx(i, x, y, z) = p(x1 - x2)
  dfy(i, x, y, z) = p(y1 - y2)
  dfz(i, x, y, z) = p(z1 - z2)
  set xrange [-L / 2.0 : L + L / 2.0]
  set yrange [-L / 2.0 : L + L / 2.0]
  set zrange [-L / 2.0 : L + L / 2.0]
  set object 1 polygon from 0, 0, 0 to L, 0, 0 to L, L, 0 to 0, L, 0 to 0, 0, 0
  set object 2 polygon from 0, 0, L to L, 0, L to L, L, L to 0, L, L to 0, 0, L
  set object 3 polygon from 0, 0, 0 to L, 0, 0 to L, 0, L to 0, 0, L to 0, 0, 0
  set object 4 polygon from 0, L, 0 to L, L, 0 to L, L, L to 0, L, L to 0, L, 0
  set object 5 polygon from 0, 0, 0 to 0, 0, L to 0, L, L to 0, L, 0 to 0, 0, 0
  set object 6 polygon from L, 0, 0 to L, 0, L to L, L, L to L, L, 0 to L, 0, 0
  splot for [dx = -L : L : L] for [dy = -L : L : L] for [dz = -L : L : L] \
    'run-latest/polys.data' using \
    ($2 + dx) : ($3 + dy) : ($4 + dz) every ::0::0 \
    with points linetype 1 pointtype 7, \
    for [dx = -L : L : L] for [dy = -L : L : L] for [dz = -L : L : L] \
    '' using \
    (fx($1, $2 + dx, $3 + dy, $4 + dz)) : \
    (fy($1, $2 + dx, $3 + dy, $4 + dz)) : \
    (fz($1, $2 + dx, $3 + dy, $4 + dz)) : \
    (dfx($1, $2 + dx, $3 + dy, $4 + dz)) : \
    (dfy($1, $2 + dx, $3 + dy, $4 + dz)) : \
    (dfz($1, $2 + dx, $3 + dy, $4 + dz)) \
    with vectors linetype 1 filled
} else {
  dfx(i, x, y, z) = x1 - x2
  dfy(i, x, y, z) = y1 - y2
  dfz(i, x, y, z) = z1 - z2
  set xrange [-L / 2.0 : L / 2.0]
  set yrange [-L / 2.0 : L / 2.0]
  set zrange [-L / 2.0 : L / 2.0]
  unset object 1
  unset object 2
  unset object 3
  unset object 4
  unset object 5
  unset object 6
  splot 'run-latest/polys.data' using 2 : 3 : 4 every ::0::0 \
    with points linetype 1 pointtype 7, \
    '' using \
    (fx($1, $2, $3, $4)) : (fy($1, $2, $3, $4)) : (fz($1, $2, $3, $4)) : \
    (dfx($1, $2, $3, $4)) : (dfy($1, $2, $3, $4)) : (dfz($1, $2, $3, $4)) \
    with vectors linetype 1 filled
}
