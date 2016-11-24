# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$i_P$'
set ylabel '$E$'
# set ytics nomirror
# set y2tics
# set y2label '$X$'
plot 'qho-tde.data' using 1 : 2 with lines linetype 1 title '$E$', \
     'qho-tde.data' using 1 : 3 with lines linetype 2 title '$\langle E\rangle$', \
     'qho-tde.data' using 1 : ($3 - $4) with lines linetype 2 notitle, \
     'qho-tde.data' using 1 : ($3 + $4) with lines linetype 2 notitle
