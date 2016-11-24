# set terminal epslatex
# set output 'epslatex.tex'
set xlabel '$i$'
set ylabel '$?$'
set ytics nomirror
set y2tics
set y2label '$\tilde ?$'
plot 'qho-tde.data' using 1 : 2 axes x1y1 with lines title '$E$', \
     'qho-tde.data' using 1 : 3 axes x1y1 with lines title '$\langle E\rangle$', \
     'qho-tde.data' using 1 : 4 axes x1y2 with lines title '$\tilde E$'
