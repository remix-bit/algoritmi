reset

# we should also increase font size...
set xlabel "t"
set ylabel "N(t)"
set title "Radioactive Decay"

set log y
plot "decay.dat"
replot 1000 * exp(-0.01 * x) with lines
