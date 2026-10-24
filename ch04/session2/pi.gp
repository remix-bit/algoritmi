reset

set log x

plot "pi.dat"
replot "pi.dat" using 1:3 with lines
replot 1 / sqrt(x)
