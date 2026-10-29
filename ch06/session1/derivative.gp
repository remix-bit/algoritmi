reset

set log x
set log y

plot "derivative.dat"
replot "derivative.dat" using 1:3
replot "derivative.dat" using 1:4
replot "derivative.dat" using 1:5
