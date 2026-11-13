reset

set term qt 1 size 500,400
set title "phi for various initial guesses"
set xlabel "r"
set ylabel "\phi"
set yrange [-10:15]
set grid

plot "poisson.dat" index 0 using 1:2 with lines title "s = 0"
replot "poisson.dat" index 1 using 1:2 with lines title "s = 0.2"
replot "poisson.dat" index 2 using 1:2 with lines title "s = 0.4"
replot "poisson.dat" index 3 using 1:2 with lines title "s = 0.6"
replot "poisson.dat" index 4 using 1:2 with lines title "s = 0.8"
replot "poisson.dat" index 5 using 1:2 with lines title "s = 1.0"

set term qt 2 size 500,400
set title "residual against s"
set xlabel "s"
set ylabel "res"
set yrange [-10:90]

plot "poisson.dat" index 6 using 1:2 with lines

set term qt 3 size 500,400
set title "final result"
set xlabel "r"
set ylabel "phi"
set yrange [-5:3]

plot "poisson.dat" index 7 using 1:2 title "numerical result"
replot 1 - 0.5 * (x + 2) * exp(-x) title "analytical result"
