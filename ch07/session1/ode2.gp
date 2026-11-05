reset

set xrange [-4:4]
set yrange [-4:4]


set xlabel "x"
set ylabel "y"

plot "ode2.dat" using 2:3 title "Euler" w lp
replot "ode2.dat" using 4:5 title "RK2" w lp
replot "ode2.dat" using 6:7 title "RK4" w lp
