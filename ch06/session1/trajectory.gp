reset

plot "trajectory.dat" with lines title "position"
replot "trajectory.dat" using 1:3 with lines title "velocity"
replot "trajectory.dat" using 1:4 with lines title "acceleration"
