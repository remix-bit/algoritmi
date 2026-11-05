reset

plot "convstudy.dat" using 1:2 title "Euler"
replot "convstudy.dat" using 1:3 title "RK2"
replot "convstudy.dat" using 1:4 title "RK4"
