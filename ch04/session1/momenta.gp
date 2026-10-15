set term qt size 600,400
set log x
set log y
set xlabel "N"  font ",14"
set ylabel "err" font ",14"

plot "prn_momenta.dat" using 1:2 title "k=1"
replot "prn_momenta.dat" using 1:3 title "k=2"
replot 1.0/sqrt(x) title "1/sqrt(N)"
