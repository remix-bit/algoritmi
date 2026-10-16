reset

# we should also increase font size...
set xlabel "x"
set ylabel "f(x)"
set title "Gaussian Estimate through Neumann Method"

plot "neumann_gaussian.dat"
replot exp(-0.5 * (x / 0.5)**2) / (0.5 * sqrt(2 * pi)) with lines
