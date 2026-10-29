reset

set xrange [-5:80]

plot sqrt(x) + asin(sqrt(x/250)) + asin(sqrt(x/80)) - pi
replot sqrt(x) + asin(sqrt(x/250)) + asin(sqrt(x/80)) - 2*pi
replot sqrt(x) + asin(sqrt(x/250)) + asin(sqrt(x/80)) - 3*pi

replot 0
