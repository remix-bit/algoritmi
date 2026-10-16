reset

set parametric
set trange [0:2*pi]

fx(t) = cos(t)
fy(t) = sin(t)
plot fx(t), fy(t)
unset parametric
