int Bisection(double (*)(double), double, double, double, double, double &);

int FalsePosition(double (*)(double), double, double, 
                  double, double, double &);

int Secant(double (*)(double), double, double,
           double, double, double &);

int Newton(double (*)(double), double (*)(double), 
           double, double, double, double, double &);

int Bracket(double (*)(double), double, double,
            int, double*, double*, int&);
