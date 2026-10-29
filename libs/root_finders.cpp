#include "root_finders.h"
#include <iostream>

#ifndef ROOTFIND_DEBUG
#define ROOTFIND_DEBUG 0
#endif

int Bisection(double (*func) (double), double a, double b,
              double xtol, double ytol, double &zero) {
  double fmid;      // = func evaluated at midpoint
                    // (zero == midpoint)
  int counter = 0;

  // check if a is a root 
  double fa = func(a);
  if (fabs(fa) <= ytol) {
    zero = a;
    std::cout << "a = " << a << " is a zero: succesful bisection" << std::endl;
    return 0;
  }

  // check if b is a root
  double fb = func(b);
  if (fabs(fb) <= ytol) {
    zero = b;
    std::cout << "b = " << b << " is a zero: succesful bisection" << std::endl;
    return 0;
  }

  // check if [a,b] is a correct bracket
  // (with the update rules it is guaranteed that f(a) * f(b) > 0 always,
  // after this check)
  if (fa * fb > 0) { 
    std::cout << "bracket has f(a) * f(b) > 0: unsuccesful bisection" << std::endl; 
    return 1;
  } 

  
  // while del > xtol, and f(zero) > ytol
  do {
    if (counter > 128) {
      std::cout << "too many iterations: unsuccesful bisection" << std::endl;
      return 1;
    }

    // applying bisection algorithm
    counter++;
    zero = 0.5 * (a + b);
    fmid = func(zero);

    #if ROOTFIND_DEBUG
    std::cout << "Bisection(): k = " << counter
         << "; [a, b] = [" << a << ", " << b << "]; " 
         << "xm = " << zero
         << "; dx = " << fabs(b-a)
         << "; fm = " << fmid << std::endl;
    #endif

    // if midpoint has the same sign of b  update b
    // else                                update a
    if (fa * fmid < 0) {
      b = zero;
      fb = fmid;
    } else {  // if fmid == 0 the ytol will take care of it
      a = zero;
      fa = fmid;
    }
  } while ((fabs(b-a) > xtol) && (fabs(fmid) > ytol)); 

  return 0;
}

int FalsePosition(double (*func) (double), double a, double b,
              double xtol, double ytol, double &zero) {
  double del = fabs(b-a); // = distance between updated points
  double fmid;            // = func evaluated at intermediate point
                          // (zero == intermediate point)
  int counter = 0;

  // check if a is a root
  double fa = func(a);
  if (fa == 0) {
    zero = a;
    std::cout << "a = " << a << " is a zero: succesful false position" 
              << std::endl;
    return 0;
  }

  // check if b is a root
  double fb = func(b);
  if (fb == 0) {
    zero = b;
    std::cout << "b = " << b << " is a zero: succesful false position" 
              << std::endl;
    return 0;
  }

  // check if [a,b] is a correct bracket
  // (with the update rules it is guaranteed that f(a) * f(b) > 0 always,
  // after this check)
  if (fa * fb > 0) { 
    std::cout << "bracket with f(a) * f(b) > 0 : " 
         << "unsuccesful false position" << std::endl;
    return 1;
  } 

  // while del > xtol, and f(zero) > ytol
  do {
    if (counter > 128) {
      std::cout << "too many iterations: unsuccesful false position" 
                << std::endl;
      return 1;
    }

    // applying false position algorithm
    counter++;
    zero = (a * fb - b * fa) / (fb - fa);
    fmid = func(zero);

    #if ROOTFIND_DEBUG
    std::cout << "FalsePosition(): k = " << counter
         << "; [a, b] = [" << a << ", " << b << "]; " 
         << "xm = " << zero
         << "; |del| = " << del
         << "; fm = " << fmid << std::endl;
    #endif

    // if midpoint has the same sign of b  update b
    // else                                update a
    if (fa * fmid < 0) {
      del = fabs(b - zero);
      b   = zero;
      fb  = fmid;
    } else { // if fmid == 0 than ytol will take care of it
      del = fabs(a - zero);
      a   = zero;
      fa  = fmid;
    }
  } while (del > xtol && fabs(fmid) > ytol);

  return 0;
}

int Secant(double (*func) (double), double a, double b, 
           double xtol, double ytol, double &zero) {
  double del    = b - a; // = x_k+1 - x_k (distance update)
  double olddel;         // = x_k - x_k-1 (needed to perform conv checks)
  int counter   = 1;

  // check if a is a root
  double fa = func(a);
  if (fa == 0) {
    zero = a;
    std::cout << "a = " << a << " is a zero: succesful secant" << std::endl;
    return 0;
  }

  // check if b is a root
  double fb = func(b);
  if (fb == 0) {
    zero = b;
    std::cout << "b = " << b << " is a zero: succesful secant" << std::endl;
    return 0;
  }

  while ((fabs(del) > xtol) && (fabs(fb) > ytol)) {
    if (counter > 128) {
      std::cout << "too many iterations: unsuccesful secant" << std::endl;
      return 1;
    }

    if (fabs(fb - fa) < 1e-18) {
      std::cout << "secant line is almost horizontal (div/0): "
           << "unsuccesful secant" << std::endl;
      return 1;
    }

    #if ROOTFIND_DEBUG
    std::cout << "Secant(): k = " << counter
         << "; xa = " << a << ", xb = " << b 
         << "; dx = " << del << std::endl;
    #endif

    // applying secant algorithm
    counter++;
    olddel = del;                       // saving previous del
    del    = fb * (b - a) / (fb - fa);  // compute new del with rule
    a      = b;                         // a    <- x_k
    fa     = fb;                        // f(a) <- f(x_k)
    b     -= del;                       // b    <- x_k+1 = x_k - del
    fb     = func(b);                   // f(b) <- f(x_k+1)

    if (fabs(del) > fabs(olddel) * 5.0) {
      std::cout << "intervals started increasing in size: "
           << "unsuccesful secant" << std::endl;
      return 1;
    }
  }
  // saving last x_k+1 to zero
  zero = b;

  return 0;
}

int Newton(double (*func) (double), double(*derivative) (double),
           double a, double b, 
           double xtol, double ytol, double &zero) {
  int counter = 0;

  // init first guess as midpoint of interval
  zero = 0.5 * (a + b);
  double fxc  = func(zero);
  double dfxc = derivative(zero);

  double del = fxc / derivative(zero); // = x_k+1 - x_k (update distance)
  double olddel;                       // = x_k - x_k-1 (for conv checks)
 
  while ((fabs(del) > xtol) && (fabs(fxc) > ytol)) {
    if (counter > 128) {
      std::cout << "too many iterations: unsuccesful Newton" << std::endl;
      return 1;
    }

      if (fabs(dfxc) < 1e-18) {
      std::cout << "tangent line is almost horizontal (div/0): "
           << "unsuccesful Newton" << std::endl;
      return 1;
    }

    // applying Newton algorithm
    counter++;
    olddel = del;                      // saving last update distance
    del    = fxc / dfxc; 
    zero  -= del;
    fxc    = func(zero);
    dfxc   = derivative(zero);

    if (fabs(del) > fabs(olddel) * 5.0) {
      std::cout << "intervals started increasing in size: "
           << "unsuccesful Newton" << std::endl;
      return 1;
    }

    #if ROOTFIND_DEBUG
    std::cout << "Newton(): k = " << counter
         << "; xc = " << zero 
         << "; dx = " << del << std::endl;
    #endif
  }

  return 0;
}

int Bracket(double (*func) (double), double a, double b,
            int n, double *start, double *end, int &possibleroots) {
  double len = b - a;
  double del = len / n;

  // init of the extremes for the first bracket
  double c = a;
  double d = a + del;

  // init of function values for the first bracket
  double fc = func(c);
  double fd = func(d);

  possibleroots = 0;

  for (int i = 0; i < n + 1; i++) {
    #if ROOTFIND_DEBUG
    std::cout << "(c,d) = (" << c << ", " << d << ");  "
         << "(fc, fd) = (" << fc << ", " << fd  << ");  " << std::endl;
    #endif

    if (fc * fd < 0) {
      start[possibleroots] = c;
      end[possibleroots]   = d;
      possibleroots++;
      #if ROOTFIND_DEBUG
      std::cout << " ^ added!" << std::endl;
      #endif
    }

    // update to extremes of next bracket
    c = a + del * i;
    d = a + del * (i + 1);

    // update function values at extremes of next bracket
    // (only one evaluation is needed!)
    fc = fd;
    fd = func(d);
  }

  return 0;
}
