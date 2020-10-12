#ifndef LEGENDRE
#define LEGENDRE

void legendre_pt_val_and_deriv(int N, double x, double& f, double& df, int& ierr);
void findroot_legendre_newton(int N, double xi, double& x, int& ierr);

#endif