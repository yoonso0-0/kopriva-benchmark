#ifndef LEGENDRE
#define LEGENDRE

template <typename T, size_t N>
void legendre_gauss_nodes_weights(T (&nodes)[N], T (&weights)[N])
{
    std::cout << "template function !!" << std::endl;
    return;
    
    // Chebyshev Gauss points (Kopriva p.63)
    // ...

}

void legendre_findroot_newton(int N, double xi, double& x, int& ierr);
void legendre_pt_val_and_deriv(int N, double x, double& f, double& df, int& ierr);

#endif