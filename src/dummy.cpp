#include "QDMEwald.hpp"

#include <vector>

int main() {
  Kokkos::initialize();
  {
    // testing the initialization with different data types
    QDMEwald<double> qdme_d(1.0, 1.0, 1.0, 1.0);
    QDMEwald<float> qdme_f(1.0f, 1.0f, 1.0f, 1.0f);
    QDMEwald<long double> qdme_ld(1.0l, 1.0l, 1.0l, 1.0l);

    // testing system size
    double l = 1.00;

    // testing particle positions
    std::vector<double> xyz{0.00, 0.00, 0.00, 0.50, 0.00, 0.00, 0.00, 0.50,
                            0.00, 0.50, 0.50, 0.00, 0.00, 0.00, 0.50, 0.50,
                            0.00, 0.50, 0.00, 0.50, 0.50, 0.50, 0.50, 0.50};

    // testing charges
    std::vector<double> q{-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0};

    // testing dipole momenta
    std::vector<double> d{-0.5, 0.3,  -0.1, -0.4, 0.4,  -0.2, -0.3, 0.5,
                          -0.3, -0.2, 0.4,  -0.4, -0.1, 0.3,  -0.5, -0.0,
                          0.2,  -0.4, 0.1,  0.1,  -0.3, 0.2,  0.0,  -0.2};

    std::vector<double> Q{
        -0.5, 0.3,  -0.1, -0.1, 0.5,  -0.3, 0.1,  0.1,  -0.3, -0.4, 0.4,  -0.2,
        0.0,  0.4,  -0.4, 0.2,  0.0,  -0.2, -0.3, 0.5,  -0.3, 0.1,  0.3,  -0.5,
        0.3,  -0.1, -0.1, -0.2, 0.4,  -0.4, 0.2,  0.2,  -0.4, 0.4,  -0.2, -0.0,
        -0.1, 0.3,  -0.5, 0.3,  0.1,  -0.3, 0.5,  -0.3, 0.1,  0.0,  0.2,  -0.4,
        0.4,  0.0,  -0.2, 0.4,  -0.4, 0.2,  0.1,  0.1,  -0.3, 0.5,  -0.1, -0.1,
        0.3,  -0.5, 0.3,  0.2,  0.0,  -0.2, 0.4,  -0.2, 0.0,  0.2,  -0.4, 0.4};

    qdme_d.compute(xyz, q, d, Q);
  }
  Kokkos::finalize();
}
