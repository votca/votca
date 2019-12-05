#include "QDMEwald.hpp"

#include <iostream>
#include <vector>

int main(int argc, char** argv) {
  Kokkos::initialize();
  {

    // testing system size
    const int cry_l = (argc >= 2) ? atoi(argv[1]) : 16;
    const int N = cry_l * cry_l * cry_l;
    const double l = (double)cry_l / 2.0;

    // testing particle positions
    std::vector<double> xyz(3 * N);

    // testing charges
    std::vector<double> q(N);

    // testing dipole momenta
    std::vector<double> d(3 * N);
    std::vector<double> Q(9 * N);

    for (int i = 0; i < N; ++i) {
      int ix = i % cry_l;
      int iy = (i % (cry_l * cry_l)) / cry_l;
      int iz = i / (cry_l * cry_l);

      xyz.at(3 * i + 0) = (double)(ix)*0.5;
      xyz.at(3 * i + 1) = (double)(iy)*0.5;
      xyz.at(3 * i + 2) = (double)(iz)*0.5;

      q.at(i) = (double)((ix + iy + iz) % 2) ? 1.0 : -1.0;
    }

    for (auto i : d) i = 0.0;

    for (auto i : Q) i = 0.0;

    double alpha = 1.02113246946;
    double r_max = 3.64;
    double k_max = std::pow(7.59093986701, 2.0);

    // testing the initialization with different data types
    QDMEwald<double> qdme_d(alpha, k_max, r_max, l);
    QDMEwald<float> qdme_f((float)alpha, (float)k_max, (float)r_max, (float)l);
    QDMEwald<long double> qdme_ld((long double)alpha, (long double)k_max,
                                  (long double)r_max, (long double)l);

    qdme_d.compute(xyz, q, d, Q);

    std::cout << qdme_d.get_energy() << std::endl;
  }
  Kokkos::finalize();
}
