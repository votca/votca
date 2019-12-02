#include "QDMEwald.hpp"

#include <iostream>
#include <vector>

int main() {
  Kokkos::initialize();
  {

    // testing system size
    constexpr int cry_l = 16;
    constexpr int N = cry_l * cry_l * cry_l;
    constexpr double l = (double)cry_l / 2.0;

    // testing particle positions
    std::vector<double> xyz(3 * N);

    // testing charges
    std::vector<double> q(N);

    // testing dipole momenta
    std::vector<double> d(3 * N);
    std::vector<double> Q(9 * N);

    for (int i = 0; i < N; ++i) {
      xyz.at(3 * i + 0) = (double)(i % cry_l) * 0.5;
      xyz.at(3 * i + 1) = (double)((i % (cry_l * cry_l)) / cry_l) * 0.5;
      xyz.at(3 * i + 2) = (double)(i / (cry_l * cry_l)) * 0.5;

      q.at(i) = -1 + (double)(i % 2) * 2.0;

      std::cout << i << ": " << xyz.at(3 * i + 0) << " " << xyz.at(3 * i + 1)
                << " " << xyz.at(3 * i + 2) << " " << q.at(i) << " "
                << std::endl;
    }

    for (auto i : d) i = 0.0;

    for (auto i : Q) i = 0.0;

    double alpha = 1.02113246946;
    double r_max = 3.64;
    double k_max = 7.59093986701;

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
