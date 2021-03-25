#include "votca/xtp/pmdecomposition.h"

namespace votca {
namespace xtp {
    void PMDecomposition::compute() {
    Eigen::MatrixXd mo_coeff = orbitals.MOs().eigenvectors();
    std::cout << mo_coeff << std::endl;
}
}
}