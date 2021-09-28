#include "votca/xtp/pmdecomposition.h"
#include "votca/xtp/aomatrix.h"
#include <votca/tools/eigenio_matrixmarket.h>

namespace votca {
namespace xtp {
void PMDecomposition::compute() {
  Eigen::MatrixXd mo_coeff = orbitals.MOs().eigenvectors().leftCols(orbitals.getNumberOfAlphaElectrons());
  QMMolecule mol = orbitals.QMAtoms();
  basis.Load(orbitals.getDFTbasisName());
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Eigen::MatrixXd S = overlap.Matrix();
  double diff_D = 10000;
  Index i = 1;
  while (diff_D > 1e-6) {
    XTP_LOG(Log::error, log) << "Iteration: " << i << std::flush;
    Eigen::MatrixXd uppertriangular = orbitalselections(mo_coeff, S);
    Index maxrow, maxcol;
    diff_D = uppertriangular.maxCoeff(&maxrow, &maxcol);
    XTP_LOG(Log::error, log) << maxrow << " " << maxcol << std::flush;
    XTP_LOG(Log::error, log) << diff_D << std::flush;
    Eigen::MatrixXd max_orbs(mo_coeff.rows(), 2);
    max_orbs << mo_coeff.col(maxrow), mo_coeff.col(maxcol);
    Eigen::MatrixXd new_orbs(mo_coeff.rows(), 2);
    new_orbs = rotatedorbitals(max_orbs, maxrow, maxcol);
    update_maximums(mo_coeff, maxrow, maxcol, new_orbs);

    i +=1;
  }
  orbitals.MOs().eigenvectors().leftCols(orbitals.getNumberOfAlphaElectrons()) = mo_coeff;
  Eigen::MatrixXd check = orbitals.MOs().eigenvectors().transpose() * S * orbitals.MOs().eigenvectors();
  XTP_LOG(Log::error, log) << check << std::flush;
  orbitals.WriteToCpt("pm_orbitals.orb");
  votca::tools::EigenIO_MatrixMarket::WriteMatrix(
         "test_pm_orbitals.mm", mo_coeff);
}

// Eigen::MatrixXd PMDecomposition::columnwise(const Eigen::MatrixXd &S, Eigen::VectorXd &v) {
//   Eigen::MatrixXd a(S.rows(), S.cols());
//   for (int p = 0; p < S.rows(); p++) {
//     a.col(p) = v(p) * S.col(p);
//   }
//   return a;
// }

// Eigen::MatrixXd PMDecomposition::rowwise(const Eigen::MatrixXd &S, Eigen::VectorXd &v) {
//   Eigen::MatrixXd a(S.rows(), S.cols());
//   for (int p = 0; p < S.rows(); p++) {
//     a.row(p) = v(p) * S.row(p);
//   }
//   return a;
// }

Eigen::MatrixXd PMDecomposition::rotatedorbitals(Eigen::MatrixXd &maxorbs, Index s, Index t) {
  Eigen::VectorXd vec1, vec2, new_vec1, new_vec2;
  double gam, sin_gamma, cos_gamma;
  Eigen::MatrixXd neworbitals(maxorbs.rows(), 2);
  vec1 = maxorbs.col(0);
  vec2 = maxorbs.col(1);
  gam = 0.25 * asin(B(s,t) / sqrt((A(s,t) * A(s,t)) + (B(s,t) * B(s,t))));
  sin_gamma = std::sin(gam);
  cos_gamma = std::cos(gam);
  new_vec1 = (cos_gamma * vec1) + (sin_gamma * vec2);
  new_vec2 = -1 * (sin_gamma * vec1) + (cos_gamma * vec2);
  neworbitals << new_vec1, new_vec2;
  XTP_LOG(Log::error, log) << sin_gamma << std::flush;
  return neworbitals;
}

//Function to select n(n-1)/2 orbitals and process Ast and Bst
Eigen::MatrixXd PMDecomposition::orbitalselections(Eigen::MatrixXd &m, const Eigen::MatrixXd &S) {
  Eigen::VectorXd req_vec1, req_vec2;
  Eigen::RowVectorXd spt, sps, tpt;
  Eigen::MatrixXd req_mat(m.rows(), 2), a, b, c, d, e, f;
  Eigen::MatrixXd zeromatrix = Eigen::MatrixXd::Zero(m.cols(), m.cols());
  A = Eigen::MatrixXd::Zero(m.cols(), m.cols());
  B = Eigen::MatrixXd::Zero(m.cols(), m.cols());
  for (int s = 0; s < m.cols(); s++) {
    for (int t = 0; t < m.cols(); t++) {
      if (t > s) {
        req_vec1 = m.col(s);
        req_vec2 = m.col(t);
        a = S*req_vec1.asDiagonal();
        b = S*req_vec2.asDiagonal();
        c = req_vec1.transpose()*a; //vec1.S.vec1
        d = req_vec1.transpose()*b; //term1 of eq 31
        e = req_vec1.asDiagonal()*b; //term2 of eq 31
        f = req_vec2.transpose()*b; //vec2.S.vec2
        sps = c.colwise().sum();
        tpt = f.colwise().sum();
        spt = 0.5 * (d.colwise().sum() + e.rowwise().sum().transpose());
        std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
        Index start = 0;
        double Ast = 0;
        double Bst = 0;
        for (Index atom_id = 0; atom_id < Index(numfuncpatom.size()); atom_id++)
        {
            double sps_x = sps.segment(start,numfuncpatom[atom_id]).sum();
            double spt_x = spt.segment(start,numfuncpatom[atom_id]).sum();
            double tpt_x = tpt.segment(start,numfuncpatom[atom_id]).sum();
            Ast += spt_x * spt_x - 0.25 * ((sps_x - tpt_x) * (sps_x - tpt_x));
            Bst += spt_x * (sps_x - tpt_x);
            start += numfuncpatom[atom_id];
        }
        A(s,t) = Ast;
        B(s,t) = Bst;
        double parameter = Ast + sqrt((Ast * Ast) + (Bst * Bst));
        zeromatrix(s,t) = parameter;
      }
    }
  }
  return zeromatrix;
}

void PMDecomposition::update_maximums(Eigen::MatrixXd &m, Index col1, Index col2, Eigen::MatrixXd &new_orbs) {
  m.col(col1) = new_orbs.col(0);
  m.col(col2) = new_orbs.col(1);
}
}  // namespace xtp
}  // namespace votca