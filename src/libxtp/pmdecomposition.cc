#include "votca/xtp/pmdecomposition.h"
#include "votca/xtp/aomatrix.h"

namespace votca {
namespace xtp {
void PMDecomposition::compute() {
  Eigen::MatrixXd mo_coeff = orbitals.MOs().eigenvectors();
  QMMolecule mol = orbitals.QMAtoms();
  basis.Load(orbitals.getDFTbasisName());
  aobasis.Fill(basis, mol);
  AOOverlap overlap;
  overlap.Fill(aobasis);
  Eigen::MatrixXd S = overlap.Matrix();
  double diff_D = 100;
  while (diff_D > 1) {
    Eigen::MatrixXd uppertriangular = orbitalselections(mo_coeff, S);
    XTP_LOG(Log::error, log) << uppertriangular << std::flush;
    Index maxrow, maxcol;
    diff_D = uppertriangular.maxCoeff(&maxrow, &maxcol);
    XTP_LOG(Log::error, log) << maxrow << " " << maxcol << std::flush;
    XTP_LOG(Log::error, log) << diff_D << std::flush;
    Eigen::MatrixXd max_orbs(mo_coeff.rows(), 2);
    max_orbs << mo_coeff.col(maxrow), mo_coeff.col(maxcol);
    Eigen::MatrixXd new_orbs(mo_coeff.rows(), 2);
    new_orbs = rotatedorbitals(max_orbs, maxrow, maxcol);
    update_maximums(mo_coeff, maxrow, maxcol, new_orbs);
  }
}

Eigen::MatrixXd PMDecomposition::columnwise(const Eigen::MatrixXd &S, Eigen::VectorXd &v) {
  Eigen::MatrixXd a(S.rows(), S.cols());
  for (int p = 0; p < S.rows(); p++) {
    a.col(p) = v(p) * S.col(p);
  }
  return a;
}

Eigen::MatrixXd PMDecomposition::rowwise(const Eigen::MatrixXd &S, Eigen::VectorXd &v) {
  Eigen::MatrixXd a(S.rows(), S.cols());
  for (int p = 0; p < S.rows(); p++) {
    a.row(p) = v(p) * S.row(p);
  }
  return a;
}

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
        req_mat << req_vec1, req_vec2;
        a = columnwise(S,req_vec1);
        b = columnwise(S,req_vec2);
        c = rowwise(a, req_vec2);
        d = rowwise(b, req_vec1);
        e = rowwise(a, req_vec1);
        f = rowwise(b, req_vec2);
        sps = e.colwise().sum();
        tpt = f.colwise().sum();
        spt = 0.5 * (c.colwise().sum() + d.colwise().sum());
        std::vector<Index> numfuncpatom = aobasis.getFuncPerAtom();
        Eigen::RowVectorXd sps_per_atom(numfuncpatom.size()), spt_per_atom(numfuncpatom.size()), tpt_per_atom(numfuncpatom.size());
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
        A(t,s) = Ast;
        B(t,s) = Bst;
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