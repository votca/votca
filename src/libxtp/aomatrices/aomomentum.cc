/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aomatrix3d.h>
#include <votca/xtp/aotransform.h>
namespace votca {
namespace xtp {

void AOMomentum::FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd> >& matrix,
                           const AOShell& shell_row,
                           const AOShell& shell_col) const {

  /* Calculating the AO matrix of the gradient operator requires
   * the raw overlap matrix (i.e. in unnormalized cartesians)
   * with lmax of the column shell increased by one:
   *
   *         phi_(ijk) = x^i y^j z^k exp(-beta r^2)
   * => d/dx phi_(ijk) = (i*x^(i-1) - 2beta x^(i+1)) y^j z^k exp(-beta r^2)
   *    d/dy phi_(ijk) = x^i (j*y^(j-1) - 2beta y^(j+1)) z^k exp(-beta r^2)
   *    d/dz phi_(ijk) = x^i y^j (k*z^(k-1) - 2beta z^(k+1)) exp(-beta r^2)
   *
   * i.e.:   d/dx phi_s  = d/dx phi_(000) = -2beta phi_(100) = -2beta phi_px
   *         d/dy phi_px = d/dy phi_(100) = -2beta phi_(110) = -2beta phi_dxy
   *         d/dz phi_pz = d/dz phi_(001) = phi_(000) - 2beta phi_(002)
   *                                      = phi_s     - 2beta phi_dxx
   *
   * and with that
   *         <s|d/dx|s>  = -2beta <s|px>
   *         <s|d/dy|px> = -2beta <s|dxy>
   *         <s|d/dz|pz> = <s|s> - 2beta <s|dxx>
   *         ...
   *
   */

  // shell info, only lmax tells how far to go
  int lmax_row = shell_row.getLmax();
  int lmax_col = shell_col.getLmax();

  if (lmax_col > 4) {
    throw std::runtime_error(
        "Momentum transition dipoles only implemented for S,P,D,F,G functions "
        "in DFT basis!");
  }

  // set size of internal block for recursion
  int nrows = AOTransform::getBlockSize(lmax_row);
  int ncols = AOTransform::getBlockSize(lmax_col);

  // initialize local matrix block for unnormalized cartesians
  std::array<Eigen::MatrixXd, 3> mom;
  for (int i_comp = 0; i_comp < 3; i_comp++) {
    mom[i_comp] = Eigen::MatrixXd ::Zero(nrows, ncols);
  }

  std::array<Eigen::MatrixXd, 6> scd_mom;
  for (int i_comp = 0; i_comp < 6; i_comp++) {
    scd_mom[i_comp] = Eigen::MatrixXd ::Zero(nrows, ncols);
  }

  const Eigen::Vector3d& pos_row = shell_row.getPos();
  const Eigen::Vector3d& pos_col = shell_col.getPos();
  const Eigen::Vector3d diff = pos_row - pos_col;

  double distsq = diff.squaredNorm();

  std::array<int, 165> nx = AOTransform::nx();
  std::array<int, 165> ny = AOTransform::ny();
  std::array<int, 165> nz = AOTransform::nz();
  std::array<int, 165> i_less_x = AOTransform::i_less_x();
  std::array<int, 165> i_less_y = AOTransform::i_less_y();
  std::array<int, 165> i_less_z = AOTransform::i_less_z();
  std::array<int, 120> i_more_x = AOTransform::i_more_x();
  std::array<int, 120> i_more_y = AOTransform::i_more_y();
  std::array<int, 120> i_more_z = AOTransform::i_more_z();

  for (const auto& gaussian_row : shell_row) {
    const double decay_row = gaussian_row.getDecay();

    for (const auto& gaussian_col : shell_col) {

      const double decay_col = gaussian_col.getDecay();

      const double fak = 0.5 / (decay_row + decay_col);
      const double fak2 = 2.0 * fak;
      double _exparg = fak2 * decay_row * decay_col * distsq;

      /// check if distance between positions is large, then skip step

      if (_exparg > 30.0) {
        continue;
      }

      AOOverlap overlap;
      int L_offset = 1;
      Eigen::MatrixXd ol =
          overlap.Primitive_Overlap(gaussian_row, gaussian_col, L_offset);

      double alpha2 = 2.0 * decay_row;
      double beta2 = 2.0 * decay_col;
      for (int i = 0; i < ncols; i++) {

        int nx_i = nx[i];
        int ny_i = ny[i];
        int nz_i = nz[i];
        int ilx_i = i_less_x[i];
        int ily_i = i_less_y[i];
        int ilz_i = i_less_z[i];
        int imx_i = i_more_x[i];
        int imy_i = i_more_y[i];
        int imz_i = i_more_z[i];

        for (int j = 0; j < nrows; j++) {

          mom[0](j, i) = nx_i * ol(j, ilx_i) - beta2 * ol(j, imx_i);
          mom[1](j, i) = ny_i * ol(j, ily_i) - beta2 * ol(j, imy_i);
          mom[2](j, i) = nz_i * ol(j, ilz_i) - beta2 * ol(j, imz_i);

          int nx_j = nx[j];
          int ny_j = ny[j];
          int nz_j = nz[j];
          int ilx_j = i_less_x[j];
          int ily_j = i_less_y[j];
          int ilz_j = i_less_z[j];
          int imx_j = i_more_x[j];
          int imy_j = i_more_y[j];
          int imz_j = i_more_z[j];

          scd_mom[0](j, i) =
              nx_j * (beta2 * ol(ilx_j, imx_i) - nx_i * ol(ilx_j, ilx_i)) -
              alpha2 * (beta2 * ol(imx_j, imx_i) -
                        nx_i * ol(imx_j, ilx_i));  // d2/(dxdx)
          scd_mom[1](j, i) =
              nx_j * (beta2 * ol(ilx_j, imy_i) - ny_i * ol(ilx_j, ily_i)) -
              alpha2 * (beta2 * ol(imx_j, imy_i) -
                        ny_i * ol(imx_j, ily_i));  // d2/(dxdy)
          scd_mom[2](j, i) =
              nx_j * (beta2 * ol(ilx_j, imz_i) - nz_i * ol(ilx_j, ilz_i)) -
              alpha2 * (beta2 * ol(imx_j, imz_i) -
                        nz_i * ol(imx_j, ilz_i));  // d2/(dxdz)

          scd_mom[3](j, i) =
              ny_j * (beta2 * ol(ily_j, imy_i) - ny_i * ol(ily_j, ily_i)) -
              alpha2 * (beta2 * ol(imy_j, imy_i) -
                        ny_i * ol(imy_j, ily_i));  // d2/(dydy)
          scd_mom[4](j, i) =
              ny_j * (beta2 * ol(ily_j, imz_i) - nz_i * ol(ily_j, ilz_i)) -
              alpha2 * (beta2 * ol(imy_j, imz_i) -
                        nz_i * ol(imy_j, ilz_i));  // d2/(dydz)

          scd_mom[5](j, i) =
              nz_j * (beta2 * ol(ilz_j, imz_i) - nz_i * ol(ilz_j, ilz_i)) -
              alpha2 * (beta2 * ol(imz_j, imz_i) -
                        nz_i * ol(imz_j, ilz_i));  // d2/(dzdz)
        }
      }

      Eigen::MatrixXd trafo_row = AOTransform::getTrafo(gaussian_row);
      Eigen::MatrixXd trafo_col = AOTransform::getTrafo(gaussian_col);
      // cartesian -> spherical
      for (int i = 0; i < 3; i++) {
        Eigen::MatrixXd mom_sph =
            trafo_row.transpose() *
            mom[i].bottomRightCorner(shell_row.getCartesianNumFunc(),
                                     shell_col.getCartesianNumFunc()) *
            trafo_col;
        // save to matrix
        matrix[i] += mom_sph;
      }

    }  // shell_col Gaussians
  }    // shell_row Gaussians
}

}  // namespace xtp
}  // namespace votca
