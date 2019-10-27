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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include "votca/xtp/aoshell.h"
#include "votca/xtp/aobasis.h"
#include <votca/xtp/aomatrix.h>

namespace votca {
namespace xtp {

void AOShell::normalizeContraction() {
  AOOverlap overlap;
  Eigen::MatrixXd block = overlap.FillShell(*this);
  std::vector<int> numsubshells = NumFuncSubShell(_type);
  int contraction_index = FindLmin(_type);
  int aoindex = 0;
  for (int& numsubshell : numsubshells) {
    double norm = std::sqrt(block(aoindex, aoindex));
    for (auto& gaussian : _gaussians) {
      gaussian._contraction[contraction_index] /= norm;
    }
    aoindex += numsubshell;
    contraction_index++;
  }
  return;
}

void AOShell::EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                          Eigen::Block<Eigen::MatrixX3d>& gradAOvalues,
                          const Eigen::Vector3d& grid_pos) const {

  // need position of shell
  const Eigen::Vector3d center = (grid_pos - _pos);
  const double& center_x = center.x();
  const double& center_y = center.y();
  const double& center_z = center.z();
  const double distsq = center.squaredNorm();

  // iterate over Gaussians in this shell
  for (const AOGaussianPrimitive& gaussian : _gaussians) {

    const double alpha = gaussian.getDecay();
    const Eigen::VectorXd& contractions = gaussian.getContraction();

    const double expofactor =
        gaussian.getPowfactor() * std::exp(-alpha * distsq);
    const Eigen::Vector3d second_term = -2.0 * alpha * center;
    // split combined shells
    int i_func = -1;
    int i_act;
    for (const char& single_shell : _type) {
      // single type shells
      if (single_shell == 'S') {

        i_act = i_func + 1;
        double AOvalue = contractions[0] * expofactor;
        AOvalues(i_act) += AOvalue;  // s-function
        gradAOvalues.row(i_act) +=
            second_term * AOvalue;  // gradient of s-function
        i_func++;
      } else if (single_shell == 'P') {
        const double factor = 2. * sqrt(alpha) * contractions[1] * expofactor;

        i_act = i_func + 1;                  // Y 1,0
        double AOvalue = factor * center_z;  // Y 1,0
        AOvalues(i_act) += AOvalue;
        gradAOvalues.row(i_act) += second_term * AOvalue;
        gradAOvalues(i_act, 2) += factor;

        i_act++;                      // Y 1,-1
        AOvalue = factor * center_y;  // Y 1,-1
        AOvalues(i_act) += AOvalue;
        gradAOvalues.row(i_act) += second_term * AOvalue;
        gradAOvalues(i_act, 1) += factor;

        i_act++;                      // Y 1,1
        AOvalue = factor * center_x;  // Y 1,1
        AOvalues(i_act) += AOvalue;
        gradAOvalues(i_act, 0) += factor;
        gradAOvalues.row(i_act) += second_term * AOvalue;  // y gradient

        i_func += 3;
      } else if (single_shell == 'D') {
        const double factor = 2. * alpha * contractions[2] * expofactor;
        const double factor_1 = factor / sqrt(3.);

        i_act = i_func + 1;  // Y 2,0
        double AOvalue =
            factor_1 * (3. * center_z * center_z - distsq);  // Y 2,0
        AOvalues(i_act) += AOvalue;
        Eigen::Array3d coeff(-2, -2, 4);
        gradAOvalues.row(i_act) +=
            (factor_1 * coeff * center.array()).matrix() +
            second_term * AOvalue;

        i_act++;                                        // Y 2,-1
        AOvalue = 2. * factor * (center_y * center_z);  // Y 2,-1
        AOvalues(i_act) += AOvalue;
        coeff = {0, 2 * center_z, 2 * center_y};
        gradAOvalues.row(i_act) +=
            factor * coeff.matrix() + second_term * AOvalue;

        i_act++;                                        // Y 2,1
        AOvalue = 2. * factor * (center_x * center_z);  // Y 2,1
        AOvalues(i_act) += AOvalue;
        coeff = {2 * center_z, 0, 2 * center_x};
        gradAOvalues.row(i_act) +=
            factor * coeff.matrix() + second_term * AOvalue;

        i_act++;                                        // Y 2,-2
        AOvalue = 2. * factor * (center_x * center_y);  // Y 2,-2
        AOvalues(i_act) += AOvalue;
        coeff = {2 * center_y, 2 * center_x, 0};
        gradAOvalues.row(i_act) +=
            factor * coeff.matrix() + second_term * AOvalue;

        i_act++;  // Y 2,2
        AOvalue =
            factor * (center_x * center_x - center_y * center_y);  // Y 2,2
        AOvalues(i_act) += AOvalue;
        coeff = {2 * center_x, -2 * center_y, 0};
        gradAOvalues.row(i_act) +=
            factor * coeff.matrix() + second_term * AOvalue;

        i_func += 5;
      } else if (single_shell == 'F') {
        const double factor =
            2. * pow(alpha, 1.5) * contractions[3] * expofactor;
        const double factor_1 = factor * 2. / sqrt(15.);
        const double factor_2 = factor * sqrt(2.) / sqrt(5.);
        const double factor_3 = factor * sqrt(2.) / sqrt(3.);
        const double cx_cx = center_x * center_x;
        const double cx_cy = center_x * center_y;
        const double cx_cz = center_x * center_z;
        const double cy_cy = center_y * center_y;
        const double cy_cz = center_y * center_z;
        const double cz_cz = center_z * center_z;

        i_act = i_func + 1;  // Y 3,0
        double AOvalue =
            factor_1 * center_z * (5. * cz_cz - 3. * distsq);  // Y 3,0
        AOvalues(i_act) += AOvalue;
        Eigen::Array3d coeff = {-6. * cx_cz, -6. * cy_cz,
                                3. * (3. * cz_cz - distsq)};
        gradAOvalues.row(i_act) +=
            factor_1 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                // Y 3,-1
        AOvalue = factor_2 * center_y * (5. * cz_cz - distsq);  // Y 3,-1
        AOvalues(i_act) += AOvalue;
        coeff = {-2. * cx_cy, 4. * cz_cz - cx_cx - 3. * cy_cy, 8. * cy_cz};
        gradAOvalues.row(i_act) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                // Y 3,1
        AOvalue = factor_2 * center_x * (5. * cz_cz - distsq);  // Y 3,1
        AOvalues(i_act) += AOvalue;
        coeff = {4. * cz_cz - cy_cy - 3. * cx_cx, -2. * cx_cy, 8. * cx_cz};
        gradAOvalues.row(i_act) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                 // Y 3,-2
        AOvalue = 4. * factor * center_x * center_y * center_z;  // Y 3,-2
        AOvalues(i_act) += AOvalue;
        coeff = {cy_cz, cx_cz, cx_cy};
        gradAOvalues.row(i_act) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;

        i_act++;                                             // Y 3,2
        AOvalue = 2. * factor * center_z * (cx_cx - cy_cy);  // Y 3,2
        AOvalues(i_act) += AOvalue;
        coeff = {2. * cx_cz, -2. * cy_cz, cx_cx - cy_cy};
        gradAOvalues.row(i_act) +=
            2 * factor * coeff.matrix() + second_term * AOvalue;

        i_act++;                                               // Y 3,-3
        AOvalue = factor_3 * center_y * (3. * cx_cx - cy_cy);  // Y 3,-3
        AOvalues(i_act) += AOvalue;
        coeff = {6. * cx_cy, 3. * (cx_cx - cy_cy), 0};
        gradAOvalues.row(i_act) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                               // Y 3,3
        AOvalue = factor_3 * center_x * (cx_cx - 3. * cy_cy);  // Y 3,3
        AOvalues(i_act) += AOvalue;
        coeff = {3. * (cx_cx - cy_cy), -6. * cx_cy, 0};
        gradAOvalues.row(i_act) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;

        i_func += 7;
      } else if (single_shell == 'G') {
        const double factor =
            2. / sqrt(3.) * alpha * alpha * contractions[4] * expofactor;
        const double factor_1 = factor / sqrt(35.);
        const double factor_2 = factor * 4. / sqrt(14.);
        const double factor_3 = factor * 2. / sqrt(7.);
        const double factor_4 = factor * 2. * sqrt(2.);
        const double cx_cx = center_x * center_x;
        const double cx_cy = center_x * center_y;
        const double cx_cz = center_x * center_z;
        const double cy_cy = center_y * center_y;
        const double cy_cz = center_y * center_z;
        const double cz_cz = center_z * center_z;

        i_act = i_func + 1;  // Y 4,0
        double AOvalue =
            factor_1 * (35. * cz_cz * cz_cz - 30. * cz_cz * distsq +
                        3. * distsq * distsq);  // Y 4,0
        AOvalues(i_act) += AOvalue;
        Eigen::Array3d coeff = {12. * center_x * (distsq - 5. * cz_cz),
                                12. * center_y * (distsq - 5. * cz_cz),
                                16. * center_z * (5. * cz_cz - 3. * distsq)};

        gradAOvalues.row(i_act) +=
            factor_1 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                  // Y 4,-1
        AOvalue = factor_2 * cy_cz * (7. * cz_cz - 3. * distsq);  // Y 4,-1
        AOvalues(i_act) += AOvalue;
        coeff = {(-6. * center_x * cy_cz),
                 center_z * (4. * cz_cz - 3. * cx_cx - 9. * cy_cy),
                 3. * center_y * (5. * cz_cz - distsq)};
        gradAOvalues.row(i_act) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                  // Y 4,1
        AOvalue = factor_2 * cx_cz * (7. * cz_cz - 3. * distsq);  // Y 4,1
        AOvalues(i_act) += AOvalue;
        coeff = {center_z * (4. * cz_cz - 9. * cx_cx - 3. * cy_cy),
                 (-6. * center_y * cx_cz),
                 3. * center_x * (5. * cz_cz - distsq)};
        gradAOvalues.row(i_act) +=
            factor_2 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                  // Y 4,-2
        AOvalue = 2. * factor_3 * cx_cy * (7. * cz_cz - distsq);  // Y 4,-2
        AOvalues(i_act) += AOvalue;
        coeff = {center_y * (6. * cz_cz - 3. * cx_cx - cy_cy),
                 center_x * (6. * cz_cz - cx_cx - 3. * cy_cy),
                 12. * center_z * cx_cy};
        gradAOvalues.row(i_act) +=
            2 * factor_3 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                                       // Y 4,2
        AOvalue = factor_3 * (cx_cx - cy_cy) * (7. * cz_cz - distsq);  // Y 4,2
        AOvalues(i_act) += AOvalue;
        coeff = {4. * center_x * (3. * cz_cz - cx_cx),
                 4. * center_y * (cy_cy - 3. * cz_cz),
                 12. * center_z * (cx_cx - cy_cy)};
        gradAOvalues.row(i_act) +=
            factor_3 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                            // Y 4,-3
        AOvalue = factor_4 * cy_cz * (3. * cx_cx - cy_cy);  // Y 4,-3
        AOvalues(i_act) += AOvalue;
        coeff = {6. * center_x * cy_cz, 3. * center_z * (cx_cx - cy_cy),
                 center_y * (3. * cx_cx - cy_cy)};
        gradAOvalues.row(i_act) +=
            factor_4 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                            // Y 4,3
        AOvalue = factor_4 * cx_cz * (cx_cx - 3. * cy_cy);  // Y 4,3
        AOvalues(i_act) += AOvalue;
        coeff = {3. * center_z * (cx_cx - cy_cy), (-6. * center_y * cx_cz),
                 center_x * (cx_cx - 3. * cy_cy)};
        gradAOvalues.row(i_act) +=
            factor_4 * coeff.matrix() + second_term * AOvalue;

        i_act++;                                          // Y 4,-4
        AOvalue = 4. * factor * cx_cy * (cx_cx - cy_cy);  // Y 4,-4
        AOvalues(i_act) += AOvalue;
        coeff = {center_y * (3. * cx_cx - cy_cy),
                 center_x * (cx_cx - 3. * cy_cy), 0};
        gradAOvalues.row(i_act) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;

        i_act++;  // Y 4,4
        AOvalue = factor * (cx_cx * cx_cx - 6. * cx_cx * cy_cy +
                            cy_cy * cy_cy);  // Y 4,4
        AOvalues(i_act) += AOvalue;
        coeff = {center_x * (cx_cx - 3. * cy_cy),
                 center_y * (cy_cy - 3. * cx_cx), 0};
        gradAOvalues.row(i_act) +=
            4 * factor * coeff.matrix() + second_term * AOvalue;

        i_func += 9;
      } else {
        throw std::runtime_error("Shell type:" + std::string(1, single_shell) +
                                 " not known");
      }
    }
  }  // contractions
  return;
}

void AOShell::EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>& AOvalues,
                          const Eigen::Vector3d& grid_pos) const {

  // need position of shell
  const Eigen::Vector3d center = grid_pos - _pos;
  const double center_x = center[0];
  const double center_y = center[1];
  const double center_z = center[2];
  const double distsq = center.squaredNorm();

  // iterate over Gaussians in this shell
  for (const AOGaussianPrimitive& gaussian : _gaussians) {

    const double alpha = gaussian.getDecay();
    const Eigen::VectorXd& contractions = gaussian.getContraction();

    double expofactor = gaussian.getPowfactor() * exp(-alpha * distsq);

    // split combined shells
    int i_func = -1;

    for (const char& single_shell : _type) {
      // single type shells
      if (single_shell == 'S') {
        AOvalues(i_func + 1) += contractions[0] * expofactor;  // s-function
        i_func++;
      } else if (single_shell == 'P') {
        double factor = 2. * sqrt(alpha) * contractions[1];
        AOvalues(i_func + 1) += factor * center_z * expofactor;  // Y 1,0
        AOvalues(i_func + 2) += factor * center_y * expofactor;  // Y 1,-1
        AOvalues(i_func + 3) += factor * center_x * expofactor;  // Y 1,1
        i_func += 3;
      } else if (single_shell == 'D') {
        double factor = 2. * alpha * contractions[2];
        double factor_1 = factor / sqrt(3.);
        AOvalues(i_func + 1) += factor_1 * (3. * center_z * center_z - distsq) *
                                expofactor;  // Y 2,0
        AOvalues(i_func + 2) +=
            2. * factor * (center_y * center_z) * expofactor;  // Y 2,-1
        AOvalues(i_func + 3) +=
            2. * factor * (center_x * center_z) * expofactor;  // Y 2,1
        AOvalues(i_func + 4) +=
            2. * factor * (center_x * center_y) * expofactor;  // Y 2,-2
        AOvalues(i_func + 5) += factor *
                                (center_x * center_x - center_y * center_y) *
                                expofactor;  // Y 2,2
        i_func += 5;
      } else if (single_shell == 'F') {
        double factor = 2. * pow(alpha, 1.5) * contractions[3];
        double factor_1 = factor * 2. / sqrt(15.);
        double factor_2 = factor * sqrt(2.) / sqrt(5.);
        double factor_3 = factor * sqrt(2.) / sqrt(3.);
        double cx_cx = center_x * center_x;
        double cy_cy = center_y * center_y;
        double cz_cz = center_z * center_z;

        AOvalues(i_func + 1) += factor_1 * center_z *
                                (5. * cz_cz - 3. * distsq) *
                                expofactor;  // Y 3,0
        AOvalues(i_func + 2) +=
            factor_2 * center_y * (5. * cz_cz - distsq) * expofactor;  // Y 3,-1
        AOvalues(i_func + 3) +=
            factor_2 * center_x * (5. * cz_cz - distsq) * expofactor;  // Y 3,1
        AOvalues(i_func + 4) += 4. * factor * center_x * center_y * center_z *
                                expofactor;  // Y 3,-2
        AOvalues(i_func + 5) +=
            2. * factor * center_z * (cx_cx - cy_cy) * expofactor;  // Y 3,2
        AOvalues(i_func + 6) +=
            factor_3 * center_y * (3. * cx_cx - cy_cy) * expofactor;  // Y 3,-3
        AOvalues(i_func + 7) +=
            factor_3 * center_x * (cx_cx - 3. * cy_cy) * expofactor;  // Y 3,3
        i_func += 7;
      } else if (single_shell == 'G') {
        double factor = 2. / sqrt(3.) * alpha * alpha * contractions[4];
        double factor_1 = factor / sqrt(35.);
        double factor_2 = factor * 4. / sqrt(14.);
        double factor_3 = factor * 2. / sqrt(7.);
        double factor_4 = factor * 2. * sqrt(2.);
        double cx_cx = center_x * center_x;
        double cx_cy = center_x * center_y;
        double cx_cz = center_x * center_z;
        double cy_cy = center_y * center_y;
        double cy_cz = center_y * center_z;
        double cz_cz = center_z * center_z;

        AOvalues(i_func + 1) += factor_1 *
                                (35. * cz_cz * cz_cz - 30. * cz_cz * distsq +
                                 3. * distsq * distsq) *
                                expofactor;  // Y 4,0
        AOvalues(i_func + 2) += factor_2 * cy_cz * (7. * cz_cz - 3. * distsq) *
                                expofactor;  // Y 4,-1
        AOvalues(i_func + 3) += factor_2 * cx_cz * (7. * cz_cz - 3. * distsq) *
                                expofactor;  // Y 4,1
        AOvalues(i_func + 4) += 2. * factor_3 * cx_cy * (7. * cz_cz - distsq) *
                                expofactor;  // Y 4,-2
        AOvalues(i_func + 5) += factor_3 * (cx_cx - cy_cy) *
                                (7. * cz_cz - distsq) * expofactor;  // Y 4,2
        AOvalues(i_func + 6) +=
            factor_4 * cy_cz * (3. * cx_cx - cy_cy) * expofactor;  // Y 4,-3
        AOvalues(i_func + 7) +=
            factor_4 * cx_cz * (cx_cx - 3. * cy_cy) * expofactor;  // Y 4,3
        AOvalues(i_func + 8) +=
            4. * factor * cx_cy * (cx_cx - cy_cy) * expofactor;  // Y 4,-4
        AOvalues(i_func + 9) +=
            factor * (cx_cx * cx_cx - 6. * cx_cx * cy_cy + cy_cy * cy_cy) *
            expofactor;  // Y 4,4

        i_func += 9;
      } else {
        throw std::runtime_error("Single shell type " +
                                 std::string(1, single_shell) + " not known ");
      }
    }
  }  // contractions
  return;
}

std::ostream& operator<<(std::ostream& out, const AOShell& shell) {
  out << "AtomIndex:" << shell.getAtomIndex();
  out << " Shelltype:" << shell.getType() << " Scale:" << shell.getScale()
      << " Lmax:" << shell.getLmax() << " MinDecay:" << shell.getMinDecay()
      << " Func:" << shell.getNumFunc() << "\n";
  for (const auto& gaussian : shell) {
    out << " Gaussian Decay: " << gaussian.getDecay();
    out << " Contractions:";
    for (long i = 0; i < gaussian.getContraction().size(); i++) {
      out << " " << gaussian.getContraction()[i];
    }
    out << "\n";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
