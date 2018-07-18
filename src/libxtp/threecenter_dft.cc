/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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



#include <votca/xtp/threecenter.h>
#include <votca/xtp/symmetric_matrix.h>
#include <votca/xtp/eigen.h>
namespace votca {
  namespace xtp {

    void TCMatrix_dft::Fill(AOBasis& _auxbasis, AOBasis& _dftbasis, const Eigen::MatrixXd& V_sqrtm1) {

      for (unsigned int i = 0; i < _auxbasis.AOBasisSize(); i++) {
        try {
          _matrix.push_back(Symmetric_Matrix(_dftbasis.AOBasisSize()));
        } catch (std::bad_alloc& ba) {
          std::cerr << "Basisset/aux basis too large for 3c calculation. Not enough RAM. Caught bad alloc: " << ba.what() << std::endl;
          exit(0);
        }

      }
      #pragma omp parallel for schedule(dynamic)
      for (int _is = _dftbasis.getNumofShells()-1; _is >=0; _is--) {
        const Eigen::MatrixXd V=V_sqrtm1;
        const AOShell* _dftshell = _dftbasis.getShell(_is);
        std::vector< Eigen::MatrixXd > block;
        for (int i = 0; i < _dftshell->getNumFunc(); i++) {
          int size = _dftshell->getStartIndex() + i+1;
          block.push_back(Eigen::MatrixXd::Zero(_auxbasis.AOBasisSize(), size));
        }
        FillBlock(block, _is, _dftbasis, _auxbasis);
        int offset = _dftshell->getStartIndex();
        for (unsigned i = 0; i < block.size(); ++i) {
          Eigen::MatrixXd temp =V * block[i];
          for (int mu = 0; mu < temp.rows(); ++mu) {
            for (int j = 0; j < temp.cols(); ++j) {
              _matrix[mu](i + offset, j) = temp(mu, j);
            }
          }
        }
      }
      return;
    }

    /*
     * Determines the 3-center integrals for a given shell in the aux basis
     * by calculating the 3-center overlap integral of the functions in the
     * aux shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock)
     */

    void TCMatrix_dft::FillBlock(std::vector< Eigen::MatrixXd >& _block, int shellindex, const AOBasis& dftbasis, const AOBasis& auxbasis) {
      const AOShell* left_dftshell = dftbasis.getShell(shellindex);
      tensor3d::extent_gen extents;
      int _start = left_dftshell->getStartIndex();
      // alpha-loop over the aux basis function
      for (const AOShell* shell_aux:auxbasis) {
        int _aux_start = shell_aux->getStartIndex();


        for (int _is = 0; _is <= shellindex; _is++) {

          const AOShell* _shell_col = dftbasis.getShell(_is);
          int _col_start=_shell_col->getStartIndex();
          tensor3d threec_block(extents[ range(0, shell_aux->getNumFunc()) ][ range(0, left_dftshell->getNumFunc()) ][ range(0, _shell_col->getNumFunc())]);
          for (int i = 0; i < shell_aux->getNumFunc(); ++i) {
            for (int j = 0; j < left_dftshell->getNumFunc(); ++j) {
              for (int k = 0; k < _shell_col->getNumFunc(); ++k) {
                threec_block[i][j][k] = 0.0;
              }
            }
          }

          bool nonzero = FillThreeCenterRepBlock(threec_block, shell_aux, left_dftshell, _shell_col);
          if (nonzero) {

            for (int _left = 0; _left < left_dftshell->getNumFunc(); _left++) {
              for (int _aux = 0; _aux < shell_aux->getNumFunc(); _aux++) {
                for (int _col = 0; _col < _shell_col->getNumFunc(); _col++) {
                  //symmetry
                  if ((_col_start + _col)>(_start + _left)) {
                    break;
                  }
                  _block[_left](_aux_start + _aux, _col_start + _col) = threec_block[_aux][_left][_col];
                }
              }
            }
          }
        }
      }
      return;
    }






  }
}
