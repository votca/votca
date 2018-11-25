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

namespace votca {
  namespace xtp {

    void TCMatrix_gwbse::Initialize(int basissize, int mmin, int mmax, int nmin, int nmax) {

      // here as storage indices starting from zero
      _nmin = nmin;
      _nmax = nmax;
      _ntotal = nmax - nmin + 1;
      _mmin = mmin;
      _mmax = mmax;
      _mtotal = mmax - mmin + 1;
      _basissize = basissize;

      // vector has mtotal elements
      _matrix.resize(_mtotal);

      // each element is a gwabasis-by-n matrix, initialize to zero
      for (int i = 0; i < this->msize(); i++) {
        _matrix[i] = MatrixXfd::Zero(_ntotal,_basissize);
      }
      
    }

    /*
     * Cleaning TCMatrix data and free memory
     */
    void TCMatrix_gwbse::Cleanup() {

      for (unsigned i = 0; i < _matrix.size(); i++) {
        _matrix[i].resize(0, 0);
      }
      _matrix.clear();
      return;
    } // TCMatrix::Cleanup

    /*
     * Modify 3-center matrix elements consistent with use of symmetrized 
     * Coulomb interaction. 
     */
    void TCMatrix_gwbse::MultiplyRightWithAuxMatrix(const Eigen::MatrixXd& matrix) {

#pragma omp parallel for
      for (int i_occ = 0; i_occ < _mtotal; i_occ++) {
    #if (GWBSE_DOUBLE)
          Eigen::MatrixXd temp= _matrix[ i_occ ]*matrix;
       _matrix[ i_occ ] = temp;
#else  
       const Eigen::MatrixXd m = _matrix[ i_occ ].cast<double>();
       _matrix[ i_occ ].noalias()=(m*matrix).cast<float>();
#endif
       
      }
      return;
    } 

    void TCMatrix_gwbse::Print() {

      for (int k = 0; k < _mtotal; k++) {
        std::cout <<k<<std::endl;
        std::cout <<this->_matrix[k]<< std::endl;  
      }
      return;
    }

    /*
     * Fill the 3-center object by looping over shells of GW basis set and
     * calling FillBlock, which calculates all 3-center overlap integrals
     * associated to a particular shell, convoluted with the DFT orbital
     * coefficients
     */
    void TCMatrix_gwbse::Fill(const AOBasis& gwbasis, const AOBasis& dftbasis, const Eigen::MatrixXd& dft_orbitals) {

      // loop over all shells in the GW basis and get _Mmn for that shell
#pragma omp parallel for schedule(guided)//private(_block)
      for (unsigned is = 0; is < gwbasis.getNumofShells(); is++) {
        const AOShell* shell = gwbasis.getShell(is);
        std::vector< Eigen::MatrixXd > block;
        for (int i = 0; i < _mtotal; i++) {
          block.push_back(Eigen::MatrixXd::Zero(_ntotal,shell->getNumFunc()));
        }
        // Fill block for this shell (3-center overlap with _dft_basis + multiplication with _dft_orbitals )
        FillBlock(block, shell, dftbasis, dft_orbitals);

        // put into correct position
        for (int m_level = 0; m_level < this->msize(); m_level++) {
          for (int i_gw = 0; i_gw < shell->getNumFunc(); i_gw++) {
            for (int n_level = 0; n_level < this->nsize(); n_level++) {
              _matrix[m_level]( n_level,shell->getStartIndex() + i_gw) = block[m_level](n_level,i_gw);
            } // n-th DFT orbital
          } // GW basis function in shell
        } // m-th DFT orbital
      } // shells of GW basis set
      return;
    }

    /*
     * Determines the 3-center integrals for a given shell in the GW basis
     * by calculating the 3-center overlap integral of the functions in the
     * GW shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock),
     * followed by a convolution of those with the DFT orbital coefficients 
     */

    void TCMatrix_gwbse::FillBlock(std::vector< Eigen::MatrixXd >& block, const AOShell* auxshell, const AOBasis& dftbasis, const Eigen::MatrixXd& dft_orbitals) {
      tensor3d::extent_gen extents;
      std::vector<Eigen::MatrixXd> symmstorage;
      for (int i = 0; i < auxshell->getNumFunc(); ++i) {
        symmstorage.push_back(Eigen::MatrixXd::Zero(dftbasis.AOBasisSize(), dftbasis.AOBasisSize()));
      }
      const Eigen::MatrixXd dftm = dft_orbitals.block(0, _mmin, dft_orbitals.rows(), _mtotal);
      const Eigen::MatrixXd dftn = dft_orbitals.block(0, _nmin, dft_orbitals.rows(), _ntotal);
      // alpha-loop over the "left" DFT basis function
      for (unsigned row = 0; row < dftbasis.getNumofShells(); row++) {

        const AOShell* shell_row = dftbasis.getShell(row);
        const int row_start = shell_row->getStartIndex();
        // ThreecMatrix is symmetric, restrict explicit calculation to triangular matrix
        for (unsigned col = 0; col <= row; col++) {
          const AOShell* shell_col = dftbasis.getShell(col);
          const int col_start = shell_col->getStartIndex();

          tensor3d threec_block(extents[ range(0, auxshell->getNumFunc()) ][ range(0, shell_row->getNumFunc()) ][ range(0, shell_col->getNumFunc())]);
          for (int i = 0; i < auxshell->getNumFunc(); ++i) {
            for (int j = 0; j < shell_row->getNumFunc(); ++j) {
              for (int k = 0; k < shell_col->getNumFunc(); ++k) {
                threec_block[i][j][k] = 0.0;
              }
            }
          }

          bool nonzero = FillThreeCenterRepBlock(threec_block, auxshell, shell_row, shell_col);
          if (nonzero) {
            for (int _aux = 0; _aux < auxshell->getNumFunc(); _aux++) {
              for (int _row = 0; _row < shell_row->getNumFunc(); _row++) {
                for (int _col = 0; _col < shell_col->getNumFunc(); _col++) {
                  //symmetry
                  if ((col_start + _col)>(row_start + _row)) {
                    continue;
                  }
                  symmstorage[_aux](row_start + _row, col_start + _col) = threec_block[_aux][_row][_col];
                } // ROW copy
              } // COL copy
            } // AUX copy
          }
        } // gamma-loop
      } // alpha-loop
      for (int k = 0; k < auxshell->getNumFunc(); ++k) {
        Eigen::MatrixXd& matrix = symmstorage[k];
        for (int i = 0; i < matrix.rows(); ++i) {
          for (int j = 0; j < i; ++j) {
            matrix(j, i) = matrix(i, j);
          }
        }
        Eigen::MatrixXd threec_inMo = dftn.transpose() * matrix*dftm;
        for (int i = 0; i < threec_inMo.cols(); ++i) {
          for (int j = 0; j < threec_inMo.rows(); ++j) {
            block[i](j, k) = threec_inMo(j, i);
          }
        }
      }
      return;
    }

    void TCMatrix_gwbse::Prune(int min, int max) {
        
      _matrix.resize(max + 1);
      // entries until min can be freed
      for (int i = 0; i < min; i++) {
        _matrix[i].resize(0, 0);
      }
      for (unsigned i = min; i < _matrix.size(); i++) {
        _matrix[i].conservativeResize( max + 1,Eigen::NoChange);
      }
      return;
    }




  }
}

