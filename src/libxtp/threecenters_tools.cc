/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/xtp/threecenters.h>



using namespace votca::tools;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        /*
         * Cleaning TCMatrix data and free memory
         */
        void TCMatrix::Cleanup() {

            for (unsigned _i = 0; _i < _matrix.size(); _i++) {
                _matrix[ _i ].resize(0, 0, false);
            }
            _matrix.clear();
            return;
        } // TCMatrix::Cleanup

        
        /*
         * Modify 3-center matrix elements consistent with use of symmetrized 
         * Coulomb interaction. 
         */
        void TCMatrix::Symmetrize(const ub::matrix<double>& _coulomb) {

            #pragma omp parallel for
            for (int _i_occ = 0; _i_occ < this->get_mtot(); _i_occ++) {
	      // fist cast _matrix[_i_occ] to double for efficient prod() overloading
	      ub::matrix<double> _matrix_double = _matrix[ _i_occ ];
	      _matrix[ _i_occ ] = ub::prod(_coulomb, _matrix_double);
	      
            }
            return;
        } // TCMatrix::Symmetrize
        
           void TCMatrix::Print(string _ident) {
	  //cout << "\n" << endl;
            for (int k = 0; k < this->mtotal; k++) {
                for (unsigned i = 0; i < _matrix[1].size1(); i++) {
                    for (int j = 0; j< this->ntotal; j++) {
                        cout << _ident << "[" << i + 1 << ":" << k + 1 << ":" << j + 1 << "] " << this->_matrix[k](i, j) << endl;
                    }
                }
            }
            return;
        }

        
        /*
         * Fill the 3-center object by looping over shells of GW basis set and
         * calling FillBlock, which calculates all 3-center overlap integrals
         * associated to a particular shell, convoluted with the DFT orbital
         * coefficients
         */
        void TCMatrix::Fill(const AOBasis& _gwbasis,const AOBasis& _dftbasis,const ub::matrix<double>& _dft_orbitals) {

            // loop over all shells in the GW basis and get _Mmn for that shell
            #pragma omp parallel for //private(_block)
            for ( unsigned _is= 0; _is <  _gwbasis.getNumofShells() ; _is++ ){ 
                const AOShell* _shell = _gwbasis.getShell(_is);
                int _start = _shell->getStartIndex();
                // each element is a shell_size-by-n matrix, initialize to zero
                std::vector< ub::matrix<double> > _block(this->get_mtot());
                for (int i = 0; i < this->get_mtot(); i++) {
                    _block[i] = ub::zero_matrix<double>(_shell->getNumFunc(), this->get_ntot());
                }
                // Fill block for this shell (3-center overlap with _dft_basis + multiplication with _dft_orbitals )
                FillBlock(_block, _shell, _dftbasis, _dft_orbitals);

                // put into correct position
                for (int _m_level = 0; _m_level < this->get_mtot(); _m_level++) {
                    for (int _i_gw = 0; _i_gw < _shell->getNumFunc(); _i_gw++) {
                        for (int _n_level = 0; _n_level < this->get_ntot(); _n_level++) {

                            _matrix[_m_level](_start + _i_gw, _n_level) = _block[_m_level](_i_gw, _n_level);

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
        
        void TCMatrix::FillBlock(std::vector< ub::matrix<double> >& _block,const AOShell* _shell,const AOBasis& dftbasis,const ub::matrix<double>& _dft_orbitals) {
	 
            // prepare local storage for 3-center overlap x m-orbitals
            ub::matrix<double> _imstore = ub::zero_matrix<double>(mtotal * _shell->getNumFunc(), dftbasis.AOBasisSize());
        

            // alpha-loop over the "left" DFT basis function
            for (AOBasis::AOShellIterator _row = dftbasis.firstShell(); _row != dftbasis.lastShell(); ++_row) {
                const AOShell* _shell_row = dftbasis.getShell(_row);
                int _row_start = _shell_row->getStartIndex();
                int _row_end = _row_start + _shell_row->getNumFunc();

                // get slice of _dft_orbitals for m-summation, belonging to this shell
                const ub::matrix<double>  _m_orbitals = ub::subrange(_dft_orbitals, mmin, mmax + 1, _row_start, _row_end);

                // gamma-loop over the "right" DFT basis function
                for (AOBasis::AOShellIterator _col = dftbasis.firstShell(); _col != dftbasis.lastShell(); ++_col) {
                    const AOShell* _shell_col = dftbasis.getShell(_col);
                    int _col_start = _shell_col->getStartIndex();
                    int _col_end = _col_start + _shell_col->getNumFunc();

                    // get 3-center overlap directly as _subvector
                    ub::matrix<double> _subvector = ub::zero_matrix<double>(_shell_row->getNumFunc(), _shell->getNumFunc() * _shell_col->getNumFunc());
                    bool nonzero = FillThreeCenterRepBlock(_subvector, _shell, _shell_row, _shell_col);

                    // if this contributes, multiply _subvector with _dft_orbitals and place in _imstore
                    if (nonzero) {

                        ub::matrix<double> _temp = ub::prod(_m_orbitals, _subvector);
                     
                        // put _temp into _imstore
                        for (unsigned _m_level = 0; _m_level < _temp.size1(); _m_level++) {
                            for (int _i_gw = 0; _i_gw < _shell->getNumFunc(); _i_gw++) {
                                int _ridx = _shell->getNumFunc() * (_m_level - this->mmin) + _i_gw;

                                for (int _cidx = _col_start; _cidx < _col_end; _cidx++) {
                                    int _tidx = _shell_col->getNumFunc() * (_i_gw) + _cidx - _col_start;
                                    _imstore(_ridx, _cidx) += _temp(_m_level, _tidx);
                                } // index magic
                            } // GW basis function in shell
                        } // m-level
                    } // IF: adding contribution                        
                } // gamma-loop
            } // alpha-loop


            // get transposed slice of _dft_orbitals
            const ub::matrix<double> _n_orbitals = ub::trans(ub::subrange(_dft_orbitals, nmin, nmax + 1, 0, _dft_orbitals.size2()));

            // Now, finally multiply _imstore with _n_orbitals
            ub::matrix<double> _temp = ub::prod(_imstore, _n_orbitals);
            //ub::matrix<real_gwbse> _temp = ub::prod(_imstore, _n_orbitals);

            // and put it into the block it belongs to
            for (int _m_level = 0; _m_level < mtotal; _m_level++) {
                for (int _i_gw = 0; _i_gw < _shell->getNumFunc(); _i_gw++) {
                    int _midx = _shell->getNumFunc() *(_m_level - mmin) + _i_gw;
                    for (int _n_level = 0; _n_level < ntotal; _n_level++) {
                        _block[_m_level](_i_gw, _n_level) = _temp(_midx, _n_level);
                    } // n-level
                } // GW basis function in shell
            } // m-level
            return;
        } // TCMatrix::FillBlock

        
        void TCMatrix::Prune ( int _basissize, int min, int max){

            int size1 = _matrix[0].size1();           
            // vector needs only max entries
            _matrix.resize( max + 1 );
            // entries until min can be freed
            for ( int i = 0; i < min ; i++){
                _matrix[i].resize(0,0,false);
            }
            for ( unsigned i=min; i < _matrix.size(); i++){
                _matrix[i].resize(size1,max+1);
            }
            return;
        }
        
 


    }
}

