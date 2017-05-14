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



#include <votca/xtp/threecenters.h>



using namespace votca::tools;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        /*
         * Cleaning TCMatrix_dft data and free memory
         */
        void TCMatrix_dft::Cleanup() {

            for (unsigned _i = 0; _i < _matrix.size(); _i++) {
                _matrix[ _i ].resize(0, 0, false);
            }
            _matrix.clear();
            return;
        } // TCMatrix_dft::Cleanup

      
  
        void TCMatrix_dft::Fill(AOBasis& _auxbasis, AOBasis& _dftbasis) {

            for (unsigned int i=0; i< _auxbasis.AOBasisSize(); i++){
                 try{
                _matrix.push_back(ub::symmetric_matrix<double>(_dftbasis.AOBasisSize()));   
                }
                catch(std::bad_alloc& ba){
                    std::cerr << "Basisset/aux basis too large for 3c calculation. Not enough RAM. Caught bad alloc: " << ba.what() << endl;
                    exit(0);
                }
                
              
                for(unsigned k=0;k<_matrix[i].size1();k++){
                    for(unsigned j=0;j<_matrix[i].size2();j++){
                
                    _matrix[i](k,j)=0.0;
                }
                }
            }
            // loop over all shells in the GW basis and get _Mmn for that shell
            #pragma omp parallel for //private(_block)
            for ( unsigned _is= 0; _is <  _auxbasis.getNumofShells() ; _is++ ){
          
                const AOShell* _shell = _auxbasis.getShell(_is);

                // Fill block for this shell (3-center overlap with _dft_basis )
                FillBlock(_shell, _dftbasis);
               
                
            } // shells of aux basis set
            return;
        } // TCMatrix_dft::Fill


        
        /*
         * Determines the 3-center integrals for a given shell in the aux basis
         * by calculating the 3-center overlap integral of the functions in the
         * aux shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock)
         */
        
        void TCMatrix_dft::FillBlock(const AOShell* _shell, const AOBasis& dftbasis) {


            int _start = _shell->getStartIndex();

            // alpha-loop over the "left" DFT basis function
            for (AOBasis::AOShellIterator _row = dftbasis.firstShell(); _row != dftbasis.lastShell(); ++_row) {
                const AOShell* _shell_row = dftbasis.getShell(_row);
                int _row_start = _shell_row->getStartIndex();

                // gamma-loop over the "right" DFT basis function with symmetry
                for (AOBasis::AOShellIterator _col = dftbasis.firstShell(); _col <= _row; ++_col) {
                    const AOShell* _shell_col = dftbasis.getShell(_col);
                    int _col_start = _shell_col->getStartIndex();

                    // get 3-center overlap directly as _subvector
                    ub::matrix<double> _subvector = ub::zero_matrix<double>(_shell_row->getNumFunc(), _shell->getNumFunc() * _shell_col->getNumFunc());

                    bool nonzero = FillThreeCenterRepBlock(_subvector, _shell, _shell_row, _shell_col);



                    if (nonzero) {
                        // and put it into the block it belongs to
                        // functions in ONE AUXshell
                        for (int _aux = 0; _aux < _shell->getNumFunc(); _aux++) {
                            // column in ONE DFTshell
                            for (int _col = 0; _col < _shell_col->getNumFunc(); _col++) {
                                int _index = _shell_col->getNumFunc() * _aux + _col;

                                for (int _row = 0; _row < _shell_row->getNumFunc(); _row++) {
                                    //symmetry
                                    if ((_col_start + _col)>(_row_start + _row)) {
                                        continue;
                                    }
                                   
                                    _matrix[_start + _aux](_row_start + _row, _col_start + _col) = _subvector(_row, _index);

                                } // ROW copy
                            } // COL copy
                        } // AUX copy
                    }
                } // DFT col
            } // DFT row
            return;
        } // TCMatrix_dft::FillBlock


        
 


    }
}
