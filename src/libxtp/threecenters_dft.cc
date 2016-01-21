/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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


using namespace std;
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

        } // TCMatrix_dft::Cleanup

      
  
        void TCMatrix_dft::Fill(AOBasis& _auxbasis, AOBasis& _dftbasis) {

            for (int i=0; i< _auxbasis._AOBasisSize; i++){
                _matrix.push_back(ub::zero_matrix<double>(_dftbasis._AOBasisSize, _dftbasis._AOBasisSize));        
            }
           /*
            std::vector< int > _limits;
            int _temp=0;
            for (vector< AOShell* >::iterator _is = _auxbasis.firstShell(); _is != _auxbasis.lastShell(); _is++) {
                _limits.push_back(_temp);
                 
                _temp+=(*_is)->getNumFunc();
               
                cout << " StartIndex " << (*_is)->getStartIndex() << endl;
            }
            
            for ( int i=0;i<_limits.size();i++){
                cout <<"grenze " << _limits[i]<< endl;
            }
            exit(0);
            */ 
            // loop over all shells in the GW basis and get _Mmn for that shell
            //#pragma omp parallel for //private(_block)
            for ( unsigned _is= 0; _is <  _auxbasis._aoshells.size() ; _is++ ){
            // for (vector< AOShell* >::iterator _is = _gwbasis.firstShell(); _is != _gwbasis.lastShell(); _is++) {
                //cout << " act threads: " << omp_get_thread_num( ) << " total threads " << omp_get_num_threads( ) << " max threads " << omp_get_max_threads( ) <<endl;
                AOShell* _shell = _auxbasis.getShell(_is);
               

                //cout << _is << " shell " << _shell->getType() << " from " << _shell->getStartIndex() << " to " << _shell->getStartIndex() + _shell->getNumFunc() << endl;
                
                
                // Fill block for this shell (3-center overlap with _dft_basis )
                FillBlock(_shell, _dftbasis);
               
                
            } // shells of aux basis set
            //for (int i=0;i<_matrix.size();i++){
                  //  cout<< "Matrix(" << i<<"):\n"<<_matrix[i]<< endl;
                
                //}
        } // TCMatrix_dft::Fill


        
        /*
         * Determines the 3-center integrals for a given shell in the aux basis
         * by calculating the 3-center overlap integral of the functions in the
         * aux shell with ALL functions in the DFT basis set (FillThreeCenterOLBlock)
         */
        
        void TCMatrix_dft::FillBlock(AOShell* _shell, AOBasis& dftbasis) {
	  //void TCMatrix_dft::FillBlock(std::vector< ub::matrix<float> >& _block, AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals) {



           int _start=_shell->getStartIndex();
           //cout << " FB start " << _start << endl;
            // alpha-loop over the "left" DFT basis function
            for (vector< AOShell* >::iterator _row = dftbasis.firstShell(); _row != dftbasis.lastShell(); _row++) {
                AOShell* _shell_row = dftbasis.getShell(_row);
                int _row_start = _shell_row->getStartIndex();
                //int _row_end = _row_start + _shell_row->getNumFunc();

                // gamma-loop over the "right" DFT basis function
                for (vector< AOShell* >::iterator _col = dftbasis.firstShell(); _col != dftbasis.lastShell(); _col++) {
                    AOShell* _shell_col = dftbasis.getShell(_col);
                    int _col_start = _shell_col->getStartIndex();
                    //int _col_end = _col_start + _shell_col->getNumFunc();

                    // get 3-center overlap directly as _subvector
                    ub::matrix<double> _subvector = ub::zero_matrix<double>(_shell_row->getNumFunc(), _shell->getNumFunc() * _shell_col->getNumFunc());
                    //ub::matrix<float> _subvector = ub::zero_matrix<float>(_shell_row->getNumFunc(), _shell->getNumFunc() * _shell_col->getNumFunc());
                    
                    //bool nonzero = FillThreeCenterOLBlock(_subvector, _shell, _shell_row, _shell_col);
                    bool nonzero=FillThreeCenterRepBlock(_subvector, _shell, _shell_row, _shell_col);
                    //cout << "subvector " <<_subvector<< endl;
                    /* test 
                    for (int j=0;j<_subvector.size1();j++){
                    for (int i=0;i<_subvector.size1();i++){
                        
                        _subvector(i,_shell->getNumFunc()*j+i)=1;
                        
                    }}
                    */
                    
                    
                  /*  for (int i = 0; i < _subvector.size1(); i++){
                        for (int j = 0; j < _subvector.size2(); j++){
                            
                            
                            cout << "SV " << _subvector(i,j) << endl;
                            
                        }
                        
                        
                    } */
                    
                    if (nonzero) {
                        // and put it into the block it belongs to
                        // functions in ONE AUXshell
                        for (int _aux = 0; _aux < _shell->getNumFunc(); _aux++) {
                                // column in ONE DFTshell
                                //for (int _i_col = 0; _i_col < _shell_col->getNumFunc(); _i_col++) {
                for (int _col = 0; _col < _shell_col->getNumFunc(); _col++) {
                                    //int _index=_shell_col->getNumFunc() * _aux + _i_col;
                                    int _index=_shell_col->getNumFunc() * _aux + _col;
                                    
                                             //for (int _i_row = 0; _i_row < _shell_row->getNumFunc(); _i_row++) {
                    for (int _row = 0; _row < _shell_row->getNumFunc(); _row++) {
                                                 
                                                //cout << "MAGIC " << _start+_aux << " : " << _row_start + _i_row << " : " << _col_start + _i_col << endl;
                                                _matrix[_start+_aux](_row_start + _row, _col_start + _col) = _subvector(_row, _index);
                                                
                                                } // ROW copy
                                } // COL copy
                        } // AUX copy
                   }
                } // DFT col
            } // DFT row
           
            //       exit(0);
        } // TCMatrix_dft::FillBlock


        
 


    }
}