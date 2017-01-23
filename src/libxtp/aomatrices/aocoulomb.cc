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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/linalg.h>


using namespace votca::tools;

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;

    
    int AOCoulomb::getExtraBlockSize(int _lmax_row, int _lmax_col){
        int _block_size = _lmax_col + _lmax_row +1;
        return _block_size;
    }

    void AOCoulomb::FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, AOBasis* ecp) {

            // shell info, only lmax tells how far to go
            const int _lmax_row = _shell_row->getLmax();
            const int _lmax_col = _shell_col->getLmax();

            // set size of internal block for recursion
            int _nrows = this->getBlockSize(_lmax_row);
            int _ncols = this->getBlockSize(_lmax_col);
            int _nextra = this->getExtraBlockSize(_lmax_row, _lmax_col);
            int _l_sum = _lmax_row + _lmax_col;
           // int _ma_dim = this->getBlockSize(_l_sum);
            
            int nmax=20; // This is hardcoded using getBlocksize leads to problems the if clauses are not that restrictive and so if you do use a smaller array it might lead to problems
             if(_lmax_row>3 ||_lmax_col>3){
                 nmax=35;
             }
            
            // get shell positions
            const vec& _pos_row = _shell_row->getPos();
            const vec& _pos_col = _shell_col->getPos();
            const vec _diff = _pos_row - _pos_col;
            double _distsq = (_diff.getX() * _diff.getX()) + (_diff.getY() * _diff.getY()) + (_diff.getZ() * _diff.getZ());
            
            const double pi = boost::math::constants::pi<double>();
             // some helpers
            std::vector<double> _wmp;
            std::vector<double> _wmq;
            _wmp.resize(3);
            _wmq.resize(3);
            
         

           
            
            typedef std::vector< AOGaussianPrimitive* >::iterator GaussianIterator;
        // iterate over Gaussians in this _shell_row
            for ( GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
                const double& _decay_row = (*itr)->decay;
            
                for ( GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
                    
                     // get decay constants 
                        const double& _decay_col = (*itc)->decay;

                       
                        // check if distance between postions is big, then skip step   
                        //double _exparg = _fakac2 * _decay_row * _decay_col *_distsq;
                        // if ( _exparg > 30.0 ) { continue; } //!!!!!CUTOFF not applicable to AOCoulomb (at least not like this...)
                    
                                    // get a multi dimensional array
                         
                         //ma_type _cou(boost::extents[_nrows][_ncols][_nextra]);
                         ma_type _cou(boost::extents[nmax][nmax][_nextra]);
                         
                                    // initialize to zero_cou[0][0][i] 
                           //       for(index i = 0; i != _nrows; ++i) {
                           //  for(index j = 0; j != _ncols; ++j){
                           for (index i = 0; i != nmax; ++i) {
                               for (index j = 0; j != nmax; ++j) {
                                   for (index k = 0; k != _nextra; ++k) {
                                       _cou[i][j][k] = 0.0;
                                   }
                               }
                           }

                       

            // some helpers
            const double _fakac = 0.5 / (_decay_row + _decay_col);
            const double _fakac2 = 1. / (_decay_row + _decay_col);
            const double _faka = 0.5 / _decay_row;
            const double _faka2 = 2.0 * _faka;

            const double _fakaca = _decay_row / (_decay_row + _decay_col);
            const double _fakaac = _decay_row / (_decay_row + _decay_col);
            
            const double _fakac3 = 3.0 * _fakac;
            //const double _fakac4 = 4.0 * _fakac;
            const double _fakc = 0.5 / _decay_col;
            const double _fakc2 = 2.0 * _fakc;


            
            //bool ident;
            
            //if ( )
            
            if (sqrt(_distsq) > 0.01 ){
            _wmp[0] = _fakac2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_row.getX();
            _wmp[1] = _fakac2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_row.getY();
            _wmp[2] = _fakac2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_row.getZ();

            _wmq[0] = _fakac2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_col.getX();
            _wmq[1] = _fakac2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_col.getY();
            _wmq[2] = _fakac2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_col.getZ();

            } else {
            _wmp[0] = 0.0;
            _wmp[1] = 0.0;
            _wmp[2] = 0.0;

            _wmq[0] = 0.0;
            _wmq[1] = 0.0;
            _wmq[2] = 0.0;
                

            }
            const double _T = _fakaca * _decay_col * _distsq;

        

            

            double _fak = 2.0 * pow(pi, 2.5) / (_decay_row * _decay_col * sqrt(_decay_row + _decay_col));
            _fak = _fak * pow(4.0 * _decay_row * _decay_col / (pi * pi), 0.75);

            std::vector<double> _FmT(_nextra, 0.0); // that size needs to be checked!
            // call xint01(FmT,8,T,u_lower)
            XIntegrate(_FmT, _T);

            // get initial data from _FmT -> s-s element
            for (index i = 0; i != _nextra; ++i) {
                _cou[0][0][i] = _fak * _FmT[i];
            }


            if (_l_sum >= 1) {
                // p-s || s-p

                // p-s elements
                if (_lmax_row > 0) {
                    _cou[1][0][0] = _wmp[0] * _cou[0][0][1]; // !p-s
                    _cou[2][0][0] = _wmp[1] * _cou[0][0][1]; // !p-s
                    _cou[3][0][0] = _wmp[2] * _cou[0][0][1]; // !p-s  
                }

                // s-p elements
                if (_lmax_col > 0) {
                    _cou[0][1][0] = _wmq[0] * _cou[0][0][1]; //!s-p
                    _cou[0][2][0] = _wmq[1] * _cou[0][0][1]; //!s-p
                    _cou[0][3][0] = _wmq[2] * _cou[0][0][1]; //!s-p    
                }
            }

            if (_l_sum >= 2) {
                // d-s || s-d || p-p

                // p-s-1, needed in d-s, d-p, f-s (_lmax_row > 1 )
                _cou[1][0][1] = _wmp[0] * _cou[0][0][2];
                _cou[2][0][1] = _wmp[1] * _cou[0][0][2];
                _cou[3][0][1] = _wmp[2] * _cou[0][0][2];

                // s-p-1, needed in s-d, p-p, p-d, d-p ( _lmax_col > 0 )
                _cou[0][1][1] = _wmq[0] * _cou[0][0][2];
                _cou[0][2][1] = _wmq[1] * _cou[0][0][2];
                _cou[0][3][1] = _wmq[2] * _cou[0][0][2];

                // d-s elements, req.: p-s(1)
                if (_lmax_row > 1) {
                    _cou[4][0][0] = _wmp[1] * _cou[1][0][1];
                    _cou[5][0][0] = _wmp[2] * _cou[1][0][1];
                    _cou[7][0][0] = _wmp[0] * _cou[1][0][1] + _faka * (_cou[0][0][0] - _fakaac * _cou[0][0][1]);
                    _cou[8][0][0] = _wmp[1] * _cou[2][0][1] + _faka * (_cou[0][0][0] - _fakaac * _cou[0][0][1]);
                    _cou[6][0][0] = _wmp[2] * _cou[2][0][1];
                    _cou[9][0][0] = _wmp[2] * _cou[3][0][1] + _faka * (_cou[0][0][0] - _fakaac * _cou[0][0][1]);
                }

                // s-d elements req.: s-p(1)
                if (_lmax_col > 1) {


                    _cou[0][7][0] = _wmq[0] * _cou[0][1][1] + _fakc * (_cou[0][0][0] - _fakaca * _cou[0][0][1]);
                    _cou[0][4][0] = _wmq[1] * _cou[0][1][1];
                    _cou[0][5][0] = _wmq[2] * _cou[0][1][1];
                    _cou[0][8][0] = _wmq[1] * _cou[0][2][1] + _fakc * (_cou[0][0][0] - _fakaca * _cou[0][0][1]);
                    _cou[0][6][0] = _wmq[2] * _cou[0][2][1];
                    _cou[0][9][0] = _wmq[2] * _cou[0][3][1] + _fakc * (_cou[0][0][0] - _fakaca * _cou[0][0][1]);
                }

                // p-p elements req.: s-p(1)
                if (_lmax_row > 0 && _lmax_col > 0) {

                    //cout << "Setting p-p elements" << endl;
                    _cou[1][1][0] = _wmp[0] * _cou[0][1][1] + _fakac * _cou[0][0][1];
                    _cou[2][1][0] = _wmp[1] * _cou[0][1][1];
                    _cou[3][1][0] = _wmp[2] * _cou[0][1][1];
                    _cou[1][2][0] = _wmp[0] * _cou[0][2][1];
                    _cou[2][2][0] = _wmp[1] * _cou[0][2][1] + _fakac * _cou[0][0][1];
                    _cou[3][2][0] = _wmp[2] * _cou[0][2][1];
                    _cou[1][3][0] = _wmp[0] * _cou[0][3][1];
                    _cou[2][3][0] = _wmp[1] * _cou[0][3][1];
                    _cou[3][3][0] = _wmp[2] * _cou[0][3][1] + _fakac * _cou[0][0][1];
                }
            }



            if (_l_sum >= 3) {

                // p-s-2 -> d-s-1 | d-p-1 -> f-s, f-p, f-d ( _lmax_row > 2 )
                _cou[7][0][1] = _wmp[0] * _cou[1][0][2] + _faka * (_cou[0][0][1] - _fakaac * _cou[0][0][2]);
                _cou[4][0][1] = _wmp[1] * _cou[1][0][2];
                _cou[5][0][1] = _wmp[2] * _cou[1][0][2];
                _cou[8][0][1] = _wmp[1] * _cou[2][0][2] + _faka * (_cou[0][0][1] - _fakaac * _cou[0][0][2]);
                _cou[6][0][1] = _wmp[2] * _cou[2][0][2];
                _cou[9][0][1] = _wmp[2] * _cou[3][0][2] + _faka * (_cou[0][0][1] - _fakaac * _cou[0][0][2]);

                // s-p-2 -> p-d, d-d, d-p 
                _cou[0][1][2] = _wmq[0] * _cou[0][0][3];
                _cou[0][2][2] = _wmq[1] * _cou[0][0][3];
                _cou[0][3][2] = _wmq[2] * _cou[0][0][3];

                // d-s-1 (only g )

                // p-p-1
                _cou[1][1][1] = _wmp[0] * _cou[0][1][2] + _fakac * _cou[0][0][2];
                _cou[2][1][1] = _wmp[1] * _cou[0][1][2];
                _cou[3][1][1] = _wmp[2] * _cou[0][1][2];
                _cou[1][2][1] = _wmp[0] * _cou[0][2][2];
                _cou[2][2][1] = _wmp[1] * _cou[0][2][2] + _fakac * _cou[0][0][2];
                _cou[3][2][1] = _wmp[2] * _cou[0][2][2];
                _cou[1][3][1] = _wmp[0] * _cou[0][3][2];
                _cou[2][3][1] = _wmp[1] * _cou[0][3][2];
                _cou[3][3][1] = _wmp[2] * _cou[0][3][2] + _fakac * _cou[0][0][2];

                // s-d-1
                _cou[0][7][1] = _wmq[0] * _cou[0][1][2] + _fakc * (_cou[0][0][1] - _fakaca * _cou[0][0][2]);
                _cou[0][4][1] = _wmq[1] * _cou[0][1][2];
                _cou[0][5][1] = _wmq[2] * _cou[0][1][2];
                _cou[0][8][1] = _wmq[1] * _cou[0][2][2] + _fakc * (_cou[0][0][1] - _fakaca * _cou[0][0][2]);
                _cou[0][6][1] = _wmq[2] * _cou[0][2][2];
                _cou[0][9][1] = _wmq[2] * _cou[0][3][2] + _fakc * (_cou[0][0][1] - _fakaca * _cou[0][0][2]);

                // f-s
                if (_lmax_row > 2) {
                    _cou[13][0][0] = _wmp[0] * _cou[4][0][1] + _faka * (_cou[2][0][0] - _fakaac * _cou[2][0][1]);
                    _cou[14][0][0] = _wmp[1] * _cou[4][0][1] + _faka * (_cou[1][0][0] - _fakaac * _cou[1][0][1]);
                    _cou[19][0][0] = _wmp[2] * _cou[4][0][1];
                    _cou[15][0][0] = _wmp[0] * _cou[5][0][1] + _faka * (_cou[3][0][0] - _fakaac * _cou[3][0][1]);
                    _cou[16][0][0] = _wmp[2] * _cou[5][0][1] + _faka * (_cou[1][0][0] - _fakaac * _cou[1][0][1]);
                    _cou[17][0][0] = _wmp[1] * _cou[6][0][1] + _faka * (_cou[3][0][0] - _fakaac * _cou[3][0][1]);
                    _cou[18][0][0] = _wmp[2] * _cou[6][0][1] + _faka * (_cou[2][0][0] - _fakaac * _cou[2][0][1]);
                    _cou[10][0][0] = _wmp[0] * _cou[7][0][1] + _faka2 * (_cou[1][0][0] - _fakaac * _cou[1][0][1]);
                    _cou[11][0][0] = _wmp[1] * _cou[8][0][1] + _faka2 * (_cou[2][0][0] - _fakaac * _cou[2][0][1]);
                    _cou[12][0][0] = _wmp[2] * _cou[9][0][1] + _faka2 * (_cou[3][0][0] - _fakaac * _cou[3][0][1]);
                }

                // d-p
                if (_lmax_row > 1 && _lmax_col > 0) {
                    _cou[4][1][0] = _wmp[1] * _cou[1][1][1];
                    _cou[4][2][0] = _wmp[1] * _cou[1][2][1] + _fakac * _cou[1][0][1];
                    _cou[4][3][0] = _wmp[1] * _cou[1][3][1];

                    _cou[5][1][0] = _wmp[2] * _cou[1][1][1];
                    _cou[5][2][0] = _wmp[2] * _cou[1][2][1];
                    _cou[5][3][0] = _wmp[2] * _cou[1][3][1] + _fakac * _cou[1][0][1];

                    _cou[6][1][0] = _wmp[2] * _cou[2][1][1];
                    _cou[6][2][0] = _wmp[2] * _cou[2][2][1];
                    _cou[6][3][0] = _wmp[2] * _cou[2][3][1] + _fakac * _cou[2][0][1];

                    _cou[7][1][0] = _wmp[0] * _cou[1][1][1] + _faka * (_cou[0][1][0] - _fakaac * _cou[0][1][1]) + _fakac * _cou[1][0][1];
                    _cou[7][2][0] = _wmp[0] * _cou[1][2][1] + _faka * (_cou[0][2][0] - _fakaac * _cou[0][2][1]);
                    _cou[7][3][0] = _wmp[0] * _cou[1][3][1] + _faka * (_cou[0][3][0] - _fakaac * _cou[0][3][1]);

                    _cou[8][1][0] = _wmp[1] * _cou[2][1][1] + _faka * (_cou[0][1][0] - _fakaac * _cou[0][1][1]);
                    _cou[8][2][0] = _wmp[1] * _cou[2][2][1] + _faka * (_cou[0][2][0] - _fakaac * _cou[0][2][1]) + _fakac * _cou[2][0][1];
                    _cou[8][3][0] = _wmp[1] * _cou[2][3][1] + _faka * (_cou[0][3][0] - _fakaac * _cou[0][3][1]);

                    _cou[9][1][0] = _wmp[2] * _cou[3][1][1] + _faka * (_cou[0][1][0] - _fakaac * _cou[0][1][1]);
                    _cou[9][2][0] = _wmp[2] * _cou[3][2][1] + _faka * (_cou[0][2][0] - _fakaac * _cou[0][2][1]);
                    _cou[9][3][0] = _wmp[2] * _cou[3][3][1] + _faka * (_cou[0][3][0] - _fakaac * _cou[0][3][1]) + _fakac * _cou[3][0][1];
                }

                // p-d
                if (_lmax_row > 0 && _lmax_col > 1) {
                    _cou[1][4][0] = _wmp[0] * _cou[0][4][1] + _fakac * _cou[0][2][1];
                    _cou[2][4][0] = _wmp[1] * _cou[0][4][1] + _fakac * _cou[0][1][1];
                    _cou[3][4][0] = _wmp[2] * _cou[0][4][1];
                    _cou[1][5][0] = _wmp[0] * _cou[0][5][1] + _fakac * _cou[0][3][1];
                    _cou[2][5][0] = _wmp[1] * _cou[0][5][1];
                    _cou[3][5][0] = _wmp[2] * _cou[0][5][1] + _fakac * _cou[0][1][1];
                    _cou[1][6][0] = _wmp[0] * _cou[0][6][1];
                    _cou[2][6][0] = _wmp[1] * _cou[0][6][1] + _fakac * _cou[0][3][1];
                    _cou[3][6][0] = _wmp[2] * _cou[0][6][1] + _fakac * _cou[0][2][1];
                    _cou[1][7][0] = _wmp[0] * _cou[0][7][1] + _fakac2 * _cou[0][1][1];
                    _cou[2][7][0] = _wmp[1] * _cou[0][7][1];
                    _cou[3][7][0] = _wmp[2] * _cou[0][7][1];
                    _cou[1][8][0] = _wmp[0] * _cou[0][8][1];
                    _cou[2][8][0] = _wmp[1] * _cou[0][8][1] + _fakac2 * _cou[0][2][1];
                    _cou[3][8][0] = _wmp[2] * _cou[0][8][1];
                    _cou[1][9][0] = _wmp[0] * _cou[0][9][1];
                    _cou[2][9][0] = _wmp[1] * _cou[0][9][1];
                    _cou[3][9][0] = _wmp[2] * _cou[0][9][1] + _fakac2 * _cou[0][3][1];
                }

                // s-f
                if (_lmax_col > 2) {
                    _cou[0][10][0] = _wmq[0] * _cou[0][7][1] + _fakc2 * (_cou[0][1][0] - _fakaca * _cou[0][1][1]);
                    _cou[0][11][0] = _wmq[1] * _cou[0][8][1] + _fakc2 * (_cou[0][2][0] - _fakaca * _cou[0][2][1]);
                    _cou[0][12][0] = _wmq[2] * _cou[0][9][1] + _fakc2 * (_cou[0][3][0] - _fakaca * _cou[0][3][1]);
                    _cou[0][13][0] = _wmq[0] * _cou[0][4][1] + _fakc * (_cou[0][2][0] - _fakaca * _cou[0][2][1]);
                    _cou[0][14][0] = _wmq[1] * _cou[0][4][1] + _fakc * (_cou[0][1][0] - _fakaca * _cou[0][1][1]);
                    _cou[0][19][0] = _wmq[2] * _cou[0][4][1];
                    _cou[0][15][0] = _wmq[0] * _cou[0][5][1] + _fakc * (_cou[0][3][0] - _fakaca * _cou[0][3][1]);
                    _cou[0][16][0] = _wmq[2] * _cou[0][5][1] + _fakc * (_cou[0][1][0] - _fakaca * _cou[0][1][1]);
                    _cou[0][17][0] = _wmq[1] * _cou[0][6][1] + _fakc * (_cou[0][3][0] - _fakaca * _cou[0][3][1]);
                    _cou[0][18][0] = _wmq[2] * _cou[0][6][1] + _fakc * (_cou[0][2][0] - _fakaca * _cou[0][2][1]);
                }

            }

            if (_l_sum >= 4) {

                // p-s-3 (for g)

                // s-p-3
                _cou[0][1][3] = _wmq[0] * _cou[0][0][4];
                _cou[0][2][3] = _wmq[1] * _cou[0][0][4];
                _cou[0][3][3] = _wmq[2] * _cou[0][0][4];

                // d-s-2 (for g)

                // p-p-2
                _cou[1][1][2] = _wmp[0] * _cou[0][1][3] + _fakac * _cou[0][0][3];
                _cou[2][1][2] = _wmp[1] * _cou[0][1][3];
                _cou[3][1][2] = _wmp[2] * _cou[0][1][3];
                _cou[1][2][2] = _wmp[0] * _cou[0][2][3];
                _cou[2][2][2] = _wmp[1] * _cou[0][2][3] + _fakac * _cou[0][0][3];
                _cou[3][2][2] = _wmp[2] * _cou[0][2][3];
                _cou[1][3][2] = _wmp[0] * _cou[0][3][3];
                _cou[2][3][2] = _wmp[1] * _cou[0][3][3];
                _cou[3][3][2] = _wmp[2] * _cou[0][3][3] + _fakac * _cou[0][0][3];

                // s-d-2
                _cou[0][7][2] = _wmq[0] * _cou[0][1][3] + _fakc * (_cou[0][0][2] - _fakaca * _cou[0][0][3]);
                _cou[0][4][2] = _wmq[1] * _cou[0][1][3];
                _cou[0][5][2] = _wmq[2] * _cou[0][1][3];
                _cou[0][8][2] = _wmq[1] * _cou[0][2][3] + _fakc * (_cou[0][0][2] - _fakaca * _cou[0][0][3]);
                _cou[0][6][2] = _wmq[2] * _cou[0][2][3];
                _cou[0][9][2] = _wmq[2] * _cou[0][3][3] + _fakc * (_cou[0][0][2] - _fakaca * _cou[0][0][3]);

                // f-s-1 (only g))

                // d-p-1
                _cou[7][1][1] = _wmp[0] * _cou[1][1][2] + _faka * (_cou[0][1][1] - _fakaac * _cou[0][1][2]) + _fakac * _cou[1][0][2];
                _cou[4][1][1] = _wmp[1] * _cou[1][1][2];
                _cou[5][1][1] = _wmp[2] * _cou[1][1][2];
                _cou[7][2][1] = _wmp[0] * _cou[1][2][2] + _faka * (_cou[0][2][1] - _fakaac * _cou[0][2][2]);
                _cou[4][2][1] = _wmp[1] * _cou[1][2][2] + _fakac * _cou[1][0][2];
                _cou[5][2][1] = _wmp[2] * _cou[1][2][2];
                _cou[7][3][1] = _wmp[0] * _cou[1][3][2] + _faka * (_cou[0][3][1] - _fakaac * _cou[0][3][2]);
                _cou[4][3][1] = _wmp[1] * _cou[1][3][2];
                _cou[5][3][1] = _wmp[2] * _cou[1][3][2] + _fakac * _cou[1][0][2];
                _cou[8][1][1] = _wmp[1] * _cou[2][1][2] + _faka * (_cou[0][1][1] - _fakaac * _cou[0][1][2]);
                _cou[6][1][1] = _wmp[2] * _cou[2][1][2];
                _cou[8][2][1] = _wmp[1] * _cou[2][2][2] + _faka * (_cou[0][2][1] - _fakaac * _cou[0][2][2]) + _fakac * _cou[2][0][2];
                _cou[6][2][1] = _wmp[2] * _cou[2][2][2];
                _cou[8][3][1] = _wmp[1] * _cou[2][3][2] + _faka * (_cou[0][3][1] - _fakaac * _cou[0][3][2]);
                _cou[6][3][1] = _wmp[2] * _cou[2][3][2] + _fakac * _cou[2][0][2];
                _cou[9][1][1] = _wmp[2] * _cou[3][1][2] + _faka * (_cou[0][1][1] - _fakaac * _cou[0][1][2]);
                _cou[9][2][1] = _wmp[2] * _cou[3][2][2] + _faka * (_cou[0][2][1] - _fakaac * _cou[0][2][2]);
                _cou[9][3][1] = _wmp[2] * _cou[3][3][2] + _faka * (_cou[0][3][1] - _fakaac * _cou[0][3][2]) + _fakac * _cou[3][0][2];

                // p-d-1
                _cou[1][4][1] = _wmp[0] * _cou[0][4][2] + _fakac * _cou[0][2][2];
                _cou[2][4][1] = _wmp[1] * _cou[0][4][2] + _fakac * _cou[0][1][2];
                _cou[3][4][1] = _wmp[2] * _cou[0][4][2];
                _cou[1][5][1] = _wmp[0] * _cou[0][5][2] + _fakac * _cou[0][3][2];
                _cou[2][5][1] = _wmp[1] * _cou[0][5][2];
                _cou[3][5][1] = _wmp[2] * _cou[0][5][2] + _fakac * _cou[0][1][2];
                _cou[1][6][1] = _wmp[0] * _cou[0][6][2];
                _cou[2][6][1] = _wmp[1] * _cou[0][6][2] + _fakac * _cou[0][3][2];
                _cou[3][6][1] = _wmp[2] * _cou[0][6][2] + _fakac * _cou[0][2][2];
                _cou[1][7][1] = _wmp[0] * _cou[0][7][2] + _fakac2 * _cou[0][1][2];
                _cou[2][7][1] = _wmp[1] * _cou[0][7][2];
                _cou[3][7][1] = _wmp[2] * _cou[0][7][2];
                _cou[1][8][1] = _wmp[0] * _cou[0][8][2];
                _cou[2][8][1] = _wmp[1] * _cou[0][8][2] + _fakac2 * _cou[0][2][2];
                _cou[3][8][1] = _wmp[2] * _cou[0][8][2];
                _cou[1][9][1] = _wmp[0] * _cou[0][9][2];
                _cou[2][9][1] = _wmp[1] * _cou[0][9][2];
                _cou[3][9][1] = _wmp[2] * _cou[0][9][2] + _fakac2 * _cou[0][3][2];

                // s-f-1
                _cou[0][10][1] = _wmq[0] * _cou[0][7][2] + _fakc2 * (_cou[0][1][1] - _fakaca * _cou[0][1][2]);
                _cou[0][11][1] = _wmq[1] * _cou[0][8][2] + _fakc2 * (_cou[0][2][1] - _fakaca * _cou[0][2][2]);
                _cou[0][12][1] = _wmq[2] * _cou[0][9][2] + _fakc2 * (_cou[0][3][1] - _fakaca * _cou[0][3][2]);
                _cou[0][13][1] = _wmq[0] * _cou[0][4][2] + _fakc * (_cou[0][2][1] - _fakaca * _cou[0][2][2]);
                _cou[0][14][1] = _wmq[1] * _cou[0][4][2] + _fakc * (_cou[0][1][1] - _fakaca * _cou[0][1][2]);
                _cou[0][19][1] = _wmq[2] * _cou[0][4][2];
                _cou[0][15][1] = _wmq[0] * _cou[0][5][2] + _fakc * (_cou[0][3][1] - _fakaca * _cou[0][3][2]);
                _cou[0][16][1] = _wmq[2] * _cou[0][5][2] + _fakc * (_cou[0][1][1] - _fakaca * _cou[0][1][2]);
                _cou[0][17][1] = _wmq[1] * _cou[0][6][2] + _fakc * (_cou[0][3][1] - _fakaca * _cou[0][3][2]);
                _cou[0][18][1] = _wmq[2] * _cou[0][6][2] + _fakc * (_cou[0][2][1] - _fakaca * _cou[0][2][2]);

                // f-p
                if (_lmax_row > 2) {
                    _cou[10][1][0] = _wmp[0] * _cou[7][1][1] + _faka2 * (_cou[1][1][0] - _fakaac * _cou[1][1][1]) + _fakac * _cou[7][0][1];
                    _cou[10][2][0] = _wmp[0] * _cou[7][2][1] + _faka2 * (_cou[1][2][0] - _fakaac * _cou[1][2][1]);
                    _cou[10][3][0] = _wmp[0] * _cou[7][3][1] + _faka2 * (_cou[1][3][0] - _fakaac * _cou[1][3][1]);
                    _cou[11][1][0] = _wmp[1] * _cou[8][1][1] + _faka2 * (_cou[2][1][0] - _fakaac * _cou[2][1][1]);
                    _cou[11][2][0] = _wmp[1] * _cou[8][2][1] + _faka2 * (_cou[2][2][0] - _fakaac * _cou[2][2][1]) + _fakac * _cou[8][0][1];
                    _cou[11][3][0] = _wmp[1] * _cou[8][3][1] + _faka2 * (_cou[2][3][0] - _fakaac * _cou[2][3][1]);
                    _cou[12][1][0] = _wmp[2] * _cou[9][1][1] + _faka2 * (_cou[3][1][0] - _fakaac * _cou[3][1][1]);
                    _cou[12][2][0] = _wmp[2] * _cou[9][2][1] + _faka2 * (_cou[3][2][0] - _fakaac * _cou[3][2][1]);
                    _cou[12][3][0] = _wmp[2] * _cou[9][3][1] + _faka2 * (_cou[3][3][0] - _fakaac * _cou[3][3][1]) + _fakac * _cou[9][0][1];
                    _cou[13][1][0] = _wmp[0] * _cou[4][1][1] + _faka * (_cou[2][1][0] - _fakaac * _cou[2][1][1]) + _fakac * _cou[4][0][1];
                    _cou[13][2][0] = _wmp[0] * _cou[4][2][1] + _faka * (_cou[2][2][0] - _fakaac * _cou[2][2][1]);
                    _cou[13][3][0] = _wmp[0] * _cou[4][3][1] + _faka * (_cou[2][3][0] - _fakaac * _cou[2][3][1]);
                    _cou[14][1][0] = _wmp[1] * _cou[4][1][1] + _faka * (_cou[1][1][0] - _fakaac * _cou[1][1][1]);
                    _cou[14][2][0] = _wmp[1] * _cou[4][2][1] + _faka * (_cou[1][2][0] - _fakaac * _cou[1][2][1]) + _fakac * _cou[4][0][1];
                    _cou[14][3][0] = _wmp[1] * _cou[4][3][1] + _faka * (_cou[1][3][0] - _fakaac * _cou[1][3][1]);
                    _cou[15][1][0] = _wmp[0] * _cou[5][1][1] + _faka * (_cou[3][1][0] - _fakaac * _cou[3][1][1]) + _fakac * _cou[5][0][1];
                    _cou[15][2][0] = _wmp[0] * _cou[5][2][1] + _faka * (_cou[3][2][0] - _fakaac * _cou[3][2][1]);
                    _cou[15][3][0] = _wmp[0] * _cou[5][3][1] + _faka * (_cou[3][3][0] - _fakaac * _cou[3][3][1]);
                    _cou[16][1][0] = _wmp[2] * _cou[5][1][1] + _faka * (_cou[1][1][0] - _fakaac * _cou[1][1][1]);
                    _cou[16][2][0] = _wmp[2] * _cou[5][2][1] + _faka * (_cou[1][2][0] - _fakaac * _cou[1][2][1]);
                    _cou[16][3][0] = _wmp[2] * _cou[5][3][1] + _faka * (_cou[1][3][0] - _fakaac * _cou[1][3][1]) + _fakac * _cou[5][0][1];
                    _cou[17][1][0] = _wmp[1] * _cou[6][1][1] + _faka * (_cou[3][1][0] - _fakaac * _cou[3][1][1]);
                    _cou[17][2][0] = _wmp[1] * _cou[6][2][1] + _faka * (_cou[3][2][0] - _fakaac * _cou[3][2][1]) + _fakac * _cou[6][0][1];
                    _cou[17][3][0] = _wmp[1] * _cou[6][3][1] + _faka * (_cou[3][3][0] - _fakaac * _cou[3][3][1]);
                    _cou[18][1][0] = _wmp[2] * _cou[6][1][1] + _faka * (_cou[2][1][0] - _fakaac * _cou[2][1][1]);
                    _cou[18][2][0] = _wmp[2] * _cou[6][2][1] + _faka * (_cou[2][2][0] - _fakaac * _cou[2][2][1]);
                    _cou[18][3][0] = _wmp[2] * _cou[6][3][1] + _faka * (_cou[2][3][0] - _fakaac * _cou[2][3][1]) + _fakac * _cou[6][0][1];
                    _cou[19][1][0] = _wmp[2] * _cou[4][1][1];
                    _cou[19][2][0] = _wmp[2] * _cou[4][2][1];
                    _cou[19][3][0] = _wmp[2] * _cou[4][3][1] + _fakac * _cou[4][0][1];
                }

                // d-d
                if (_lmax_row > 1 && _lmax_col > 1) {
                    _cou[7][4][0] = _wmp[0] * _cou[1][4][1] + _faka * (_cou[0][4][0] - _fakaac * _cou[0][4][1]) + _fakac * _cou[1][2][1];
                    _cou[4][4][0] = _wmp[1] * _cou[1][4][1] + _fakac * _cou[1][1][1];
                    _cou[5][4][0] = _wmp[2] * _cou[1][4][1];
                    _cou[7][5][0] = _wmp[0] * _cou[1][5][1] + _faka * (_cou[0][5][0] - _fakaac * _cou[0][5][1]) + _fakac * _cou[1][3][1];
                    _cou[4][5][0] = _wmp[1] * _cou[1][5][1];
                    _cou[5][5][0] = _wmp[2] * _cou[1][5][1] + _fakac * _cou[1][1][1];
                    _cou[7][6][0] = _wmp[0] * _cou[1][6][1] + _faka * (_cou[0][6][0] - _fakaac * _cou[0][6][1]);
                    _cou[4][6][0] = _wmp[1] * _cou[1][6][1] + _fakac * _cou[1][3][1];
                    _cou[5][6][0] = _wmp[2] * _cou[1][6][1] + _fakac * _cou[1][2][1];
                    _cou[7][7][0] = _wmp[0] * _cou[1][7][1] + _faka * (_cou[0][7][0] - _fakaac * _cou[0][7][1]) + _fakac2 * _cou[1][1][1];
                    _cou[4][7][0] = _wmp[1] * _cou[1][7][1];
                    _cou[5][7][0] = _wmp[2] * _cou[1][7][1];
                    _cou[7][8][0] = _wmp[0] * _cou[1][8][1] + _faka * (_cou[0][8][0] - _fakaac * _cou[0][8][1]);
                    _cou[4][8][0] = _wmp[1] * _cou[1][8][1] + _fakac2 * _cou[1][2][1];
                    _cou[5][8][0] = _wmp[2] * _cou[1][8][1];
                    _cou[7][9][0] = _wmp[0] * _cou[1][9][1] + _faka * (_cou[0][9][0] - _fakaac * _cou[0][9][1]);
                    _cou[4][9][0] = _wmp[1] * _cou[1][9][1];
                    _cou[5][9][0] = _wmp[2] * _cou[1][9][1] + _fakac2 * _cou[1][3][1];
                    _cou[8][4][0] = _wmp[1] * _cou[2][4][1] + _faka * (_cou[0][4][0] - _fakaac * _cou[0][4][1]) + _fakac * _cou[2][1][1];
                    _cou[6][4][0] = _wmp[2] * _cou[2][4][1];
                    _cou[8][5][0] = _wmp[1] * _cou[2][5][1] + _faka * (_cou[0][5][0] - _fakaac * _cou[0][5][1]);
                    _cou[6][5][0] = _wmp[2] * _cou[2][5][1] + _fakac * _cou[2][1][1];
                    _cou[8][6][0] = _wmp[1] * _cou[2][6][1] + _faka * (_cou[0][6][0] - _fakaac * _cou[0][6][1]) + _fakac * _cou[2][3][1];
                    _cou[6][6][0] = _wmp[2] * _cou[2][6][1] + _fakac * _cou[2][2][1];
                    _cou[8][7][0] = _wmp[1] * _cou[2][7][1] + _faka * (_cou[0][7][0] - _fakaac * _cou[0][7][1]);
                    _cou[6][7][0] = _wmp[2] * _cou[2][7][1];
                    _cou[8][8][0] = _wmp[1] * _cou[2][8][1] + _faka * (_cou[0][8][0] - _fakaac * _cou[0][8][1]) + _fakac2 * _cou[2][2][1];
                    _cou[6][8][0] = _wmp[2] * _cou[2][8][1];
                    _cou[8][9][0] = _wmp[1] * _cou[2][9][1] + _faka * (_cou[0][9][0] - _fakaac * _cou[0][9][1]);
                    _cou[6][9][0] = _wmp[2] * _cou[2][9][1] + _fakac2 * _cou[2][3][1];
                    _cou[9][4][0] = _wmp[2] * _cou[3][4][1] + _faka * (_cou[0][4][0] - _fakaac * _cou[0][4][1]);
                    _cou[9][5][0] = _wmp[2] * _cou[3][5][1] + _faka * (_cou[0][5][0] - _fakaac * _cou[0][5][1]) + _fakac * _cou[3][1][1];
                    _cou[9][6][0] = _wmp[2] * _cou[3][6][1] + _faka * (_cou[0][6][0] - _fakaac * _cou[0][6][1]) + _fakac * _cou[3][2][1];
                    _cou[9][7][0] = _wmp[2] * _cou[3][7][1] + _faka * (_cou[0][7][0] - _fakaac * _cou[0][7][1]);
                    _cou[9][8][0] = _wmp[2] * _cou[3][8][1] + _faka * (_cou[0][8][0] - _fakaac * _cou[0][8][1]);
                    _cou[9][9][0] = _wmp[2] * _cou[3][9][1] + _faka * (_cou[0][9][0] - _fakaac * _cou[0][9][1]) + _fakac2 * _cou[3][3][1];
                }

                // p-f
                if (_lmax_col > 2) {
                    _cou[1][10][0] = _wmp[0] * _cou[0][10][1] + _fakac3 * _cou[0][7][1];
                    _cou[2][10][0] = _wmp[1] * _cou[0][10][1];
                    _cou[3][10][0] = _wmp[2] * _cou[0][10][1];
                    _cou[1][11][0] = _wmp[0] * _cou[0][11][1];
                    _cou[2][11][0] = _wmp[1] * _cou[0][11][1] + _fakac3 * _cou[0][8][1];
                    _cou[3][11][0] = _wmp[2] * _cou[0][11][1];
                    _cou[1][12][0] = _wmp[0] * _cou[0][12][1];
                    _cou[2][12][0] = _wmp[1] * _cou[0][12][1];
                    _cou[3][12][0] = _wmp[2] * _cou[0][12][1] + _fakac3 * _cou[0][9][1];
                    _cou[1][13][0] = _wmp[0] * _cou[0][13][1] + _fakac2 * _cou[0][4][1];
                    _cou[2][13][0] = _wmp[1] * _cou[0][13][1] + _fakac * _cou[0][7][1];
                    _cou[3][13][0] = _wmp[2] * _cou[0][13][1];
                    _cou[1][14][0] = _wmp[0] * _cou[0][14][1] + _fakac * _cou[0][8][1];
                    _cou[2][14][0] = _wmp[1] * _cou[0][14][1] + _fakac2 * _cou[0][4][1];
                    _cou[3][14][0] = _wmp[2] * _cou[0][14][1];
                    _cou[1][15][0] = _wmp[0] * _cou[0][15][1] + _fakac2 * _cou[0][5][1];
                    _cou[2][15][0] = _wmp[1] * _cou[0][15][1];
                    _cou[3][15][0] = _wmp[2] * _cou[0][15][1] + _fakac * _cou[0][7][1];
                    _cou[1][16][0] = _wmp[0] * _cou[0][16][1] + _fakac * _cou[0][9][1];
                    _cou[2][16][0] = _wmp[1] * _cou[0][16][1];
                    _cou[3][16][0] = _wmp[2] * _cou[0][16][1] + _fakac2 * _cou[0][5][1];
                    _cou[1][17][0] = _wmp[0] * _cou[0][17][1];
                    _cou[2][17][0] = _wmp[1] * _cou[0][17][1] + _fakac2 * _cou[0][6][1];
                    _cou[3][17][0] = _wmp[2] * _cou[0][17][1] + _fakac * _cou[0][8][1];
                    _cou[1][18][0] = _wmp[0] * _cou[0][18][1];
                    _cou[2][18][0] = _wmp[1] * _cou[0][18][1] + _fakac * _cou[0][9][1];
                    _cou[3][18][0] = _wmp[2] * _cou[0][18][1] + _fakac2 * _cou[0][6][1];
                    _cou[1][19][0] = _wmp[0] * _cou[0][19][1] + _fakac * _cou[0][6][1];
                    _cou[2][19][0] = _wmp[1] * _cou[0][19][1] + _fakac * _cou[0][5][1];
                    _cou[3][19][0] = _wmp[2] * _cou[0][19][1] + _fakac * _cou[0][4][1];
                }
            }


            if (_l_sum >= 5) {

                // p-s-4 (only g)

                // s-p-4 
                _cou[0][1][4] = _wmq[0] * _cou[0][0][5];
                _cou[0][2][4] = _wmq[1] * _cou[0][0][5];
                _cou[0][3][4] = _wmq[2] * _cou[0][0][5];

                // d-s-3 (only g)

                // p-p-3 (only g)

                // s-d-3
                _cou[0][7][3] = _wmq[0] * _cou[0][1][4] + _fakc * (_cou[0][0][3] - _fakaca * _cou[0][0][4]);
                _cou[0][4][3] = _wmq[1] * _cou[0][1][4];
                _cou[0][5][3] = _wmq[2] * _cou[0][1][4];
                _cou[0][8][3] = _wmq[1] * _cou[0][2][4] + _fakc * (_cou[0][0][3] - _fakaca * _cou[0][0][4]);
                _cou[0][6][3] = _wmq[2] * _cou[0][2][4];
                _cou[0][9][3] = _wmq[2] * _cou[0][3][4] + _fakc * (_cou[0][0][3] - _fakaca * _cou[0][0][4]);

                // p-d-2
                _cou[1][4][2] = _wmp[0] * _cou[0][4][3] + _fakac * _cou[0][2][3];
                _cou[2][4][2] = _wmp[1] * _cou[0][4][3] + _fakac * _cou[0][1][3];
                _cou[3][4][2] = _wmp[2] * _cou[0][4][3];
                _cou[1][5][2] = _wmp[0] * _cou[0][5][3] + _fakac * _cou[0][3][3];
                _cou[2][5][2] = _wmp[1] * _cou[0][5][3];
                _cou[3][5][2] = _wmp[2] * _cou[0][5][3] + _fakac * _cou[0][1][3];
                _cou[1][6][2] = _wmp[0] * _cou[0][6][3];
                _cou[2][6][2] = _wmp[1] * _cou[0][6][3] + _fakac * _cou[0][3][3];
                _cou[3][6][2] = _wmp[2] * _cou[0][6][3] + _fakac * _cou[0][2][3];
                _cou[1][7][2] = _wmp[0] * _cou[0][7][3] + _fakac2 * _cou[0][1][3];
                _cou[2][7][2] = _wmp[1] * _cou[0][7][3];
                _cou[3][7][2] = _wmp[2] * _cou[0][7][3];
                _cou[1][8][2] = _wmp[0] * _cou[0][8][3];
                _cou[2][8][2] = _wmp[1] * _cou[0][8][3] + _fakac2 * _cou[0][2][3];
                _cou[3][8][2] = _wmp[2] * _cou[0][8][3];
                _cou[1][9][2] = _wmp[0] * _cou[0][9][3];
                _cou[2][9][2] = _wmp[1] * _cou[0][9][3];
                _cou[3][9][2] = _wmp[2] * _cou[0][9][3] + _fakac2 * _cou[0][3][3];

                // s-f-2
                _cou[0][13][2] = _wmq[0] * _cou[0][4][3] + _fakc * (_cou[0][2][2] - _fakaca * _cou[0][2][3]);
                _cou[0][14][2] = _wmq[1] * _cou[0][4][3] + _fakc * (_cou[0][1][2] - _fakaca * _cou[0][1][3]);
                _cou[0][19][2] = _wmq[2] * _cou[0][4][3];
                _cou[0][15][2] = _wmq[0] * _cou[0][5][3] + _fakc * (_cou[0][3][2] - _fakaca * _cou[0][3][3]);
                _cou[0][16][2] = _wmq[2] * _cou[0][5][3] + _fakc * (_cou[0][1][2] - _fakaca * _cou[0][1][3]);
                _cou[0][17][2] = _wmq[1] * _cou[0][6][3] + _fakc * (_cou[0][3][2] - _fakaca * _cou[0][3][3]);
                _cou[0][18][2] = _wmq[2] * _cou[0][6][3] + _fakc * (_cou[0][2][2] - _fakaca * _cou[0][2][3]);
                _cou[0][10][2] = _wmq[0] * _cou[0][7][3] + _fakc2 * (_cou[0][1][2] - _fakaca * _cou[0][1][3]);
                _cou[0][11][2] = _wmq[1] * _cou[0][8][3] + _fakc2 * (_cou[0][2][2] - _fakaca * _cou[0][2][3]);
                _cou[0][12][2] = _wmq[2] * _cou[0][9][3] + _fakc2 * (_cou[0][3][2] - _fakaca * _cou[0][3][3]);

                // f-s-2 (only g)

                // d-p-2 (only g)

                // p-f-1
                _cou[1][10][1] = _wmp[0] * _cou[0][10][2] + _fakac3 * _cou[0][7][2];
                _cou[2][10][1] = _wmp[1] * _cou[0][10][2];
                _cou[3][10][1] = _wmp[2] * _cou[0][10][2];
                _cou[1][11][1] = _wmp[0] * _cou[0][11][2];
                _cou[2][11][1] = _wmp[1] * _cou[0][11][2] + _fakac3 * _cou[0][8][2];
                _cou[3][11][1] = _wmp[2] * _cou[0][11][2];
                _cou[1][12][1] = _wmp[0] * _cou[0][12][2];
                _cou[2][12][1] = _wmp[1] * _cou[0][12][2];
                _cou[3][12][1] = _wmp[2] * _cou[0][12][2] + _fakac3 * _cou[0][9][2];
                _cou[1][13][1] = _wmp[0] * _cou[0][13][2] + _fakac2 * _cou[0][4][2];
                _cou[2][13][1] = _wmp[1] * _cou[0][13][2] + _fakac * _cou[0][7][2];
                _cou[3][13][1] = _wmp[2] * _cou[0][13][2];
                _cou[1][14][1] = _wmp[0] * _cou[0][14][2] + _fakac * _cou[0][8][2];
                _cou[2][14][1] = _wmp[1] * _cou[0][14][2] + _fakac2 * _cou[0][4][2];
                _cou[3][14][1] = _wmp[2] * _cou[0][14][2];
                _cou[1][15][1] = _wmp[0] * _cou[0][15][2] + _fakac2 * _cou[0][5][2];
                _cou[2][15][1] = _wmp[1] * _cou[0][15][2];
                _cou[3][15][1] = _wmp[2] * _cou[0][15][2] + _fakac * _cou[0][7][2];
                _cou[1][16][1] = _wmp[0] * _cou[0][16][2] + _fakac * _cou[0][9][2];
                _cou[2][16][1] = _wmp[1] * _cou[0][16][2];
                _cou[3][16][1] = _wmp[2] * _cou[0][16][2] + _fakac2 * _cou[0][5][2];
                _cou[1][17][1] = _wmp[0] * _cou[0][17][2];
                _cou[2][17][1] = _wmp[1] * _cou[0][17][2] + _fakac2 * _cou[0][6][2];
                _cou[3][17][1] = _wmp[2] * _cou[0][17][2] + _fakac * _cou[0][8][2];
                _cou[1][18][1] = _wmp[0] * _cou[0][18][2];
                _cou[2][18][1] = _wmp[1] * _cou[0][18][2] + _fakac * _cou[0][9][2];
                _cou[3][18][1] = _wmp[2] * _cou[0][18][2] + _fakac2 * _cou[0][6][2];
                _cou[1][19][1] = _wmp[0] * _cou[0][19][2] + _fakac * _cou[0][6][2];
                _cou[2][19][1] = _wmp[1] * _cou[0][19][2] + _fakac * _cou[0][5][2];
                _cou[3][19][1] = _wmp[2] * _cou[0][19][2] + _fakac * _cou[0][4][2];

                // d-d-1
                _cou[7][4][1] = _wmp[0] * _cou[1][4][2] + _faka * (_cou[0][4][1] - _fakaac * _cou[0][4][2]) + _fakac * _cou[1][2][2];
                _cou[4][4][1] = _wmp[1] * _cou[1][4][2] + _fakac * _cou[1][1][2];
                _cou[5][4][1] = _wmp[2] * _cou[1][4][2];
                _cou[7][5][1] = _wmp[0] * _cou[1][5][2] + _faka * (_cou[0][5][1] - _fakaac * _cou[0][5][2]) + _fakac * _cou[1][3][2];
                _cou[4][5][1] = _wmp[1] * _cou[1][5][2];
                _cou[5][5][1] = _wmp[2] * _cou[1][5][2] + _fakac * _cou[1][1][2];
                _cou[7][6][1] = _wmp[0] * _cou[1][6][2] + _faka * (_cou[0][6][1] - _fakaac * _cou[0][6][2]);
                _cou[4][6][1] = _wmp[1] * _cou[1][6][2] + _fakac * _cou[1][3][2];
                _cou[5][6][1] = _wmp[2] * _cou[1][6][2] + _fakac * _cou[1][2][2];
                _cou[7][7][1] = _wmp[0] * _cou[1][7][2] + _faka * (_cou[0][7][1] - _fakaac * _cou[0][7][2]) + _fakac2 * _cou[1][1][2];
                _cou[4][7][1] = _wmp[1] * _cou[1][7][2];
                _cou[5][7][1] = _wmp[2] * _cou[1][7][2];
                _cou[7][8][1] = _wmp[0] * _cou[1][8][2] + _faka * (_cou[0][8][1] - _fakaac * _cou[0][8][2]);
                _cou[4][8][1] = _wmp[1] * _cou[1][8][2] + _fakac2 * _cou[1][2][2];
                _cou[5][8][1] = _wmp[2] * _cou[1][8][2];
                _cou[7][9][1] = _wmp[0] * _cou[1][9][2] + _faka * (_cou[0][9][1] - _fakaac * _cou[0][9][2]);
                _cou[4][9][1] = _wmp[1] * _cou[1][9][2];
                _cou[5][9][1] = _wmp[2] * _cou[1][9][2] + _fakac2 * _cou[1][3][2];
                _cou[8][4][1] = _wmp[1] * _cou[2][4][2] + _faka * (_cou[0][4][1] - _fakaac * _cou[0][4][2]) + _fakac * _cou[2][1][2];
                _cou[6][4][1] = _wmp[2] * _cou[2][4][2];
                _cou[8][5][1] = _wmp[1] * _cou[2][5][2] + _faka * (_cou[0][5][1] - _fakaac * _cou[0][5][2]);
                _cou[6][5][1] = _wmp[2] * _cou[2][5][2] + _fakac * _cou[2][1][2];
                _cou[8][6][1] = _wmp[1] * _cou[2][6][2] + _faka * (_cou[0][6][1] - _fakaac * _cou[0][6][2]) + _fakac * _cou[2][3][2];
                _cou[6][6][1] = _wmp[2] * _cou[2][6][2] + _fakac * _cou[2][2][2];
                _cou[8][7][1] = _wmp[1] * _cou[2][7][2] + _faka * (_cou[0][7][1] - _fakaac * _cou[0][7][2]);
                _cou[6][7][1] = _wmp[2] * _cou[2][7][2];
                _cou[8][8][1] = _wmp[1] * _cou[2][8][2] + _faka * (_cou[0][8][1] - _fakaac * _cou[0][8][2]) + _fakac2 * _cou[2][2][2];
                _cou[6][8][1] = _wmp[2] * _cou[2][8][2];
                _cou[8][9][1] = _wmp[1] * _cou[2][9][2] + _faka * (_cou[0][9][1] - _fakaac * _cou[0][9][2]);
                _cou[6][9][1] = _wmp[2] * _cou[2][9][2] + _fakac2 * _cou[2][3][2];
                _cou[9][4][1] = _wmp[2] * _cou[3][4][2] + _faka * (_cou[0][4][1] - _fakaac * _cou[0][4][2]);
                _cou[9][5][1] = _wmp[2] * _cou[3][5][2] + _faka * (_cou[0][5][1] - _fakaac * _cou[0][5][2]) + _fakac * _cou[3][1][2];
                _cou[9][6][1] = _wmp[2] * _cou[3][6][2] + _faka * (_cou[0][6][1] - _fakaac * _cou[0][6][2]) + _fakac * _cou[3][2][2];
                _cou[9][7][1] = _wmp[2] * _cou[3][7][2] + _faka * (_cou[0][7][1] - _fakaac * _cou[0][7][2]);
                _cou[9][8][1] = _wmp[2] * _cou[3][8][2] + _faka * (_cou[0][8][1] - _fakaac * _cou[0][8][2]);
                _cou[9][9][1] = _wmp[2] * _cou[3][9][2] + _faka * (_cou[0][9][1] - _fakaac * _cou[0][9][2]) + _fakac2 * _cou[3][3][2];

                // f-p-1 (only g)

                // d-f
                _cou[7][10][0] = _wmp[0] * _cou[1][10][1] + _faka * (_cou[0][10][0] - _fakaac * _cou[0][10][1]) + _fakac3 * _cou[1][7][1];
                _cou[4][10][0] = _wmp[1] * _cou[1][10][1];
                _cou[5][10][0] = _wmp[2] * _cou[1][10][1];
                _cou[7][11][0] = _wmp[0] * _cou[1][11][1] + _faka * (_cou[0][11][0] - _fakaac * _cou[0][11][1]);
                _cou[4][11][0] = _wmp[1] * _cou[1][11][1] + _fakac3 * _cou[1][8][1];
                _cou[5][11][0] = _wmp[2] * _cou[1][11][1];
                _cou[7][12][0] = _wmp[0] * _cou[1][12][1] + _faka * (_cou[0][12][0] - _fakaac * _cou[0][12][1]);
                _cou[4][12][0] = _wmp[1] * _cou[1][12][1];
                _cou[5][12][0] = _wmp[2] * _cou[1][12][1] + _fakac3 * _cou[1][9][1];
                _cou[7][13][0] = _wmp[0] * _cou[1][13][1] + _faka * (_cou[0][13][0] - _fakaac * _cou[0][13][1]) + _fakac2 * _cou[1][4][1];
                _cou[4][13][0] = _wmp[1] * _cou[1][13][1] + _fakac * _cou[1][7][1];
                _cou[5][13][0] = _wmp[2] * _cou[1][13][1];
                _cou[7][14][0] = _wmp[0] * _cou[1][14][1] + _faka * (_cou[0][14][0] - _fakaac * _cou[0][14][1]) + _fakac * _cou[1][8][1];
                _cou[4][14][0] = _wmp[1] * _cou[1][14][1] + _fakac2 * _cou[1][4][1];
                _cou[5][14][0] = _wmp[2] * _cou[1][14][1];
                _cou[7][15][0] = _wmp[0] * _cou[1][15][1] + _faka * (_cou[0][15][0] - _fakaac * _cou[0][15][1]) + _fakac2 * _cou[1][5][1];
                _cou[4][15][0] = _wmp[1] * _cou[1][15][1];
                _cou[5][15][0] = _wmp[2] * _cou[1][15][1] + _fakac * _cou[1][7][1];
                _cou[7][16][0] = _wmp[0] * _cou[1][16][1] + _faka * (_cou[0][16][0] - _fakaac * _cou[0][16][1]) + _fakac * _cou[1][9][1];
                _cou[4][16][0] = _wmp[1] * _cou[1][16][1];
                _cou[5][16][0] = _wmp[2] * _cou[1][16][1] + _fakac2 * _cou[1][5][1];
                _cou[7][17][0] = _wmp[0] * _cou[1][17][1] + _faka * (_cou[0][17][0] - _fakaac * _cou[0][17][1]);
                _cou[4][17][0] = _wmp[1] * _cou[1][17][1] + _fakac2 * _cou[1][6][1];
                _cou[5][17][0] = _wmp[2] * _cou[1][17][1] + _fakac * _cou[1][8][1];
                _cou[7][18][0] = _wmp[0] * _cou[1][18][1] + _faka * (_cou[0][18][0] - _fakaac * _cou[0][18][1]);
                _cou[4][18][0] = _wmp[1] * _cou[1][18][1] + _fakac * _cou[1][9][1];
                _cou[5][18][0] = _wmp[2] * _cou[1][18][1] + _fakac2 * _cou[1][6][1];
                _cou[7][19][0] = _wmp[0] * _cou[1][19][1] + _faka * (_cou[0][19][0] - _fakaac * _cou[0][19][1]) + _fakac * _cou[1][6][1];
                _cou[4][19][0] = _wmp[1] * _cou[1][19][1] + _fakac * _cou[1][5][1];
                _cou[5][19][0] = _wmp[2] * _cou[1][19][1] + _fakac * _cou[1][4][1];
                _cou[8][10][0] = _wmp[1] * _cou[2][10][1] + _faka * (_cou[0][10][0] - _fakaac * _cou[0][10][1]);
                _cou[6][10][0] = _wmp[2] * _cou[2][10][1];
                _cou[8][11][0] = _wmp[1] * _cou[2][11][1] + _faka * (_cou[0][11][0] - _fakaac * _cou[0][11][1]) + _fakac3 * _cou[2][8][1];
                _cou[6][11][0] = _wmp[2] * _cou[2][11][1];
                _cou[8][12][0] = _wmp[1] * _cou[2][12][1] + _faka * (_cou[0][12][0] - _fakaac * _cou[0][12][1]);
                _cou[6][12][0] = _wmp[2] * _cou[2][12][1] + _fakac3 * _cou[2][9][1];
                _cou[8][13][0] = _wmp[1] * _cou[2][13][1] + _faka * (_cou[0][13][0] - _fakaac * _cou[0][13][1]) + _fakac * _cou[2][7][1];
                _cou[6][13][0] = _wmp[2] * _cou[2][13][1];
                _cou[8][14][0] = _wmp[1] * _cou[2][14][1] + _faka * (_cou[0][14][0] - _fakaac * _cou[0][14][1]) + _fakac2 * _cou[2][4][1];
                _cou[6][14][0] = _wmp[2] * _cou[2][14][1];
                _cou[8][15][0] = _wmp[1] * _cou[2][15][1] + _faka * (_cou[0][15][0] - _fakaac * _cou[0][15][1]);
                _cou[6][15][0] = _wmp[2] * _cou[2][15][1] + _fakac * _cou[2][7][1];
                _cou[8][16][0] = _wmp[1] * _cou[2][16][1] + _faka * (_cou[0][16][0] - _fakaac * _cou[0][16][1]);
                _cou[6][16][0] = _wmp[2] * _cou[2][16][1] + _fakac2 * _cou[2][5][1];
                _cou[8][17][0] = _wmp[1] * _cou[2][17][1] + _faka * (_cou[0][17][0] - _fakaac * _cou[0][17][1]) + _fakac2 * _cou[2][6][1];
                _cou[6][17][0] = _wmp[2] * _cou[2][17][1] + _fakac * _cou[2][8][1];
                _cou[8][18][0] = _wmp[1] * _cou[2][18][1] + _faka * (_cou[0][18][0] - _fakaac * _cou[0][18][1]) + _fakac * _cou[2][9][1];
                _cou[6][18][0] = _wmp[2] * _cou[2][18][1] + _fakac2 * _cou[2][6][1];
                _cou[8][19][0] = _wmp[1] * _cou[2][19][1] + _faka * (_cou[0][19][0] - _fakaac * _cou[0][19][1]) + _fakac * _cou[2][5][1];
                _cou[6][19][0] = _wmp[2] * _cou[2][19][1] + _fakac * _cou[2][4][1];
                _cou[9][10][0] = _wmp[2] * _cou[3][10][1] + _faka * (_cou[0][10][0] - _fakaac * _cou[0][10][1]);
                _cou[9][11][0] = _wmp[2] * _cou[3][11][1] + _faka * (_cou[0][11][0] - _fakaac * _cou[0][11][1]);
                _cou[9][12][0] = _wmp[2] * _cou[3][12][1] + _faka * (_cou[0][12][0] - _fakaac * _cou[0][12][1]) + _fakac3 * _cou[3][9][1];
                _cou[9][13][0] = _wmp[2] * _cou[3][13][1] + _faka * (_cou[0][13][0] - _fakaac * _cou[0][13][1]);
                _cou[9][14][0] = _wmp[2] * _cou[3][14][1] + _faka * (_cou[0][14][0] - _fakaac * _cou[0][14][1]);
                _cou[9][15][0] = _wmp[2] * _cou[3][15][1] + _faka * (_cou[0][15][0] - _fakaac * _cou[0][15][1]) + _fakac * _cou[3][7][1];
                _cou[9][16][0] = _wmp[2] * _cou[3][16][1] + _faka * (_cou[0][16][0] - _fakaac * _cou[0][16][1]) + _fakac2 * _cou[3][5][1];
                _cou[9][17][0] = _wmp[2] * _cou[3][17][1] + _faka * (_cou[0][17][0] - _fakaac * _cou[0][17][1]) + _fakac * _cou[3][8][1];
                _cou[9][18][0] = _wmp[2] * _cou[3][18][1] + _faka * (_cou[0][18][0] - _fakaac * _cou[0][18][1]) + _fakac2 * _cou[3][6][1];
                _cou[9][19][0] = _wmp[2] * _cou[3][19][1] + _faka * (_cou[0][19][0] - _fakaac * _cou[0][19][1]) + _fakac * _cou[3][4][1];

                // f-d 
                _cou[13][4][0] = _wmp[0] * _cou[4][4][1] + _faka * (_cou[2][4][0] - _fakaac * _cou[2][4][1]) + _fakac * _cou[4][2][1];
                _cou[14][4][0] = _wmp[1] * _cou[4][4][1] + _faka * (_cou[1][4][0] - _fakaac * _cou[1][4][1]) + _fakac * _cou[4][1][1];
                _cou[19][4][0] = _wmp[2] * _cou[4][4][1];
                _cou[13][5][0] = _wmp[0] * _cou[4][5][1] + _faka * (_cou[2][5][0] - _fakaac * _cou[2][5][1]) + _fakac * _cou[4][3][1];
                _cou[14][5][0] = _wmp[1] * _cou[4][5][1] + _faka * (_cou[1][5][0] - _fakaac * _cou[1][5][1]);
                _cou[19][5][0] = _wmp[2] * _cou[4][5][1] + _fakac * _cou[4][1][1];
                _cou[13][6][0] = _wmp[0] * _cou[4][6][1] + _faka * (_cou[2][6][0] - _fakaac * _cou[2][6][1]);
                _cou[14][6][0] = _wmp[1] * _cou[4][6][1] + _faka * (_cou[1][6][0] - _fakaac * _cou[1][6][1]) + _fakac * _cou[4][3][1];
                _cou[19][6][0] = _wmp[2] * _cou[4][6][1] + _fakac * _cou[4][2][1];
                _cou[13][7][0] = _wmp[0] * _cou[4][7][1] + _faka * (_cou[2][7][0] - _fakaac * _cou[2][7][1]) + _fakac2 * _cou[4][1][1];
                _cou[14][7][0] = _wmp[1] * _cou[4][7][1] + _faka * (_cou[1][7][0] - _fakaac * _cou[1][7][1]);
                _cou[19][7][0] = _wmp[2] * _cou[4][7][1];
                _cou[13][8][0] = _wmp[0] * _cou[4][8][1] + _faka * (_cou[2][8][0] - _fakaac * _cou[2][8][1]);
                _cou[14][8][0] = _wmp[1] * _cou[4][8][1] + _faka * (_cou[1][8][0] - _fakaac * _cou[1][8][1]) + _fakac2 * _cou[4][2][1];
                _cou[19][8][0] = _wmp[2] * _cou[4][8][1];
                _cou[13][9][0] = _wmp[0] * _cou[4][9][1] + _faka * (_cou[2][9][0] - _fakaac * _cou[2][9][1]);
                _cou[14][9][0] = _wmp[1] * _cou[4][9][1] + _faka * (_cou[1][9][0] - _fakaac * _cou[1][9][1]);
                _cou[19][9][0] = _wmp[2] * _cou[4][9][1] + _fakac2 * _cou[4][3][1];
                _cou[15][4][0] = _wmp[0] * _cou[5][4][1] + _faka * (_cou[3][4][0] - _fakaac * _cou[3][4][1]) + _fakac * _cou[5][2][1];
                _cou[16][4][0] = _wmp[2] * _cou[5][4][1] + _faka * (_cou[1][4][0] - _fakaac * _cou[1][4][1]);
                _cou[15][5][0] = _wmp[0] * _cou[5][5][1] + _faka * (_cou[3][5][0] - _fakaac * _cou[3][5][1]) + _fakac * _cou[5][3][1];
                _cou[16][5][0] = _wmp[2] * _cou[5][5][1] + _faka * (_cou[1][5][0] - _fakaac * _cou[1][5][1]) + _fakac * _cou[5][1][1];
                _cou[15][6][0] = _wmp[0] * _cou[5][6][1] + _faka * (_cou[3][6][0] - _fakaac * _cou[3][6][1]);
                _cou[16][6][0] = _wmp[2] * _cou[5][6][1] + _faka * (_cou[1][6][0] - _fakaac * _cou[1][6][1]) + _fakac * _cou[5][2][1];
                _cou[15][7][0] = _wmp[0] * _cou[5][7][1] + _faka * (_cou[3][7][0] - _fakaac * _cou[3][7][1]) + _fakac2 * _cou[5][1][1];
                _cou[16][7][0] = _wmp[2] * _cou[5][7][1] + _faka * (_cou[1][7][0] - _fakaac * _cou[1][7][1]);
                _cou[15][8][0] = _wmp[0] * _cou[5][8][1] + _faka * (_cou[3][8][0] - _fakaac * _cou[3][8][1]);
                _cou[16][8][0] = _wmp[2] * _cou[5][8][1] + _faka * (_cou[1][8][0] - _fakaac * _cou[1][8][1]);
                _cou[15][9][0] = _wmp[0] * _cou[5][9][1] + _faka * (_cou[3][9][0] - _fakaac * _cou[3][9][1]);
                _cou[16][9][0] = _wmp[2] * _cou[5][9][1] + _faka * (_cou[1][9][0] - _fakaac * _cou[1][9][1]) + _fakac2 * _cou[5][3][1];
                _cou[17][4][0] = _wmp[1] * _cou[6][4][1] + _faka * (_cou[3][4][0] - _fakaac * _cou[3][4][1]) + _fakac * _cou[6][1][1];
                _cou[18][4][0] = _wmp[2] * _cou[6][4][1] + _faka * (_cou[2][4][0] - _fakaac * _cou[2][4][1]);
                _cou[17][5][0] = _wmp[1] * _cou[6][5][1] + _faka * (_cou[3][5][0] - _fakaac * _cou[3][5][1]);
                _cou[18][5][0] = _wmp[2] * _cou[6][5][1] + _faka * (_cou[2][5][0] - _fakaac * _cou[2][5][1]) + _fakac * _cou[6][1][1];
                _cou[17][6][0] = _wmp[1] * _cou[6][6][1] + _faka * (_cou[3][6][0] - _fakaac * _cou[3][6][1]) + _fakac * _cou[6][3][1];
                _cou[18][6][0] = _wmp[2] * _cou[6][6][1] + _faka * (_cou[2][6][0] - _fakaac * _cou[2][6][1]) + _fakac * _cou[6][2][1];
                _cou[17][7][0] = _wmp[1] * _cou[6][7][1] + _faka * (_cou[3][7][0] - _fakaac * _cou[3][7][1]);
                _cou[18][7][0] = _wmp[2] * _cou[6][7][1] + _faka * (_cou[2][7][0] - _fakaac * _cou[2][7][1]);
                _cou[17][8][0] = _wmp[1] * _cou[6][8][1] + _faka * (_cou[3][8][0] - _fakaac * _cou[3][8][1]) + _fakac2 * _cou[6][2][1];
                _cou[18][8][0] = _wmp[2] * _cou[6][8][1] + _faka * (_cou[2][8][0] - _fakaac * _cou[2][8][1]);
                _cou[17][9][0] = _wmp[1] * _cou[6][9][1] + _faka * (_cou[3][9][0] - _fakaac * _cou[3][9][1]);
                _cou[18][9][0] = _wmp[2] * _cou[6][9][1] + _faka * (_cou[2][9][0] - _fakaac * _cou[2][9][1]) + _fakac2 * _cou[6][3][1];
                _cou[10][4][0] = _wmp[0] * _cou[7][4][1] + _faka2 * (_cou[1][4][0] - _fakaac * _cou[1][4][1]) + _fakac * _cou[7][2][1];
                _cou[10][5][0] = _wmp[0] * _cou[7][5][1] + _faka2 * (_cou[1][5][0] - _fakaac * _cou[1][5][1]) + _fakac * _cou[7][3][1];
                _cou[10][6][0] = _wmp[0] * _cou[7][6][1] + _faka2 * (_cou[1][6][0] - _fakaac * _cou[1][6][1]);
                _cou[10][7][0] = _wmp[0] * _cou[7][7][1] + _faka2 * (_cou[1][7][0] - _fakaac * _cou[1][7][1]) + _fakac2 * _cou[7][1][1];
                _cou[10][8][0] = _wmp[0] * _cou[7][8][1] + _faka2 * (_cou[1][8][0] - _fakaac * _cou[1][8][1]);
                _cou[10][9][0] = _wmp[0] * _cou[7][9][1] + _faka2 * (_cou[1][9][0] - _fakaac * _cou[1][9][1]);
                _cou[11][4][0] = _wmp[1] * _cou[8][4][1] + _faka2 * (_cou[2][4][0] - _fakaac * _cou[2][4][1]) + _fakac * _cou[8][1][1];
                _cou[11][5][0] = _wmp[1] * _cou[8][5][1] + _faka2 * (_cou[2][5][0] - _fakaac * _cou[2][5][1]);
                _cou[11][6][0] = _wmp[1] * _cou[8][6][1] + _faka2 * (_cou[2][6][0] - _fakaac * _cou[2][6][1]) + _fakac * _cou[8][3][1];
                _cou[11][7][0] = _wmp[1] * _cou[8][7][1] + _faka2 * (_cou[2][7][0] - _fakaac * _cou[2][7][1]);
                _cou[11][8][0] = _wmp[1] * _cou[8][8][1] + _faka2 * (_cou[2][8][0] - _fakaac * _cou[2][8][1]) + _fakac2 * _cou[8][2][1];
                _cou[11][9][0] = _wmp[1] * _cou[8][9][1] + _faka2 * (_cou[2][9][0] - _fakaac * _cou[2][9][1]);
                _cou[12][4][0] = _wmp[2] * _cou[9][4][1] + _faka2 * (_cou[3][4][0] - _fakaac * _cou[3][4][1]);
                _cou[12][5][0] = _wmp[2] * _cou[9][5][1] + _faka2 * (_cou[3][5][0] - _fakaac * _cou[3][5][1]) + _fakac * _cou[9][1][1];
                _cou[12][6][0] = _wmp[2] * _cou[9][6][1] + _faka2 * (_cou[3][6][0] - _fakaac * _cou[3][6][1]) + _fakac * _cou[9][2][1];
                _cou[12][7][0] = _wmp[2] * _cou[9][7][1] + _faka2 * (_cou[3][7][0] - _fakaac * _cou[3][7][1]);
                _cou[12][8][0] = _wmp[2] * _cou[9][8][1] + _faka2 * (_cou[3][8][0] - _fakaac * _cou[3][8][1]);
                _cou[12][9][0] = _wmp[2] * _cou[9][9][1] + _faka2 * (_cou[3][9][0] - _fakaac * _cou[3][9][1]) + _fakac2 * _cou[9][3][1];


            }

            if (_l_sum >= 6) {

                // p-s-5 (only g)

                // s-p-5
                _cou[0][1][5] = _wmq[0] * _cou[0][0][6];
                _cou[0][2][5] = _wmq[1] * _cou[0][0][6];
                _cou[0][3][5] = _wmq[2] * _cou[0][0][6];

                // p-p-4

                // s-d-4
                _cou[0][7][4] = _wmq[0] * _cou[0][1][5] + _fakc * (_cou[0][0][4] - _fakaca * _cou[0][0][5]);
                _cou[0][4][4] = _wmq[1] * _cou[0][1][5];
                _cou[0][5][4] = _wmq[2] * _cou[0][1][5];
                _cou[0][8][4] = _wmq[1] * _cou[0][2][5] + _fakc * (_cou[0][0][4] - _fakaca * _cou[0][0][5]);
                _cou[0][6][4] = _wmq[2] * _cou[0][2][5];
                _cou[0][9][4] = _wmq[2] * _cou[0][3][5] + _fakc * (_cou[0][0][4] - _fakaca * _cou[0][0][5]);

                // d-s-4

                // p-d-3

                // s-f-3
                _cou[0][13][3] = _wmq[0] * _cou[0][4][4] + _fakc * (_cou[0][2][3] - _fakaca * _cou[0][2][4]);
                _cou[0][14][3] = _wmq[1] * _cou[0][4][4] + _fakc * (_cou[0][1][3] - _fakaca * _cou[0][1][4]);
                _cou[0][19][3] = _wmq[2] * _cou[0][4][4];
                _cou[0][15][3] = _wmq[0] * _cou[0][5][4] + _fakc * (_cou[0][3][3] - _fakaca * _cou[0][3][4]);
                _cou[0][16][3] = _wmq[2] * _cou[0][5][4] + _fakc * (_cou[0][1][3] - _fakaca * _cou[0][1][4]);
                _cou[0][17][3] = _wmq[1] * _cou[0][6][4] + _fakc * (_cou[0][3][3] - _fakaca * _cou[0][3][4]);
                _cou[0][18][3] = _wmq[2] * _cou[0][6][4] + _fakc * (_cou[0][2][3] - _fakaca * _cou[0][2][4]);
                _cou[0][10][3] = _wmq[0] * _cou[0][7][4] + _fakc2 * (_cou[0][1][3] - _fakaca * _cou[0][1][4]);
                _cou[0][11][3] = _wmq[1] * _cou[0][8][4] + _fakc2 * (_cou[0][2][3] - _fakaca * _cou[0][2][4]);
                _cou[0][12][3] = _wmq[2] * _cou[0][9][4] + _fakc2 * (_cou[0][3][3] - _fakaca * _cou[0][3][4]);

                // d-p-3

                // f-s-3

                // p-f-2
                _cou[1][10][2] = _wmp[0] * _cou[0][10][3] + _fakac3 * _cou[0][7][3];
                _cou[2][10][2] = _wmp[1] * _cou[0][10][3];
                _cou[3][10][2] = _wmp[2] * _cou[0][10][3];
                _cou[1][11][2] = _wmp[0] * _cou[0][11][3];
                _cou[2][11][2] = _wmp[1] * _cou[0][11][3] + _fakac3 * _cou[0][8][3];
                _cou[3][11][2] = _wmp[2] * _cou[0][11][3];
                _cou[1][12][2] = _wmp[0] * _cou[0][12][3];
                _cou[2][12][2] = _wmp[1] * _cou[0][12][3];
                _cou[3][12][2] = _wmp[2] * _cou[0][12][3] + _fakac3 * _cou[0][9][3];
                _cou[1][13][2] = _wmp[0] * _cou[0][13][3] + _fakac2 * _cou[0][4][3];
                _cou[2][13][2] = _wmp[1] * _cou[0][13][3] + _fakac * _cou[0][7][3];
                _cou[3][13][2] = _wmp[2] * _cou[0][13][3];
                _cou[1][14][2] = _wmp[0] * _cou[0][14][3] + _fakac * _cou[0][8][3];
                _cou[2][14][2] = _wmp[1] * _cou[0][14][3] + _fakac2 * _cou[0][4][3];
                _cou[3][14][2] = _wmp[2] * _cou[0][14][3];
                _cou[1][15][2] = _wmp[0] * _cou[0][15][3] + _fakac2 * _cou[0][5][3];
                _cou[2][15][2] = _wmp[1] * _cou[0][15][3];
                _cou[3][15][2] = _wmp[2] * _cou[0][15][3] + _fakac * _cou[0][7][3];
                _cou[1][16][2] = _wmp[0] * _cou[0][16][3] + _fakac * _cou[0][9][3];
                _cou[2][16][2] = _wmp[1] * _cou[0][16][3];
                _cou[3][16][2] = _wmp[2] * _cou[0][16][3] + _fakac2 * _cou[0][5][3];
                _cou[1][17][2] = _wmp[0] * _cou[0][17][3];
                _cou[2][17][2] = _wmp[1] * _cou[0][17][3] + _fakac2 * _cou[0][6][3];
                _cou[3][17][2] = _wmp[2] * _cou[0][17][3] + _fakac * _cou[0][8][3];
                _cou[1][18][2] = _wmp[0] * _cou[0][18][3];
                _cou[2][18][2] = _wmp[1] * _cou[0][18][3] + _fakac * _cou[0][9][3];
                _cou[3][18][2] = _wmp[2] * _cou[0][18][3] + _fakac2 * _cou[0][6][3];
                _cou[1][19][2] = _wmp[0] * _cou[0][19][3] + _fakac * _cou[0][6][3];
                _cou[2][19][2] = _wmp[1] * _cou[0][19][3] + _fakac * _cou[0][5][3];
                _cou[3][19][2] = _wmp[2] * _cou[0][19][3] + _fakac * _cou[0][4][3];
                // d-d-2

                // f-p-2

                // d-f-1
                _cou[7][10][1] = _wmp[0] * _cou[1][10][2] + _faka * (_cou[0][10][1] - _fakaac * _cou[0][10][2]) + _fakac3 * _cou[1][7][2];
                _cou[4][10][1] = _wmp[1] * _cou[1][10][2];
                _cou[5][10][1] = _wmp[2] * _cou[1][10][2];
                _cou[7][11][1] = _wmp[0] * _cou[1][11][2] + _faka * (_cou[0][11][1] - _fakaac * _cou[0][11][2]);
                _cou[4][11][1] = _wmp[1] * _cou[1][11][2] + _fakac3 * _cou[1][8][2];
                _cou[5][11][1] = _wmp[2] * _cou[1][11][2];
                _cou[7][12][1] = _wmp[0] * _cou[1][12][2] + _faka * (_cou[0][12][1] - _fakaac * _cou[0][12][2]);
                _cou[4][12][1] = _wmp[1] * _cou[1][12][2];
                _cou[5][12][1] = _wmp[2] * _cou[1][12][2] + _fakac3 * _cou[1][9][2];
                _cou[7][13][1] = _wmp[0] * _cou[1][13][2] + _faka * (_cou[0][13][1] - _fakaac * _cou[0][13][2]) + _fakac2 * _cou[1][4][2];
                _cou[4][13][1] = _wmp[1] * _cou[1][13][2] + _fakac * _cou[1][7][2];
                _cou[5][13][1] = _wmp[2] * _cou[1][13][2];
                _cou[7][14][1] = _wmp[0] * _cou[1][14][2] + _faka * (_cou[0][14][1] - _fakaac * _cou[0][14][2]) + _fakac * _cou[1][8][2];
                _cou[4][14][1] = _wmp[1] * _cou[1][14][2] + _fakac2 * _cou[1][4][2];
                _cou[5][14][1] = _wmp[2] * _cou[1][14][2];
                _cou[7][15][1] = _wmp[0] * _cou[1][15][2] + _faka * (_cou[0][15][1] - _fakaac * _cou[0][15][2]) + _fakac2 * _cou[1][5][2];
                _cou[4][15][1] = _wmp[1] * _cou[1][15][2];
                _cou[5][15][1] = _wmp[2] * _cou[1][15][2] + _fakac * _cou[1][7][2];
                _cou[7][16][1] = _wmp[0] * _cou[1][16][2] + _faka * (_cou[0][16][1] - _fakaac * _cou[0][16][2]) + _fakac * _cou[1][9][2];
                _cou[4][16][1] = _wmp[1] * _cou[1][16][2];
                _cou[5][16][1] = _wmp[2] * _cou[1][16][2] + _fakac2 * _cou[1][5][2];
                _cou[7][17][1] = _wmp[0] * _cou[1][17][2] + _faka * (_cou[0][17][1] - _fakaac * _cou[0][17][2]);
                _cou[4][17][1] = _wmp[1] * _cou[1][17][2] + _fakac2 * _cou[1][6][2];
                _cou[5][17][1] = _wmp[2] * _cou[1][17][2] + _fakac * _cou[1][8][2];
                _cou[7][18][1] = _wmp[0] * _cou[1][18][2] + _faka * (_cou[0][18][1] - _fakaac * _cou[0][18][2]);
                _cou[4][18][1] = _wmp[1] * _cou[1][18][2] + _fakac * _cou[1][9][2];
                _cou[5][18][1] = _wmp[2] * _cou[1][18][2] + _fakac2 * _cou[1][6][2];
                _cou[7][19][1] = _wmp[0] * _cou[1][19][2] + _faka * (_cou[0][19][1] - _fakaac * _cou[0][19][2]) + _fakac * _cou[1][6][2];
                _cou[4][19][1] = _wmp[1] * _cou[1][19][2] + _fakac * _cou[1][5][2];
                _cou[5][19][1] = _wmp[2] * _cou[1][19][2] + _fakac * _cou[1][4][2];
                _cou[8][10][1] = _wmp[1] * _cou[2][10][2] + _faka * (_cou[0][10][1] - _fakaac * _cou[0][10][2]);
                _cou[6][10][1] = _wmp[2] * _cou[2][10][2];
                _cou[8][11][1] = _wmp[1] * _cou[2][11][2] + _faka * (_cou[0][11][1] - _fakaac * _cou[0][11][2]) + _fakac3 * _cou[2][8][2];
                _cou[6][11][1] = _wmp[2] * _cou[2][11][2];
                _cou[8][12][1] = _wmp[1] * _cou[2][12][2] + _faka * (_cou[0][12][1] - _fakaac * _cou[0][12][2]);
                _cou[6][12][1] = _wmp[2] * _cou[2][12][2] + _fakac3 * _cou[2][9][2];
                _cou[8][13][1] = _wmp[1] * _cou[2][13][2] + _faka * (_cou[0][13][1] - _fakaac * _cou[0][13][2]) + _fakac * _cou[2][7][2];
                _cou[6][13][1] = _wmp[2] * _cou[2][13][2];
                _cou[8][14][1] = _wmp[1] * _cou[2][14][2] + _faka * (_cou[0][14][1] - _fakaac * _cou[0][14][2]) + _fakac2 * _cou[2][4][2];
                _cou[6][14][1] = _wmp[2] * _cou[2][14][2];
                _cou[8][15][1] = _wmp[1] * _cou[2][15][2] + _faka * (_cou[0][15][1] - _fakaac * _cou[0][15][2]);
                _cou[6][15][1] = _wmp[2] * _cou[2][15][2] + _fakac * _cou[2][7][2];
                _cou[8][16][1] = _wmp[1] * _cou[2][16][2] + _faka * (_cou[0][16][1] - _fakaac * _cou[0][16][2]);
                _cou[6][16][1] = _wmp[2] * _cou[2][16][2] + _fakac2 * _cou[2][5][2];
                _cou[8][17][1] = _wmp[1] * _cou[2][17][2] + _faka * (_cou[0][17][1] - _fakaac * _cou[0][17][2]) + _fakac2 * _cou[2][6][2];
                _cou[6][17][1] = _wmp[2] * _cou[2][17][2] + _fakac * _cou[2][8][2];
                _cou[8][18][1] = _wmp[1] * _cou[2][18][2] + _faka * (_cou[0][18][1] - _fakaac * _cou[0][18][2]) + _fakac * _cou[2][9][2];
                _cou[6][18][1] = _wmp[2] * _cou[2][18][2] + _fakac2 * _cou[2][6][2];
                _cou[8][19][1] = _wmp[1] * _cou[2][19][2] + _faka * (_cou[0][19][1] - _fakaac * _cou[0][19][2]) + _fakac * _cou[2][5][2];
                _cou[6][19][1] = _wmp[2] * _cou[2][19][2] + _fakac * _cou[2][4][2];
                _cou[9][10][1] = _wmp[2] * _cou[3][10][2] + _faka * (_cou[0][10][1] - _fakaac * _cou[0][10][2]);
                _cou[9][11][1] = _wmp[2] * _cou[3][11][2] + _faka * (_cou[0][11][1] - _fakaac * _cou[0][11][2]);
                _cou[9][12][1] = _wmp[2] * _cou[3][12][2] + _faka * (_cou[0][12][1] - _fakaac * _cou[0][12][2]) + _fakac3 * _cou[3][9][2];
                _cou[9][13][1] = _wmp[2] * _cou[3][13][2] + _faka * (_cou[0][13][1] - _fakaac * _cou[0][13][2]);
                _cou[9][14][1] = _wmp[2] * _cou[3][14][2] + _faka * (_cou[0][14][1] - _fakaac * _cou[0][14][2]);
                _cou[9][15][1] = _wmp[2] * _cou[3][15][2] + _faka * (_cou[0][15][1] - _fakaac * _cou[0][15][2]) + _fakac * _cou[3][7][2];
                _cou[9][16][1] = _wmp[2] * _cou[3][16][2] + _faka * (_cou[0][16][1] - _fakaac * _cou[0][16][2]) + _fakac2 * _cou[3][5][2];
                _cou[9][17][1] = _wmp[2] * _cou[3][17][2] + _faka * (_cou[0][17][1] - _fakaac * _cou[0][17][2]) + _fakac * _cou[3][8][2];
                _cou[9][18][1] = _wmp[2] * _cou[3][18][2] + _faka * (_cou[0][18][1] - _fakaac * _cou[0][18][2]) + _fakac2 * _cou[3][6][2];
                _cou[9][19][1] = _wmp[2] * _cou[3][19][2] + _faka * (_cou[0][19][1] - _fakaac * _cou[0][19][2]) + _fakac * _cou[3][4][2];

                // f-d-1

                // f-f
                _cou[13][10][0] = _wmp[0] * _cou[4][10][1] + _faka * (_cou[2][10][0] - _fakaac * _cou[2][10][1]) + _fakac3 * _cou[4][7][1];
                _cou[14][10][0] = _wmp[1] * _cou[4][10][1] + _faka * (_cou[1][10][0] - _fakaac * _cou[1][10][1]);
                _cou[19][10][0] = _wmp[2] * _cou[4][10][1];
                _cou[13][11][0] = _wmp[0] * _cou[4][11][1] + _faka * (_cou[2][11][0] - _fakaac * _cou[2][11][1]);
                _cou[14][11][0] = _wmp[1] * _cou[4][11][1] + _faka * (_cou[1][11][0] - _fakaac * _cou[1][11][1]) + _fakac3 * _cou[4][8][1];
                _cou[19][11][0] = _wmp[2] * _cou[4][11][1];
                _cou[13][12][0] = _wmp[0] * _cou[4][12][1] + _faka * (_cou[2][12][0] - _fakaac * _cou[2][12][1]);
                _cou[14][12][0] = _wmp[1] * _cou[4][12][1] + _faka * (_cou[1][12][0] - _fakaac * _cou[1][12][1]);
                _cou[19][12][0] = _wmp[2] * _cou[4][12][1] + _fakac3 * _cou[4][9][1];
                _cou[13][13][0] = _wmp[0] * _cou[4][13][1] + _faka * (_cou[2][13][0] - _fakaac * _cou[2][13][1]) + _fakac2 * _cou[4][4][1];
                _cou[14][13][0] = _wmp[1] * _cou[4][13][1] + _faka * (_cou[1][13][0] - _fakaac * _cou[1][13][1]) + _fakac * _cou[4][7][1];
                _cou[19][13][0] = _wmp[2] * _cou[4][13][1];
                _cou[13][14][0] = _wmp[0] * _cou[4][14][1] + _faka * (_cou[2][14][0] - _fakaac * _cou[2][14][1]) + _fakac * _cou[4][8][1];
                _cou[14][14][0] = _wmp[1] * _cou[4][14][1] + _faka * (_cou[1][14][0] - _fakaac * _cou[1][14][1]) + _fakac2 * _cou[4][4][1];
                _cou[19][14][0] = _wmp[2] * _cou[4][14][1];
                _cou[13][15][0] = _wmp[0] * _cou[4][15][1] + _faka * (_cou[2][15][0] - _fakaac * _cou[2][15][1]) + _fakac2 * _cou[4][5][1];
                _cou[14][15][0] = _wmp[1] * _cou[4][15][1] + _faka * (_cou[1][15][0] - _fakaac * _cou[1][15][1]);
                _cou[19][15][0] = _wmp[2] * _cou[4][15][1] + _fakac * _cou[4][7][1];
                _cou[13][16][0] = _wmp[0] * _cou[4][16][1] + _faka * (_cou[2][16][0] - _fakaac * _cou[2][16][1]) + _fakac * _cou[4][9][1];
                _cou[14][16][0] = _wmp[1] * _cou[4][16][1] + _faka * (_cou[1][16][0] - _fakaac * _cou[1][16][1]);
                _cou[19][16][0] = _wmp[2] * _cou[4][16][1] + _fakac2 * _cou[4][5][1];
                _cou[13][17][0] = _wmp[0] * _cou[4][17][1] + _faka * (_cou[2][17][0] - _fakaac * _cou[2][17][1]);
                _cou[14][17][0] = _wmp[1] * _cou[4][17][1] + _faka * (_cou[1][17][0] - _fakaac * _cou[1][17][1]) + _fakac2 * _cou[4][6][1];
                _cou[19][17][0] = _wmp[2] * _cou[4][17][1] + _fakac * _cou[4][8][1];
                _cou[13][18][0] = _wmp[0] * _cou[4][18][1] + _faka * (_cou[2][18][0] - _fakaac * _cou[2][18][1]);
                _cou[14][18][0] = _wmp[1] * _cou[4][18][1] + _faka * (_cou[1][18][0] - _fakaac * _cou[1][18][1]) + _fakac * _cou[4][9][1];
                _cou[19][18][0] = _wmp[2] * _cou[4][18][1] + _fakac2 * _cou[4][6][1];
                _cou[13][19][0] = _wmp[0] * _cou[4][19][1] + _faka * (_cou[2][19][0] - _fakaac * _cou[2][19][1]) + _fakac * _cou[4][6][1];
                _cou[14][19][0] = _wmp[1] * _cou[4][19][1] + _faka * (_cou[1][19][0] - _fakaac * _cou[1][19][1]) + _fakac * _cou[4][5][1];
                _cou[19][19][0] = _wmp[2] * _cou[4][19][1] + _fakac * _cou[4][4][1];
                _cou[15][10][0] = _wmp[0] * _cou[5][10][1] + _faka * (_cou[3][10][0] - _fakaac * _cou[3][10][1]) + _fakac3 * _cou[5][7][1];
                _cou[16][10][0] = _wmp[2] * _cou[5][10][1] + _faka * (_cou[1][10][0] - _fakaac * _cou[1][10][1]);
                _cou[15][11][0] = _wmp[0] * _cou[5][11][1] + _faka * (_cou[3][11][0] - _fakaac * _cou[3][11][1]);
                _cou[16][11][0] = _wmp[2] * _cou[5][11][1] + _faka * (_cou[1][11][0] - _fakaac * _cou[1][11][1]);
                _cou[15][12][0] = _wmp[0] * _cou[5][12][1] + _faka * (_cou[3][12][0] - _fakaac * _cou[3][12][1]);
                _cou[16][12][0] = _wmp[2] * _cou[5][12][1] + _faka * (_cou[1][12][0] - _fakaac * _cou[1][12][1]) + _fakac3 * _cou[5][9][1];
                _cou[15][13][0] = _wmp[0] * _cou[5][13][1] + _faka * (_cou[3][13][0] - _fakaac * _cou[3][13][1]) + _fakac2 * _cou[5][4][1];
                _cou[16][13][0] = _wmp[2] * _cou[5][13][1] + _faka * (_cou[1][13][0] - _fakaac * _cou[1][13][1]);
                _cou[15][14][0] = _wmp[0] * _cou[5][14][1] + _faka * (_cou[3][14][0] - _fakaac * _cou[3][14][1]) + _fakac * _cou[5][8][1];
                _cou[16][14][0] = _wmp[2] * _cou[5][14][1] + _faka * (_cou[1][14][0] - _fakaac * _cou[1][14][1]);
                _cou[15][15][0] = _wmp[0] * _cou[5][15][1] + _faka * (_cou[3][15][0] - _fakaac * _cou[3][15][1]) + _fakac2 * _cou[5][5][1];
                _cou[16][15][0] = _wmp[2] * _cou[5][15][1] + _faka * (_cou[1][15][0] - _fakaac * _cou[1][15][1]) + _fakac * _cou[5][7][1];
                _cou[15][16][0] = _wmp[0] * _cou[5][16][1] + _faka * (_cou[3][16][0] - _fakaac * _cou[3][16][1]) + _fakac * _cou[5][9][1];
                _cou[16][16][0] = _wmp[2] * _cou[5][16][1] + _faka * (_cou[1][16][0] - _fakaac * _cou[1][16][1]) + _fakac2 * _cou[5][5][1];
                _cou[15][17][0] = _wmp[0] * _cou[5][17][1] + _faka * (_cou[3][17][0] - _fakaac * _cou[3][17][1]);
                _cou[16][17][0] = _wmp[2] * _cou[5][17][1] + _faka * (_cou[1][17][0] - _fakaac * _cou[1][17][1]) + _fakac * _cou[5][8][1];
                _cou[15][18][0] = _wmp[0] * _cou[5][18][1] + _faka * (_cou[3][18][0] - _fakaac * _cou[3][18][1]);
                _cou[16][18][0] = _wmp[2] * _cou[5][18][1] + _faka * (_cou[1][18][0] - _fakaac * _cou[1][18][1]) + _fakac2 * _cou[5][6][1];
                _cou[15][19][0] = _wmp[0] * _cou[5][19][1] + _faka * (_cou[3][19][0] - _fakaac * _cou[3][19][1]) + _fakac * _cou[5][6][1];
                _cou[16][19][0] = _wmp[2] * _cou[5][19][1] + _faka * (_cou[1][19][0] - _fakaac * _cou[1][19][1]) + _fakac * _cou[5][4][1];
                _cou[17][10][0] = _wmp[1] * _cou[6][10][1] + _faka * (_cou[3][10][0] - _fakaac * _cou[3][10][1]);
                _cou[18][10][0] = _wmp[2] * _cou[6][10][1] + _faka * (_cou[2][10][0] - _fakaac * _cou[2][10][1]);
                _cou[17][11][0] = _wmp[1] * _cou[6][11][1] + _faka * (_cou[3][11][0] - _fakaac * _cou[3][11][1]) + _fakac3 * _cou[6][8][1];
                _cou[18][11][0] = _wmp[2] * _cou[6][11][1] + _faka * (_cou[2][11][0] - _fakaac * _cou[2][11][1]);
                _cou[17][12][0] = _wmp[1] * _cou[6][12][1] + _faka * (_cou[3][12][0] - _fakaac * _cou[3][12][1]);
                _cou[18][12][0] = _wmp[2] * _cou[6][12][1] + _faka * (_cou[2][12][0] - _fakaac * _cou[2][12][1]) + _fakac3 * _cou[6][9][1];
                _cou[17][13][0] = _wmp[1] * _cou[6][13][1] + _faka * (_cou[3][13][0] - _fakaac * _cou[3][13][1]) + _fakac * _cou[6][7][1];
                _cou[18][13][0] = _wmp[2] * _cou[6][13][1] + _faka * (_cou[2][13][0] - _fakaac * _cou[2][13][1]);
                _cou[17][14][0] = _wmp[1] * _cou[6][14][1] + _faka * (_cou[3][14][0] - _fakaac * _cou[3][14][1]) + _fakac2 * _cou[6][4][1];
                _cou[18][14][0] = _wmp[2] * _cou[6][14][1] + _faka * (_cou[2][14][0] - _fakaac * _cou[2][14][1]);
                _cou[17][15][0] = _wmp[1] * _cou[6][15][1] + _faka * (_cou[3][15][0] - _fakaac * _cou[3][15][1]);
                _cou[18][15][0] = _wmp[2] * _cou[6][15][1] + _faka * (_cou[2][15][0] - _fakaac * _cou[2][15][1]) + _fakac * _cou[6][7][1];
                _cou[17][16][0] = _wmp[1] * _cou[6][16][1] + _faka * (_cou[3][16][0] - _fakaac * _cou[3][16][1]);
                _cou[18][16][0] = _wmp[2] * _cou[6][16][1] + _faka * (_cou[2][16][0] - _fakaac * _cou[2][16][1]) + _fakac2 * _cou[6][5][1];
                _cou[17][17][0] = _wmp[1] * _cou[6][17][1] + _faka * (_cou[3][17][0] - _fakaac * _cou[3][17][1]) + _fakac2 * _cou[6][6][1];
                _cou[18][17][0] = _wmp[2] * _cou[6][17][1] + _faka * (_cou[2][17][0] - _fakaac * _cou[2][17][1]) + _fakac * _cou[6][8][1];
                _cou[17][18][0] = _wmp[1] * _cou[6][18][1] + _faka * (_cou[3][18][0] - _fakaac * _cou[3][18][1]) + _fakac * _cou[6][9][1];
                _cou[18][18][0] = _wmp[2] * _cou[6][18][1] + _faka * (_cou[2][18][0] - _fakaac * _cou[2][18][1]) + _fakac2 * _cou[6][6][1];
                _cou[17][19][0] = _wmp[1] * _cou[6][19][1] + _faka * (_cou[3][19][0] - _fakaac * _cou[3][19][1]) + _fakac * _cou[6][5][1];
                _cou[18][19][0] = _wmp[2] * _cou[6][19][1] + _faka * (_cou[2][19][0] - _fakaac * _cou[2][19][1]) + _fakac * _cou[6][4][1];
                _cou[10][10][0] = _wmp[0] * _cou[7][10][1] + _faka2 * (_cou[1][10][0] - _fakaac * _cou[1][10][1]) + _fakac3 * _cou[7][7][1];
                _cou[10][11][0] = _wmp[0] * _cou[7][11][1] + _faka2 * (_cou[1][11][0] - _fakaac * _cou[1][11][1]);
                _cou[10][12][0] = _wmp[0] * _cou[7][12][1] + _faka2 * (_cou[1][12][0] - _fakaac * _cou[1][12][1]);
                _cou[10][13][0] = _wmp[0] * _cou[7][13][1] + _faka2 * (_cou[1][13][0] - _fakaac * _cou[1][13][1]) + _fakac2 * _cou[7][4][1];
                _cou[10][14][0] = _wmp[0] * _cou[7][14][1] + _faka2 * (_cou[1][14][0] - _fakaac * _cou[1][14][1]) + _fakac * _cou[7][8][1];
                _cou[10][15][0] = _wmp[0] * _cou[7][15][1] + _faka2 * (_cou[1][15][0] - _fakaac * _cou[1][15][1]) + _fakac2 * _cou[7][5][1];
                _cou[10][16][0] = _wmp[0] * _cou[7][16][1] + _faka2 * (_cou[1][16][0] - _fakaac * _cou[1][16][1]) + _fakac * _cou[7][9][1];
                _cou[10][17][0] = _wmp[0] * _cou[7][17][1] + _faka2 * (_cou[1][17][0] - _fakaac * _cou[1][17][1]);
                _cou[10][18][0] = _wmp[0] * _cou[7][18][1] + _faka2 * (_cou[1][18][0] - _fakaac * _cou[1][18][1]);
                _cou[10][19][0] = _wmp[0] * _cou[7][19][1] + _faka2 * (_cou[1][19][0] - _fakaac * _cou[1][19][1]) + _fakac * _cou[7][6][1];
                _cou[11][10][0] = _wmp[1] * _cou[8][10][1] + _faka2 * (_cou[2][10][0] - _fakaac * _cou[2][10][1]);
                _cou[11][11][0] = _wmp[1] * _cou[8][11][1] + _faka2 * (_cou[2][11][0] - _fakaac * _cou[2][11][1]) + _fakac3 * _cou[8][8][1];
                _cou[11][12][0] = _wmp[1] * _cou[8][12][1] + _faka2 * (_cou[2][12][0] - _fakaac * _cou[2][12][1]);
                _cou[11][13][0] = _wmp[1] * _cou[8][13][1] + _faka2 * (_cou[2][13][0] - _fakaac * _cou[2][13][1]) + _fakac * _cou[8][7][1];
                _cou[11][14][0] = _wmp[1] * _cou[8][14][1] + _faka2 * (_cou[2][14][0] - _fakaac * _cou[2][14][1]) + _fakac2 * _cou[8][4][1];
                _cou[11][15][0] = _wmp[1] * _cou[8][15][1] + _faka2 * (_cou[2][15][0] - _fakaac * _cou[2][15][1]);
                _cou[11][16][0] = _wmp[1] * _cou[8][16][1] + _faka2 * (_cou[2][16][0] - _fakaac * _cou[2][16][1]);
                _cou[11][17][0] = _wmp[1] * _cou[8][17][1] + _faka2 * (_cou[2][17][0] - _fakaac * _cou[2][17][1]) + _fakac2 * _cou[8][6][1];
                _cou[11][18][0] = _wmp[1] * _cou[8][18][1] + _faka2 * (_cou[2][18][0] - _fakaac * _cou[2][18][1]) + _fakac * _cou[8][9][1];
                _cou[11][19][0] = _wmp[1] * _cou[8][19][1] + _faka2 * (_cou[2][19][0] - _fakaac * _cou[2][19][1]) + _fakac * _cou[8][5][1];
                _cou[12][10][0] = _wmp[2] * _cou[9][10][1] + _faka2 * (_cou[3][10][0] - _fakaac * _cou[3][10][1]);
                _cou[12][11][0] = _wmp[2] * _cou[9][11][1] + _faka2 * (_cou[3][11][0] - _fakaac * _cou[3][11][1]);
                _cou[12][12][0] = _wmp[2] * _cou[9][12][1] + _faka2 * (_cou[3][12][0] - _fakaac * _cou[3][12][1]) + _fakac3 * _cou[9][9][1];
                _cou[12][13][0] = _wmp[2] * _cou[9][13][1] + _faka2 * (_cou[3][13][0] - _fakaac * _cou[3][13][1]);
                _cou[12][14][0] = _wmp[2] * _cou[9][14][1] + _faka2 * (_cou[3][14][0] - _fakaac * _cou[3][14][1]);
                _cou[12][15][0] = _wmp[2] * _cou[9][15][1] + _faka2 * (_cou[3][15][0] - _fakaac * _cou[3][15][1]) + _fakac * _cou[9][7][1];
                _cou[12][16][0] = _wmp[2] * _cou[9][16][1] + _faka2 * (_cou[3][16][0] - _fakaac * _cou[3][16][1]) + _fakac2 * _cou[9][5][1];
                _cou[12][17][0] = _wmp[2] * _cou[9][17][1] + _faka2 * (_cou[3][17][0] - _fakaac * _cou[3][17][1]) + _fakac * _cou[9][8][1];
                _cou[12][18][0] = _wmp[2] * _cou[9][18][1] + _faka2 * (_cou[3][18][0] - _fakaac * _cou[3][18][1]) + _fakac2 * _cou[9][6][1];
                _cou[12][19][0] = _wmp[2] * _cou[9][19][1] + _faka2 * (_cou[3][19][0] - _fakaac * _cou[3][19][1]) + _fakac * _cou[9][4][1];

            }


          
    if(_lmax_row>3 ||_lmax_col>3) {
            //cout << "g"<< endl;  
           FillgOrbitals(_wmp, _wmq, _cou, _decay_row, _decay_col,_lmax_row,_lmax_col);
          // cout << "g done"<< endl;  
    }
 
 
 
         
            // normalization and cartesian -> spherical factors
            int _ntrafo_row = _shell_row->getNumFunc() + _shell_row->getOffset();
            int _ntrafo_col = _shell_col->getNumFunc() + _shell_col->getOffset();

            //cout << " _ntrafo_row " << _ntrafo_row << ":" << _shell_row->getType() << endl;
            //cout << " _ntrafo_col " << _ntrafo_col << ":" << _shell_col->getType() << endl;
            ub::matrix<double> _trafo_row = ub::zero_matrix<double>(_ntrafo_row, _nrows);
            ub::matrix<double> _trafo_col = ub::zero_matrix<double>(_ntrafo_col, _ncols);

            // get transformation matrices including contraction coefficients
          std::vector<double> _contractions_row = (*itr)->contraction;
          std::vector<double> _contractions_col = (*itc)->contraction;

          this->getTrafo( _trafo_row, _lmax_row, _decay_row, _contractions_row);
          this->getTrafo( _trafo_col, _lmax_col, _decay_col, _contractions_col);

            // put _cou[i][j][0] into ublas matrix
            ub::matrix<double> _coumat = ub::zero_matrix<double>(_nrows, _ncols);
            for (unsigned i = 0; i < _coumat.size1(); i++) {
                for (unsigned j = 0; j < _coumat.size2(); j++) {
                    _coumat(i, j) = _cou[i][j][0];
                }
            }

            ub::matrix<double> _cou_tmp = ub::prod(_trafo_row, _coumat);
            ub::matrix<double> _trafo_col_tposed = ub::trans(_trafo_col);
            ub::matrix<double> _cou_sph = ub::prod(_cou_tmp, _trafo_col_tposed);
            // save to _matrix
            for (unsigned i = 0; i < _matrix.size1(); i++) {
                for (unsigned j = 0; j < _matrix.size2(); j++) {
                    _matrix(i, j) += _cou_sph(i + _shell_row->getOffset(), j + _shell_col->getOffset());
                }
            }


            //_ol.clear();

                } // _shell_col Gaussians
            } // _shell_row Gaussians
            
        }    
    
    
    
   
    
    void AOCoulomb::Symmetrize( AOOverlap& _gwoverlap, AOBasis& gwbasis, AOOverlap& _gwoverlap_inverse, AOOverlap& _gwoverlap_cholesky_inverse){
       
        //Logger* pLog = opThread->getLogger();
             
        if ( gwbasis._is_stable ){
            
            // get inverse of _aooverlap
            // Inversion of the matrix using GSL (much faster than boost)
            ub::matrix<double> _overlap_copy = _gwoverlap._aomatrix;
            linalg_invert( _overlap_copy, _gwoverlap_inverse._aomatrix );
            // cout << TimeStamp() << " Inverted GW Overlap matrix " <<   endl;
            _overlap_copy.resize(0,0);
            //_gwoverlap_inverse.Print( "S^-1" );

            // getting Cholesky decomposition of AOOverlap matrix
            AOOverlap _gwoverlap_cholesky;
            // make copy of _gwoverlap, because matrix is overwritten in GSL
            _gwoverlap_cholesky._aomatrix = _gwoverlap._aomatrix;
            linalg_cholesky_decompose( _gwoverlap_cholesky._aomatrix );
            // cout << TimeStamp() << " Calculated Cholesky decomposition of GW Overlap matrix " <<  endl;
            //_gwoverlap_cholesky.Print( "ChoS" );

            // remove L^T from Cholesky
            for (unsigned i =0; i < _gwoverlap_cholesky._aomatrix.size1(); i++ ){
                for (unsigned j = i+1; j < _gwoverlap_cholesky._aomatrix.size1(); j++ ){
                    _gwoverlap_cholesky._aomatrix(i,j) = 0.0;
                }
            }
            //_gwoverlap_cholesky.Print( "ChoS_zeroed" );

            // invert L to get L^-1
            //AOOverlap _gwoverlap_cholesky_inverse;
            _gwoverlap_cholesky_inverse.Initialize(gwbasis._AOBasisSize);
            _overlap_copy = _gwoverlap_cholesky._aomatrix;
            linalg_invert( _overlap_copy , _gwoverlap_cholesky_inverse._aomatrix );
            // cout << TimeStamp() << " Inverted Cholesky of GW Overlap " <<  endl;
            //_gwoverlap_cholesky_inverse.Print( "L^-1" );
            _overlap_copy.resize(0,0);

   
            
            
            // calculate V' = L^-1 V (L^-1)^T
            ub::matrix<double> _temp ( gwbasis._AOBasisSize, gwbasis._AOBasisSize);
            //_temp = ub::prod( _gwoverlap_cholesky_inverse._aomatrix , _gwcoulomb._aomatrix );
            _temp = ub::prod( _gwoverlap_cholesky_inverse._aomatrix , this->_aomatrix );


            // boost standard, nesting prod and trans is superslow
            //this->_aomatrix = ub::prod( _temp, ub::trans(_gwoverlap_cholesky_inverse._aomatrix ));
            ub::matrix<double> _gwoverlap_cholesky_inverse_transposed = ub::trans(_gwoverlap_cholesky_inverse._aomatrix );
            this->_aomatrix = ub::prod( _temp, _gwoverlap_cholesky_inverse_transposed);
            
            // cout << TimeStamp() << " Multiplied GW Coulomb with L^-1 and (L^-1)^T " <<  endl;
            // this->Print( "CouSu" );

            ub::vector<double>                  _eigenvalues;
            ub::matrix<double>                  _eigenvectors;

            // get eigenvectors and eigenvalues of V'
            // LA_Eigenvalues( this->_aomatrix , _eigenvalues, _eigenvectors);
            linalg_eigenvalues( this->_aomatrix , _eigenvalues, _eigenvectors);
            // calc sqrt(V')
            _temp.clear();
            //cout << _eigenvalues<<endl;
            for ( int i = 0; i  < gwbasis._AOBasisSize; i++ ){

                if ( _eigenvalues(i) < 0.0 ) {
                    cout << "Warning: negative eigenvalue!" << endl;
                    _eigenvalues(i) = 0.0;
                }
                for ( int j = 0; j < gwbasis._AOBasisSize; j++){
                    _temp(i,j) = _eigenvectors(j,i) * sqrt(_eigenvalues(i));
                }
            }
            
            this->_aomatrix = ub::prod(_eigenvectors, _temp);
            // cout << TimeStamp() << " Calculated sqrt(V') matrix " <<  endl;
            // this->Print( "CouEV" );

            // multiply with L from the left and L+ from the right
            _temp = ub::prod( _gwoverlap_cholesky._aomatrix , this->_aomatrix );
            
            
            ub::matrix<double> _gwoverlap_cholesky_transposed = ub::trans( _gwoverlap_cholesky._aomatrix );
            this->_aomatrix = ub::prod( _temp ,_gwoverlap_cholesky_transposed);
            
            
            // cout << TimeStamp() << " Coulomb matrix sqrt'ed " <<  endl;
            // this->Print( "CouSqrt" );
            // multiply _gwcoulomb with _gwoverlap_inverse
            this->_aomatrix = ub::prod( this->_aomatrix , _gwoverlap_inverse._aomatrix );
            // cout << TimeStamp() << " Final Coulomb matrix  " <<  endl;
            // this->Print( " COUfinal ");
        }
        
    }
    
    
    
    void AOCoulomb::Symmetrize_DFT( AOOverlap& _gwoverlap, AOBasis& gwbasis, AOOverlap& _gwoverlap_inverse, AOOverlap& _gwoverlap_cholesky_inverse){
        
        //Logger* pLog = opThread->getLogger();
             
        if ( gwbasis._is_stable ){
            cout << "doing something " << endl;
            // get inverse of _aooverlap
            
            
            // Inversion of the matrix using GSL (much faster than boost)
            ub::matrix<double> _overlap_copy = _gwoverlap._aomatrix;
            linalg_invert( _overlap_copy, _gwoverlap_inverse._aomatrix );
            // cout << TimeStamp() << " Inverted GW Overlap matrix " <<   endl;
            _overlap_copy.resize(0,0);
            //_gwoverlap_inverse.Print( "S^-1" );

            ub::matrix<double> Vcou_backup = this->_aomatrix;
            
            /******* SKIPPING CHOLESKY STUFF ********/
            /*
            // getting Cholesky decomposition of AOOverlap matrix
            AOOverlap _gwoverlap_cholesky;
            // make copy of _gwoverlap, because matrix is overwritten in GSL
            _gwoverlap_cholesky._aomatrix = _gwoverlap._aomatrix;
            linalg_cholesky_decompose( _gwoverlap_cholesky._aomatrix );
            // cout << TimeStamp() << " Calculated Cholesky decomposition of GW Overlap matrix " <<  endl;
            //_gwoverlap_cholesky.Print( "ChoS" );

            // remove L^T from Cholesky
            for (int i =0; i < _gwoverlap_cholesky._aomatrix.size1(); i++ ){
                for (int j = i+1; j < _gwoverlap_cholesky._aomatrix.size1(); j++ ){
                    _gwoverlap_cholesky._aomatrix(i,j) = 0.0;
                }
            }
            //_gwoverlap_cholesky.Print( "ChoS_zeroed" );

            // invert L to get L^-1
            //AOOverlap _gwoverlap_cholesky_inverse;
            _gwoverlap_cholesky_inverse.Initialize(gwbasis._AOBasisSize);
            _overlap_copy = _gwoverlap_cholesky._aomatrix;
            linalg_invert( _overlap_copy , _gwoverlap_cholesky_inverse._aomatrix );
            // cout << TimeStamp() << " Inverted Cholesky of GW Overlap " <<  endl;
            //_gwoverlap_cholesky_inverse.Print( "L^-1" );
            _overlap_copy.resize(0,0);

   
            
            
            // calculate V' = L^-1 V (L^-1)^T
            ub::matrix<double> _temp ( gwbasis._AOBasisSize, gwbasis._AOBasisSize);
            //_temp = ub::prod( _gwoverlap_cholesky_inverse._aomatrix , _gwcoulomb._aomatrix );
            _temp = ub::prod( _gwoverlap_cholesky_inverse._aomatrix , this->_aomatrix );


            // boost standard, nesting prod and trans is superslow
            //this->_aomatrix = ub::prod( _temp, ub::trans(_gwoverlap_cholesky_inverse._aomatrix ));
            ub::matrix<double> _gwoverlap_cholesky_inverse_transposed = ub::trans(_gwoverlap_cholesky_inverse._aomatrix );
            this->_aomatrix = ub::prod( _temp, _gwoverlap_cholesky_inverse_transposed);
            
            // cout << TimeStamp() << " Multiplied GW Coulomb with L^-1 and (L^-1)^T " <<  endl;
            // this->Print( "CouSu" );
*/
            
            
            ub::vector<double>                  _eigenvalues;
            ub::matrix<double>                  _eigenvectors;

            // get eigenvectors and eigenvalues of V'
            // LA_Eigenvalues( this->_aomatrix , _eigenvalues, _eigenvectors);
            linalg_eigenvalues( this->_aomatrix , _eigenvalues, _eigenvectors);
            // calc sqrt(V')
            ub::matrix<double> _temp ( gwbasis._AOBasisSize, gwbasis._AOBasisSize);
            _temp.clear();
            for ( int i = 0; i  < gwbasis._AOBasisSize; i++ ){

                if ( _eigenvalues(i) < 0.0 ) {
                    cout << "Warning: negative eigenvalue!" << endl;
                    _eigenvalues(i) = 0.0;
                }
                for ( int j = 0; j < gwbasis._AOBasisSize; j++){
                    _temp(i,j) = _eigenvectors(j,i) * sqrt(_eigenvalues(i));
                }
            }
            
            this->_aomatrix = ub::prod(_eigenvectors, _temp);
            // cout << TimeStamp() << " Calculated sqrt(V') matrix " <<  endl;
            // this->Print( "CouEV" );

            /****** SKIPPING AGAIN ******/
            /*
            // multiply with L from the left and L+ from the right
            _temp = ub::prod( _gwoverlap_cholesky._aomatrix , this->_aomatrix );
            
            
            ub::matrix<double> _gwoverlap_cholesky_transposed = ub::trans( _gwoverlap_cholesky._aomatrix );
            this->_aomatrix = ub::prod( _temp ,_gwoverlap_cholesky_transposed);
            
            
            // cout << TimeStamp() << " Coulomb matrix sqrt'ed " <<  endl;
            // this->Print( "CouSqrt" );
             * 
             * 
             */
            ub::matrix<double> Vcou =ub::prod(this->_aomatrix,this->_aomatrix);
            for ( int i =0; i < this->Dimension(); i++ ){
            for ( int j =0; j < this->Dimension(); j++ ){
                
                cout <<  i << " : "  << j << " = " << Vcou(i,j) << " vs " << Vcou_backup(i,j) << endl;
                
            }                
            }
            
            
            
            // multiply _gwcoulomb with _gwoverlap_inverse
            this->_aomatrix = ub::prod( this->_aomatrix , _gwoverlap_inverse._aomatrix );
            // cout << TimeStamp() << " Final Coulomb matrix  " <<  endl;
            // this->Print( " COUfinal ");
        }
        
    }
    
}}

