/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICEN_olE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/aomatrix.h>

#include <votca/ctp/aobasis.h>
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

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;

    
   int AOCoulomb::getExtraBlockSize(int _lmax_row, int _lmax_col){
        int _block_size = _lmax_col + _lmax_row +1;
        return _block_size;
    }

        void AOCoulomb::FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, AOBasis* ecp) {

            // shell info, only lmax tells how far to go
            int _lmax_row = _shell_row->getLmax();
            int _lmax_col = _shell_col->getLmax();

            // set size of internal block for recursion
            int _nrows = this->getBlockSize(_lmax_row);
            int _ncols = this->getBlockSize(_lmax_col);
            int _nextra = this->getExtraBlockSize(_lmax_row, _lmax_col);
            int _l_sum = _lmax_row + _lmax_col;
            int _ma_dim = this->getBlockSize(_l_sum);
            
            // get shell positions
            const vec& _pos_row = _shell_row->getPos();
            const vec& _pos_col = _shell_col->getPos();
            const vec _diff = _pos_row - _pos_col;
            double _distsq = (_diff.getX() * _diff.getX()) + (_diff.getY() * _diff.getY()) + (_diff.getZ() * _diff.getZ());
            
            const double pi = boost::math::constants::pi<double>();
             // some helpers
            vector<double> _wmp;
            vector<double> _wmq;
            
            
         

           
            
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
        // iterate over Gaussians in this _shell_row
            for ( GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
                const double& _decay_row = (*itr)->decay;
            
                for ( GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
                    
                     // get decay constants 
                        const double& _decay_col = (*itc)->decay;

                        double _fakac = 0.5 / (_decay_row + _decay_col);
                        double _fakac2 = 1. / (_decay_row + _decay_col);
                        // check if distance between postions is big, then skip step   
                        //double _exparg = _fakac2 * _decay_row * _decay_col *_distsq;
                        // if ( _exparg > 30.0 ) { continue; } //!!!!!CUTOFF not applicable to AOCoulomb (at least not like this...)
                    
                                    // get a multi dimensional array
                         typedef boost::multi_array<double, 3> ma_type;
                         //ma_type _cou(boost::extents[_nrows][_ncols][_nextra]);
                         ma_type _cou(boost::extents[_ma_dim][_ma_dim][_nextra]);
                         typedef ma_type::index index;
                                    // initialize to zero_cou[0][0][i] 
                           //       for(index i = 0; i != _nrows; ++i) {
                           //  for(index j = 0; j != _ncols; ++j){
                           for (index i = 0; i != _ma_dim; ++i) {
                               for (index j = 0; j != _ma_dim; ++j) {
                                   for (index k = 0; k != _nextra; ++k) {
                                       _cou[i][j][k] = 0.0;
                                   }
                               }
                           }

                       

            // some helpers
            const double _faka = 0.5 / _decay_row;
            const double _faka2 = 2.0 * _faka;
            //const double _faka3 = 3.0 * _faka;
            const double _fakaca = _decay_row / (_decay_row + _decay_col);
            const double _fakaac = _decay_row / (_decay_row + _decay_col);
            
            const double _fakac3 = 3.0 * _fakac;
            //const double _fakac4 = 4.0 * _fakac;
            const double _fakc = 0.5 / _decay_col;
            const double _fakc2 = 2.0 * _fakc;
            //const double _fakc3 = 3.0 * _fakc;
            //const double _fakca = _fakac;
            //const double _fakca2 = _fakac2;
            //const double _fakca3 = _fakac3;
            //const double _fakca4 = _fakac4;

            _wmp.resize(3);
            _wmq.resize(3);
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

            vector<double> _FmT(_nextra, 0.0); // that size needs to be checked!
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

            /* The following is for all elements involving g-orbitals, but I'm tired 
             * of this typing, can be done at some other point, if g-orbitals are 
             * actually needed....


       if(abs(lmax1)+abs(lmax2).ge.4) then
          _cou[1][0][3) = _wmp[0]*_cou[0][0][4)
          _cou[2][0][3) = _wmp[1]*_cou[0][0][4)
          _cou[3][0][3) = _wmp[2]*_cou[0][0][4)

          _cou[7][0][2] = _wmp[0]*_cou[1][0][3) & 
               +    faka * (_cou[0][0][2) - fakaac * _cou[0][0][3) )
          _cou[4][0][2] = _wmp[1]*_cou[1][0][3)
          _cou[5][0][2) = _wmp[2]*_cou[1][0][3)
          _cou[8][0][2) = _wmp[1]*_cou[2][0][3) & 
               +    faka * (_cou[0][0][2) - fakaac * _cou[0][0][3) )
          _cou[6][0][2) = _wmp[2]*_cou[2][0][3)
          _cou[9][0][2) = _wmp[2]*_cou[3][0][3) & 
               +    faka * (_cou[0][0][2) - fakaac * _cou[0][0][3) )

  
 
          _cou[10][0][1) = _wmp[0]*_cou[7][0][2) & 
               +    faka2 * (_cou[1][0][1) - fakaac * _cou[1][0][2) )
          _cou[11][0][1) = _wmp[1]*_cou[8][0][2) & 
               +    faka2 * (_cou[2][0][1) - fakaac * _cou[2][0][2) )
          _cou[12][0][1) = _wmp[2]*_cou[9][0][2) & 
               +    faka2 * (_cou[3][0][1) - fakaac * _cou[3][0][2) )
          _cou[13][0][1) = _wmp[0]*_cou[4][0][2) & 
               +    faka * (_cou[2][0][1) - fakaac * _cou[2][0][2) )
          _cou[14][0][1) = _wmp[1]*_cou[4][0][2) & 
               +    faka * (_cou[1][0][1) - fakaac * _cou[1][0][2) )
          _cou[19][0][1) = _wmp[2]*_cou[4][0][2)
          _cou[15][0][1) = _wmp[0]*_cou[5][0][2) & 
               +    faka * (_cou[3][0][1) - fakaac * _cou[3][0][2) )
          _cou[16][0][1) = _wmp[2]*_cou[5][0][2) & 
               +    faka * (_cou[1][0][1) - fakaac * _cou[1][0][2) )
          _cou[17][0][1) = _wmp[1]*_cou[6][0][2) & 
               +    faka * (_cou[3][0][1) - fakaac * _cou[3][0][2) )
          _cou[18][0][1) = _wmp[2]*_cou[6][0][2) & 
               +    faka * (_cou[2][0][1) - fakaac * _cou[2][0][2) )

       endif

       if(abs(lmax1)+abs(lmax2).ge.5) then
          _cou[1][0][4) = _wmp[0]*_cou[0][0][5)
          _cou[2][0][4) = _wmp[1]*_cou[0][0][5)
          _cou[3][0][4) = _wmp[2]*_cou[0][0][5)

          _cou[7][0][3) = _wmp[0]*_cou[1][0][4) &
               +    faka * (_cou[0][0][3) - fakaac * _cou[0][0][4) )
          _cou[4][0][3) = _wmp[1]*_cou[1][0][4)
          _cou[5][0][3) = _wmp[2]*_cou[1][0][4)
          _cou[8][0][3) = _wmp[1]*_cou[2][0][4) & 
               +    faka * (_cou[0][0][3) - fakaac * _cou[0][0][4) )
          _cou[6][0][3) = _wmp[2]*_cou[2][0][4)
          _cou[9][0][3) = _wmp[2]*_cou[3][0][4) & 
               +    faka * (_cou[0][0][3) - fakaac * _cou[0][0][4) )

          _cou[1][1][3) = _wmp[0]*_cou[0][1][4) + fakac * _cou[0][0][4)
          _cou[2][1][3) = _wmp[1]*_cou[0][1][4)
          _cou[3][1][3) = _wmp[2]*_cou[0][1][4)
          _cou[1][2][3) = _wmp[0]*_cou[0][2][4)
          _cou[2][2][3) = _wmp[1]*_cou[0][2][4) + fakac * _cou[0][0][4)
          _cou[3][2][3) = _wmp[2]*_cou[0][2][4)
          _cou[1][3][3) = _wmp[0]*_cou[0][3][4)
          _cou[2][3][3) = _wmp[1]*_cou[0][3][4)
          _cou[3][3][3) = _wmp[2]*_cou[0][3][4) + fakac * _cou[0][0][4)




          _cou[13][0][2) = _wmp[0]*_cou[4][0][3) & 
               +    faka * (_cou[2][0][2) - fakaac * _cou[2][0][3) )
          _cou[14][0][2) = _wmp[1]*_cou[4][0][3) & 
               +    faka * (_cou[1][0][2) - fakaac * _cou[1][0][3) )
          _cou[19][0][2) = _wmp[2]*_cou[4][0][3)
          _cou[15][0][2) = _wmp[0]*_cou[5][0][3) & 
               +    faka * (_cou[3][0][2) - fakaac * _cou[3][0][3) )
          _cou[16][0][2) = _wmp[2]*_cou[5][0][3) & 
               +    faka * (_cou[1][0][2) - fakaac * _cou[1][0][3) )
          _cou[17][0][2) = _wmp[1]*_cou[6][0][3) & 
               +    faka * (_cou[3][0][2) - fakaac * _cou[3][0][3) )
          _cou[18][0][2) = _wmp[2]*_cou[6][0][3) & 
               +    faka * (_cou[2][0][2) - fakaac * _cou[2][0][3) )
          _cou[10][0][2) = _wmp[0]*_cou[7][0][3) & 
               +    faka2 * (_cou[1][0][2) - fakaac * _cou[1][0][3) )
          _cou[11][0][2) = _wmp[1]*_cou[8][0][3) & 
               +    faka2 * (_cou[2][0][2) - fakaac * _cou[2][0][3) )
          _cou[12][0][2) = _wmp[2]*_cou[9][0][3) & 
               +    faka2 * (_cou[3][0][2) - fakaac * _cou[3][0][3) )

          _cou[7][1][2) = _wmp[0]*_cou[1][1][3) & 
               +  faka * (_cou[0][1][2) - fakaac * _cou[0][1][3) ) + fakac * _cou[1][0][3)
          _cou[4][1][2) = _wmp[1]*_cou[1][1][3)
          _cou[5][1][2) = _wmp[2]*_cou[1][1][3)
          _cou[7][2][2) = _wmp[0]*_cou[1][2][3) & 
               +    faka * (_cou[0][2][2) - fakaac * _cou[0][2][3) )
          _cou[4][2][2) = _wmp[1]*_cou[1][2][3) + fakac * _cou[1][0][3)
          _cou[5][2][2) = _wmp[2]*_cou[1][2][3)
          _cou[7][3][2) = _wmp[0]*_cou[1][3][3) & 
               +    faka * (_cou[0][3][2) - fakaac * _cou[0][3][3) )
          _cou[4][3][2) = _wmp[1]*_cou[1][3][3)
          _cou[5][3][2) = _wmp[2]*_cou[1][3][3) + fakac * _cou[1][0][3)
          _cou[8][1][2) = _wmp[1]*_cou[2][1][3) & 
               +    faka * (_cou[0][1][2) - fakaac * _cou[0][1][3) )
          _cou[6][1][2) = _wmp[2]*_cou[2][1][3)
          _cou[8][2][2) = _wmp[1]*_cou[2][2][3) & 
               +  faka * (_cou[0][2][2) - fakaac * _cou[0][2][3) ) + fakac * _cou[2][0][3)
          _cou[6][2][2) = _wmp[2]*_cou[2][2][3)
          _cou[8][3][2) = _wmp[1]*_cou[2][3][3) & 
               +    faka * (_cou[0][3][2) - fakaac * _cou[0][3][3) )
          _cou[6][3][2) = _wmp[2]*_cou[2][3][3) + fakac * _cou[2][0][3)
          _cou[9][1][2) = _wmp[2]*_cou[3][1][3) & 
               +    faka * (_cou[0][1][2) - fakaac * _cou[0][1][3) )
          _cou[9][2][2) = _wmp[2]*_cou[3][2][3) & 
               +    faka * (_cou[0][2][2) - fakaac * _cou[0][2][3) )
          _cou[9][3][2) = _wmp[2]*_cou[3][3][3) & 
               +  faka * (_cou[0][3][2) - fakaac * _cou[0][3][3) ) + fakac * _cou[3][0][3)




          _cou[13][1][1) = _wmp[0]*_cou[4][1][2) & 
               +  faka * (_cou[2][1][1) - fakaac * _cou[2][1][2) ) + fakac * _cou[4][0][2)
          _cou[14][1][1) = _wmp[1]*_cou[4][1][2) & 
               +    faka * (_cou[1][1][1) - fakaac * _cou[1][1][2) )
          _cou[19][1][1) = _wmp[2]*_cou[4][1][2)
          _cou[13][2][1) = _wmp[0]*_cou[4][2][2) & 
               +    faka * (_cou[2][2][1) - fakaac * _cou[2][2][2) )
          _cou[14][2][1) = _wmp[1]*_cou[4][2][2) & 
               +  faka * (_cou[1][2][1) - fakaac * _cou[1][2][2) ) + fakac * _cou[4][0][2)
          _cou[19][2][1) = _wmp[2]*_cou[4][2][2)
          _cou[13][3][1) = _wmp[0]*_cou[4][3][2) & 
               +    faka * (_cou[2][3][1) - fakaac * _cou[2][3][2) )
          _cou[14][3][1) = _wmp[1]*_cou[4][3][2) & 
               +    faka * (_cou[1][3][1) - fakaac * _cou[1][3][2) )
          _cou[19][3][1) = _wmp[2]*_cou[4][3][2) + fakac * _cou[4][0][2)
          _cou[15][1][1) = _wmp[0]*_cou[5][1][2) & 
               +  faka * (_cou[3][1][1) - fakaac * _cou[3][1][2) ) + fakac * _cou[5][0][2)
          _cou[16][1][1) = _wmp[2]*_cou[5][1][2) & 
               +    faka * (_cou[1][1][1) - fakaac * _cou[1][1][2) )
          _cou[15][2][1) = _wmp[0]*_cou[5][2][2) & 
               +    faka * (_cou[3][2][1) - fakaac * _cou[3][2][2) )
          _cou[16][2][1) = _wmp[2]*_cou[5][2][2) & 
               +    faka * (_cou[1][2][1) - fakaac * _cou[1][2][2) )
          _cou[15][3][1) = _wmp[0]*_cou[5][3][2) & 
               +    faka * (_cou[3][3][1) - fakaac * _cou[3][3][2) )
          _cou[16][3][1) = _wmp[2]*_cou[5][3][2) & 
               +  faka * (_cou[1][3][1) - fakaac * _cou[1][3][2) ) + fakac * _cou[5][0][2)
          _cou[17][1][1) = _wmp[1]*_cou[6][1][2) & 
               +    faka * (_cou[3][1][1) - fakaac * _cou[3][1][2) )
          _cou[18][1][1) = _wmp[2]*_cou[6][1][2) & 
               +    faka * (_cou[2][1][1) - fakaac * _cou[2][1][2) )
          _cou[17][2][1) = _wmp[1]*_cou[6][2][2) & 
               +  faka * (_cou[3][2][1) - fakaac * _cou[3][2][2) ) + fakac * _cou[6][0][2)
          _cou[18][2][1) = _wmp[2]*_cou[6][2][2) & 
               +    faka * (_cou[2][2][1) - fakaac * _cou[2][2][2) )
          _cou[17][3][1) = _wmp[1]*_cou[6][3][2) & 
               +    faka * (_cou[3][3][1) - fakaac * _cou[3][3][2) )
          _cou[18][3][1) = _wmp[2]*_cou[6][3][2) & 
               +  faka * (_cou[2][3][1) - fakaac * _cou[2][3][2) ) + fakac * _cou[6][0][2)
          _cou[10][1][1) = _wmp[0]*_cou[7][1][2) & 
               + faka2 * (_cou[1][1][1) - fakaac * _cou[1][1][2) ) + fakac * _cou[7][0][2)
          _cou[10][2][1) = _wmp[0]*_cou[7][2][2) & 
               +    faka2 * (_cou[1][2][1) - fakaac * _cou[1][2][2) )
          _cou[10][3][1) = _wmp[0]*_cou[7][3][2) & 
               +    faka2 * (_cou[1][3][1) - fakaac * _cou[1][3][2) )
          _cou[11][1][1) = _wmp[1]*_cou[8][1][2) & 
               +    faka2 * (_cou[2][1][1) - fakaac * _cou[2][1][2) )
          _cou[11][2][1) = _wmp[1]*_cou[8][2][2) & 
               + faka2 * (_cou[2][2][1) - fakaac * _cou[2][2][2) ) + fakac * _cou[8][0][2)
          _cou[11][3][1) = _wmp[1]*_cou[8][3][2) & 
               +    faka2 * (_cou[2][3][1) - fakaac * _cou[2][3][2) )
          _cou[12][1][1) = _wmp[2]*_cou[9][1][2) & 
               +    faka2 * (_cou[3][1][1) - fakaac * _cou[3][1][2) )
          _cou[12][2][1) = _wmp[2]*_cou[9][2][2) & 
               +    faka2 * (_cou[3][2][1) - fakaac * _cou[3][2][2) )
          _cou[12][3][1) = _wmp[2]*_cou[9][3][2) & 
               + faka2 * (_cou[3][3][1) - fakaac * _cou[3][3][2) ) + fakac * _cou[9][0][2)



       endif

       if(abs(lmax1)+abs(lmax2).ge.6) then
          _cou[1][0][5) = _wmp[0]*_cou[0][0][6)
          _cou[2][0][5) = _wmp[1]*_cou[0][0][6)
          _cou[3][0][5) = _wmp[2]*_cou[0][0][6)

 
             * 
             * 
  
       
          _cou[1][1][4) = _wmp[0]*_cou[0][1][5) + fakac * _cou[0][0][5)
          _cou[2][1][4) = _wmp[1]*_cou[0][1][5)
          _cou[3][1][4) = _wmp[2]*_cou[0][1][5)
          _cou[1][2][4) = _wmp[0]*_cou[0][2][5)
          _cou[2][2][4) = _wmp[1]*_cou[0][2][5) + fakac * _cou[0][0][5)
          _cou[3][2][4) = _wmp[2]*_cou[0][2][5)
          _cou[1][3][4) = _wmp[0]*_cou[0][3][5)
          _cou[2][3][4) = _wmp[1]*_cou[0][3][5)
          _cou[3][3][4) = _wmp[2]*_cou[0][3][5) + fakac * _cou[0][0][5)
          _cou[7][0][4) = _wmp[0]*_cou[1][0][5) & 
               +    faka * (_cou[0][0][4) - fakaac * _cou[0][0][5) )
          _cou[4][0][4) = _wmp[1]*_cou[1][0][5)
          _cou[5][0][4) = _wmp[2]*_cou[1][0][5)
          _cou[8][0][4) = _wmp[1]*_cou[2][0][5) & 
               +    faka * (_cou[0][0][4) - fakaac * _cou[0][0][5) )
          _cou[6][0][4) = _wmp[2]*_cou[2][0][5)
          _cou[9][0][4) = _wmp[2]*_cou[3][0][5) & 
               +    faka * (_cou[0][0][4) - fakaac * _cou[0][0][5) )

             * 
             * 
             * 
             * 
 
          _cou[1][4][3) = _wmp[0]*_cou[0][4][4) + fakac * _cou[0][2][4)
          _cou[2][4][3) = _wmp[1]*_cou[0][4][4) + fakac * _cou[0][1][4)
          _cou[3][4][3) = _wmp[2]*_cou[0][4][4)
          _cou[1][5][3) = _wmp[0]*_cou[0][5][4) + fakac * _cou[0][3][4)
          _cou[2][5][3) = _wmp[1]*_cou[0][5][4)
          _cou[3][5][3) = _wmp[2]*_cou[0][5][4) + fakac * _cou[0][1][4)
          _cou[1][6][3) = _wmp[0]*_cou[0][6][4)
          _cou[2][6][3) = _wmp[1]*_cou[0][6][4) + fakac * _cou[0][3][4)
          _cou[3][6][3) = _wmp[2]*_cou[0][6][4) + fakac * _cou[0][2][4)
          _cou[1][7][3) = _wmp[0]*_cou[0][7][4) + fakac2 * _cou[0][1][4)
          _cou[2][7][3) = _wmp[1]*_cou[0][7][4)
          _cou[3][7][3) = _wmp[2]*_cou[0][7][4)
          _cou[1][8][3) = _wmp[0]*_cou[0][8][4)
          _cou[2][8][3) = _wmp[1]*_cou[0][8][4) + fakac2 * _cou[0][2][4)
          _cou[3][8][3) = _wmp[2]*_cou[0][8][4)
          _cou[1][9][3) = _wmp[0]*_cou[0][9][4)
          _cou[2][9][3) = _wmp[1]*_cou[0][9][4)
          _cou[3][9][3) = _wmp[2]*_cou[0][9][4) + fakac2 * _cou[0][3][4)
          _cou[7][1][3) = _wmp[0]*_cou[1][1][4) & 
               +  faka * (_cou[0][1][3) - fakaac * _cou[0][1][4) ) + fakac * _cou[1][0][4)
          _cou[4][1][3) = _wmp[1]*_cou[1][1][4)
          _cou[5][1][3) = _wmp[2]*_cou[1][1][4)
          _cou[7][2][3) = _wmp[0]*_cou[1][2][4) & 
               +    faka * (_cou[0][2][3) - fakaac * _cou[0][2][4) )
          _cou[4][2][3) = _wmp[1]*_cou[1][2][4) + fakac * _cou[1][0][4)
          _cou[5][2][3) = _wmp[2]*_cou[1][2][4)
          _cou[7][3][3) = _wmp[0]*_cou[1][3][4) & 
               +    faka * (_cou[0][3][3) - fakaac * _cou[0][3][4) )
          _cou[4][3][3) = _wmp[1]*_cou[1][3][4)
          _cou[5][3][3) = _wmp[2]*_cou[1][3][4) + fakac * _cou[1][0][4)
          _cou[8][1][3) = _wmp[1]*_cou[2][1][4) & 
               +    faka * (_cou[0][1][3) - fakaac * _cou[0][1][4) )
          _cou[6][1][3) = _wmp[2]*_cou[2][1][4)
          _cou[8][2][3) = _wmp[1]*_cou[2][2][4) & 
               +  faka * (_cou[0][2][3) - fakaac * _cou[0][2][4) ) + fakac * _cou[2][0][4)
          _cou[6][2][3) = _wmp[2]*_cou[2][2][4)
          _cou[8][3][3) = _wmp[1]*_cou[2][3][4) & 
               +    faka * (_cou[0][3][3) - fakaac * _cou[0][3][4) )
          _cou[6][3][3) = _wmp[2]*_cou[2][3][4) + fakac * _cou[2][0][4)
          _cou[9][1][3) = _wmp[2]*_cou[3][1][4) & 
               +    faka * (_cou[0][1][3) - fakaac * _cou[0][1][4) )
          _cou[9][2][3) = _wmp[2]*_cou[3][2][4) & 
               +    faka * (_cou[0][2][3) - fakaac * _cou[0][2][4) )
          _cou[9][3][3) = _wmp[2]*_cou[3][3][4) & 
               +  faka * (_cou[0][3][3) - fakaac * _cou[0][3][4) ) + fakac * _cou[3][0][4)
          _cou[13][0][3) = _wmp[0]*_cou[4][0][4) & 
               +    faka * (_cou[2][0][3) - fakaac * _cou[2][0][4) )
          _cou[14][0][3) = _wmp[1]*_cou[4][0][4) & 
               +    faka * (_cou[1][0][3) - fakaac * _cou[1][0][4) )
          _cou[19][0][3) = _wmp[2]*_cou[4][0][4)
          _cou[15][0][3) = _wmp[0]*_cou[5][0][4) & 
               +    faka * (_cou[3][0][3) - fakaac * _cou[3][0][4) )
          _cou[16][0][3) = _wmp[2]*_cou[5][0][4) & 
               +    faka * (_cou[1][0][3) - fakaac * _cou[1][0][4) )
          _cou[17][0][3) = _wmp[1]*_cou[6][0][4) & 
               +    faka * (_cou[3][0][3) - fakaac * _cou[3][0][4) )
          _cou[18][0][3) = _wmp[2]*_cou[6][0][4) & 
               +    faka * (_cou[2][0][3) - fakaac * _cou[2][0][4) )
          _cou[10][0][3) = _wmp[0]*_cou[7][0][4) & 
               +    faka2 * (_cou[1][0][3) - fakaac * _cou[1][0][4) )
          _cou[11][0][3) = _wmp[1]*_cou[8][0][4) & 
               +    faka2 * (_cou[2][0][3) - fakaac * _cou[2][0][4) )
          _cou[12][0][3) = _wmp[2]*_cou[9][0][4) & 
               +    faka2 * (_cou[3][0][3) - fakaac * _cou[3][0][4) )

          _cou[7][4][2) = _wmp[0]*_cou[1][4][3) & 
               +  faka * (_cou[0][4][2) - fakaac * _cou[0][4][3) ) + fakac * _cou[1][2][3)
          _cou[4][4][2) = _wmp[1]*_cou[1][4][3) + fakac * _cou[1][1][3)
          _cou[5][4][2) = _wmp[2]*_cou[1][4][3)
          _cou[7][5][2) = _wmp[0]*_cou[1][5][3) & 
               +  faka * (_cou[0][5][2) - fakaac * _cou[0][5][3) ) + fakac * _cou[1][3][3)
          _cou[4][5][2) = _wmp[1]*_cou[1][5][3)
          _cou[5][5][2) = _wmp[2]*_cou[1][5][3) + fakac * _cou[1][1][3)
          _cou[7][6][2) = _wmp[0]*_cou[1][6][3) & 
               +    faka * (_cou[0][6][2) - fakaac * _cou[0][6][3) )
          _cou[4][6][2) = _wmp[1]*_cou[1][6][3) + fakac * _cou[1][3][3)
          _cou[5][6][2) = _wmp[2]*_cou[1][6][3) + fakac * _cou[1][2][3)
          _cou[7][7][2) = _wmp[0]*_cou[1][7][3) & 
               + faka * (_cou[0][7][2) - fakaac * _cou[0][7][3) ) + fakac2 * _cou[1][1][3)
          _cou[4][7][2) = _wmp[1]*_cou[1][7][3)
          _cou[5][7][2) = _wmp[2]*_cou[1][7][3)
          _cou[7][8][2) = _wmp[0]*_cou[1][8][3) & 
               +    faka * (_cou[0][8][2) - fakaac * _cou[0][8][3) )
          _cou[4][8][2) = _wmp[1]*_cou[1][8][3) + fakac2 * _cou[1][2][3)
          _cou[5][8][2) = _wmp[2]*_cou[1][8][3)
          _cou[7][9][2) = _wmp[0]*_cou[1][9][3) & 
               +    faka * (_cou[0][9][2) - fakaac * _cou[0][9][3) )
          _cou[4][9][2) = _wmp[1]*_cou[1][9][3)
          _cou[5][9][2) = _wmp[2]*_cou[1][9][3) + fakac2 * _cou[1][3][3)
          _cou[8][4][2) = _wmp[1]*_cou[2][4][3) & 
               +  faka * (_cou[0][4][2) - fakaac * _cou[0][4][3) ) + fakac * _cou[2][1][3)
          _cou[6][4][2) = _wmp[2]*_cou[2][4][3)
          _cou[8][5][2) = _wmp[1]*_cou[2][5][3) & 
               +    faka * (_cou[0][5][2) - fakaac * _cou[0][5][3) )
          _cou[6][5][2) = _wmp[2]*_cou[2][5][3) + fakac * _cou[2][1][3)
          _cou[8][6][2) = _wmp[1]*_cou[2][6][3) & 
               +  faka * (_cou[0][6][2) - fakaac * _cou[0][6][3) ) + fakac * _cou[2][3][3)
          _cou[6][6][2) = _wmp[2]*_cou[2][6][3) + fakac * _cou[2][2][3)
          _cou[8][7][2) = _wmp[1]*_cou[2][7][3) & 
               +    faka * (_cou[0][7][2) - fakaac * _cou[0][7][3) )
          _cou[6][7][2) = _wmp[2]*_cou[2][7][3)
          _cou[8][8][2) = _wmp[1]*_cou[2][8][3) & 
               + faka * (_cou[0][8][2) - fakaac * _cou[0][8][3) ) + fakac2 * _cou[2][2][3)
          _cou[6][8][2) = _wmp[2]*_cou[2][8][3)
          _cou[8][9][2) = _wmp[1]*_cou[2][9][3) & 
               +    faka * (_cou[0][9][2) - fakaac * _cou[0][9][3) )
          _cou[6][9][2) = _wmp[2]*_cou[2][9][3) + fakac2 * _cou[2][3][3)
          _cou[9][4][2) = _wmp[2]*_cou[3][4][3) & 
               +    faka * (_cou[0][4][2) - fakaac * _cou[0][4][3) )
          _cou[9][5][2) = _wmp[2]*_cou[3][5][3) & 
               +  faka * (_cou[0][5][2) - fakaac * _cou[0][5][3) ) + fakac * _cou[3][1][3)
          _cou[9][6][2) = _wmp[2]*_cou[3][6][3) & 
               +  faka * (_cou[0][6][2) - fakaac * _cou[0][6][3) ) + fakac * _cou[3][2][3)
          _cou[9][7][2) = _wmp[2]*_cou[3][7][3) & 
               +    faka * (_cou[0][7][2) - fakaac * _cou[0][7][3) )
          _cou[9][8][2) = _wmp[2]*_cou[3][8][3) & 
               +    faka * (_cou[0][8][2) - fakaac * _cou[0][8][3) )
          _cou[9][9][2) = _wmp[2]*_cou[3][9][3) & 
               + faka * (_cou[0][9][2) - fakaac * _cou[0][9][3) ) + fakac2 * _cou[3][3][3)
          _cou[13][1][2) = _wmp[0]*_cou[4][1][3) & 
               +  faka * (_cou[2][1][2) - fakaac * _cou[2][1][3) ) + fakac * _cou[4][0][3)
          _cou[14][1][2) = _wmp[1]*_cou[4][1][3) & 
               +    faka * (_cou[1][1][2) - fakaac * _cou[1][1][3) )
          _cou[19][1][2) = _wmp[2]*_cou[4][1][3)
          _cou[13][2][2) = _wmp[0]*_cou[4][2][3) & 
               +    faka * (_cou[2][2][2) - fakaac * _cou[2][2][3) )
          _cou[14][2][2) = _wmp[1]*_cou[4][2][3) & 
               +  faka * (_cou[1][2][2) - fakaac * _cou[1][2][3) ) + fakac * _cou[4][0][3)
          _cou[19][2][2) = _wmp[2]*_cou[4][2][3)
          _cou[13][3][2) = _wmp[0]*_cou[4][3][3) & 
               +    faka * (_cou[2][3][2) - fakaac * _cou[2][3][3) )
          _cou[14][3][2) = _wmp[1]*_cou[4][3][3) & 
               +    faka * (_cou[1][3][2) - fakaac * _cou[1][3][3) )
          _cou[19][3][2) = _wmp[2]*_cou[4][3][3) + fakac * _cou[4][0][3)
          _cou[15][1][2) = _wmp[0]*_cou[5][1][3) & 
               +  faka * (_cou[3][1][2) - fakaac * _cou[3][1][3) ) + fakac * _cou[5][0][3)
          _cou[16][1][2) = _wmp[2]*_cou[5][1][3) & 
               +    faka * (_cou[1][1][2) - fakaac * _cou[1][1][3) )
          _cou[15][2][2) = _wmp[0]*_cou[5][2][3) & 
               +    faka * (_cou[3][2][2) - fakaac * _cou[3][2][3) )
          _cou[16][2][2) = _wmp[2]*_cou[5][2][3) & 
               +    faka * (_cou[1][2][2) - fakaac * _cou[1][2][3) )
          _cou[15][3][2) = _wmp[0]*_cou[5][3][3) & 
               +    faka * (_cou[3][3][2) - fakaac * _cou[3][3][3) )
          _cou[16][3][2) = _wmp[2]*_cou[5][3][3) & 
               +  faka * (_cou[1][3][2) - fakaac * _cou[1][3][3) ) + fakac * _cou[5][0][3)
          _cou[17][1][2) = _wmp[1]*_cou[6][1][3) & 
               +    faka * (_cou[3][1][2) - fakaac * _cou[3][1][3) )
          _cou[18][1][2) = _wmp[2]*_cou[6][1][3) & 
               +    faka * (_cou[2][1][2) - fakaac * _cou[2][1][3) )
          _cou[17][2][2) = _wmp[1]*_cou[6][2][3) & 
               +  faka * (_cou[3][2][2) - fakaac * _cou[3][2][3) ) + fakac * _cou[6][0][3)
          _cou[18][2][2) = _wmp[2]*_cou[6][2][3) & 
               +    faka * (_cou[2][2][2) - fakaac * _cou[2][2][3) )
          _cou[17][3][2) = _wmp[1]*_cou[6][3][3) & 
               +    faka * (_cou[3][3][2) - fakaac * _cou[3][3][3) )
          _cou[18][3][2) = _wmp[2]*_cou[6][3][3) & 
               +  faka * (_cou[2][3][2) - fakaac * _cou[2][3][3) ) + fakac * _cou[6][0][3)
          _cou[10][1][2) = _wmp[0]*_cou[7][1][3) & 
               + faka2 * (_cou[1][1][2) - fakaac * _cou[1][1][3) ) + fakac * _cou[7][0][3)
          _cou[10][2][2) = _wmp[0]*_cou[7][2][3) & 
               +    faka2 * (_cou[1][2][2) - fakaac * _cou[1][2][3) )
          _cou[10][3][2) = _wmp[0]*_cou[7][3][3) & 
               +    faka2 * (_cou[1][3][2) - fakaac * _cou[1][3][3) )
          _cou[11][1][2) = _wmp[1]*_cou[8][1][3) & 
               +    faka2 * (_cou[2][1][2) - fakaac * _cou[2][1][3) )
          _cou[11][2][2) = _wmp[1]*_cou[8][2][3) & 
               + faka2 * (_cou[2][2][2) - fakaac * _cou[2][2][3) ) + fakac * _cou[8][0][3)
          _cou[11][3][2) = _wmp[1]*_cou[8][3][3) & 
               +    faka2 * (_cou[2][3][2) - fakaac * _cou[2][3][3) )
          _cou[12][1][2) = _wmp[2]*_cou[9][1][3) & 
               +    faka2 * (_cou[3][1][2) - fakaac * _cou[3][1][3) )
          _cou[12][2][2) = _wmp[2]*_cou[9][2][3) & 
               +    faka2 * (_cou[3][2][2) - fakaac * _cou[3][2][3) )
          _cou[12][3][2) = _wmp[2]*_cou[9][3][3) & 
               + faka2 * (_cou[3][3][2) - fakaac * _cou[3][3][3) ) + fakac * _cou[9][0][3)



         
          _cou[13][4][1) = _wmp[0]*_cou[4][4][2) & 
               +  faka * (_cou[2][4][1) - fakaac * _cou[2][4][2) ) + fakac * _cou[4][2][2)
          _cou[14][4][1) = _wmp[1]*_cou[4][4][2) & 
               +  faka * (_cou[1][4][1) - fakaac * _cou[1][4][2) ) + fakac * _cou[4][1][2)
          _cou[19][4][1) = _wmp[2]*_cou[4][4][2)
          _cou[13][5][1) = _wmp[0]*_cou[4][5][2) & 
               +  faka * (_cou[2][5][1) - fakaac * _cou[2][5][2) ) + fakac * _cou[4][3][2)
          _cou[14][5][1) = _wmp[1]*_cou[4][5][2) & 
               +    faka * (_cou[1][5][1) - fakaac * _cou[1][5][2) )
          _cou[19][5][1) = _wmp[2]*_cou[4][5][2) + fakac * _cou[4][1][2)
          _cou[13][6][1) = _wmp[0]*_cou[4][6][2) & 
               +    faka * (_cou[2][6][1) - fakaac * _cou[2][6][2) )
          _cou[14][6][1) = _wmp[1]*_cou[4][6][2) & 
               +  faka * (_cou[1][6][1) - fakaac * _cou[1][6][2) ) + fakac * _cou[4][3][2)
          _cou[19][6][1) = _wmp[2]*_cou[4][6][2) & 
               +    fakac * _cou[4][2][2)
          _cou[13][7][1) = _wmp[0]*_cou[4][7][2) & 
               + faka * (_cou[2][7][1) - fakaac * _cou[2][7][2) ) + fakac2 * _cou[4][1][2)
          _cou[14][7][1) = _wmp[1]*_cou[4][7][2) & 
               +    faka * (_cou[1][7][1) - fakaac * _cou[1][7][2) )
          _cou[19][7][1) = _wmp[2]*_cou[4][7][2)
          _cou[13][8][1) = _wmp[0]*_cou[4][8][2) & 
               +    faka * (_cou[2][8][1) - fakaac * _cou[2][8][2) )
          _cou[14][8][1) = _wmp[1]*_cou[4][8][2) & 
               + faka * (_cou[1][8][1) - fakaac * _cou[1][8][2) ) + fakac2 * _cou[4][2][2)
          _cou[19][8][1) = _wmp[2]*_cou[4][8][2)
          _cou[13][9][1) = _wmp[0]*_cou[4][9][2) & 
               +    faka * (_cou[2][9][1) - fakaac * _cou[2][9][2) )
          _cou[14][9][1) = _wmp[1]*_cou[4][9][2) & 
               +    faka * (_cou[1][9][1) - fakaac * _cou[1][9][2) )
          _cou[19][9][1) = _wmp[2]*_cou[4][9][2) + fakac2 * _cou[4][3][2)
          _cou[15][4][1) = _wmp[0]*_cou[5][4][2) & 
               +  faka * (_cou[3][4][1) - fakaac * _cou[3][4][2) ) + fakac * _cou[5][2][2)
          _cou[16][4][1) = _wmp[2]*_cou[5][4][2) & 
               +    faka * (_cou[1][4][1) - fakaac * _cou[1][4][2) )
          _cou[15][5][1) = _wmp[0]*_cou[5][5][2) & 
               +  faka * (_cou[3][5][1) - fakaac * _cou[3][5][2) ) + fakac * _cou[5][3][2)
          _cou[16][5][1) = _wmp[2]*_cou[5][5][2) & 
               +  faka * (_cou[1][5][1) - fakaac * _cou[1][5][2) ) + fakac * _cou[5][1][2)
          _cou[15][6][1) = _wmp[0]*_cou[5][6][2) & 
               +    faka * (_cou[3][6][1) - fakaac * _cou[3][6][2) )
          _cou[16][6][1) = _wmp[2]*_cou[5][6][2) & 
               +  faka * (_cou[1][6][1) - fakaac * _cou[1][6][2) ) + fakac * _cou[5][2][2)
          _cou[15][7][1) = _wmp[0]*_cou[5][7][2) & 
               + faka * (_cou[3][7][1) - fakaac * _cou[3][7][2) ) + fakac2 * _cou[5][1][2)
          _cou[16][7][1) = _wmp[2]*_cou[5][7][2) & 
               +    faka * (_cou[1][7][1) - fakaac * _cou[1][7][2) )
          _cou[15][8][1) = _wmp[0]*_cou[5][8][2) & 
               +    faka * (_cou[3][8][1) - fakaac * _cou[3][8][2) )
          _cou[16][8][1) = _wmp[2]*_cou[5][8][2) & 
               +    faka * (_cou[1][8][1) - fakaac * _cou[1][8][2) )
          _cou[15][9][1) = _wmp[0]*_cou[5][9][2) & 
               +    faka * (_cou[3][9][1) - fakaac * _cou[3][9][2) )
          _cou[16][9][1) = _wmp[2]*_cou[5][9][2) & 
               + faka * (_cou[1][9][1) - fakaac * _cou[1][9][2) ) + fakac2 * _cou[5][3][2)
          _cou[17][4][1) = _wmp[1]*_cou[6][4][2) & 
               + faka * (_cou[3][4][1) - fakaac * _cou[3][4][2) ) + fakac * _cou[6][1][2)
          _cou[18][4][1) = _wmp[2]*_cou[6][4][2) & 
               +    faka * (_cou[2][4][1) - fakaac * _cou[2][4][2) )
          _cou[17][5][1) = _wmp[1]*_cou[6][5][2) & 
               +    faka * (_cou[3][5][1) - fakaac * _cou[3][5][2) )
          _cou[18][5][1) = _wmp[2]*_cou[6][5][2) & 
               +  faka * (_cou[2][5][1) - fakaac * _cou[2][5][2) ) + fakac * _cou[6][1][2)
          _cou[17][6][1) = _wmp[1]*_cou[6][6][2) & 
               +  faka * (_cou[3][6][1) - fakaac * _cou[3][6][2) ) + fakac * _cou[6][3][2)
          _cou[18][6][1) = _wmp[2]*_cou[6][6][2) & 
               +  faka * (_cou[2][6][1) - fakaac * _cou[2][6][2) ) + fakac * _cou[6][2][2)
          _cou[17][7][1) = _wmp[1]*_cou[6][7][2) & 
               +    faka * (_cou[3][7][1) - fakaac * _cou[3][7][2) )
          _cou[18][7][1) = _wmp[2]*_cou[6][7][2) & 
               +    faka * (_cou[2][7][1) - fakaac * _cou[2][7][2) )
          _cou[17][8][1) = _wmp[1]*_cou[6][8][2) & 
               + faka * (_cou[3][8][1) - fakaac * _cou[3][8][2) ) + fakac2 * _cou[6][2][2)
          _cou[18][8][1) = _wmp[2]*_cou[6][8][2) & 
               +    faka * (_cou[2][8][1) - fakaac * _cou[2][8][2) )
          _cou[17][9][1) = _wmp[1]*_cou[6][9][2) & 
               +    faka * (_cou[3][9][1) - fakaac * _cou[3][9][2) )
          _cou[18][9][1) = _wmp[2]*_cou[6][9][2) & 
               + faka * (_cou[2][9][1) - fakaac * _cou[2][9][2) ) + fakac2 * _cou[6][3][2)
          _cou[10][4][1) = _wmp[0]*_cou[7][4][2) & 
               + faka2 * (_cou[1][4][1) - fakaac * _cou[1][4][2) ) + fakac * _cou[7][2][2)
          _cou[10][5][1) = _wmp[0]*_cou[7][5][2) & 
               + faka2 * (_cou[1][5][1) - fakaac * _cou[1][5][2) ) + fakac * _cou[7][3][2)
          _cou[10][6][1) = _wmp[0]*_cou[7][6][2) & 
               +    faka2 * (_cou[1][6][1) - fakaac * _cou[1][6][2) )
          _cou[10][7][1) = _wmp[0]*_cou[7][7][2) & 
               +faka2 * (_cou[1][7][1) - fakaac * _cou[1][7][2) ) + fakac2 * _cou[7][1][2)
          _cou[10][8][1) = _wmp[0]*_cou[7][8][2) & 
               +    faka2 * (_cou[1][8][1) - fakaac * _cou[1][8][2) )
          _cou[10][9][1) = _wmp[0]*_cou[7][9][2) & 
               +    faka2 * (_cou[1][9][1) - fakaac * _cou[1][9][2) )
          _cou[11][4][1) = _wmp[1]*_cou[8][4][2) & 
               + faka2 * (_cou[2][4][1) - fakaac * _cou[2][4][2) ) + fakac * _cou[8][1][2)
          _cou[11][5][1) = _wmp[1]*_cou[8][5][2) & 
               +    faka2 * (_cou[2][5][1) - fakaac * _cou[2][5][2) )
          _cou[11][6][1) = _wmp[1]*_cou[8][6][2) & 
               + faka2 * (_cou[2][6][1) - fakaac * _cou[2][6][2) ) + fakac * _cou[8][3][2)
          _cou[11][7][1) = _wmp[1]*_cou[8][7][2) & 
               +    faka2 * (_cou[2][7][1) - fakaac * _cou[2][7][2) )
          _cou[11][8][1) = _wmp[1]*_cou[8][8][2) & 
               +faka2 * (_cou[2][8][1) - fakaac * _cou[2][8][2) ) + fakac2 * _cou[8][2][2)
          _cou[11][9][1) = _wmp[1]*_cou[8][9][2) & 
               +    faka2 * (_cou[2][9][1) - fakaac * _cou[2][9][2) )
          _cou[12][4][1) = _wmp[2]*_cou[9][4][2) & 
               +    faka2 * (_cou[3][4][1) - fakaac * _cou[3][4][2) )
          _cou[12][5][1) = _wmp[2]*_cou[9][5][2) & 
               + faka2 * (_cou[3][5][1) - fakaac * _cou[3][5][2) ) + fakac * _cou[9][1][2)
          _cou[12][6][1) = _wmp[2]*_cou[9][6][2) & 
               + faka2 * (_cou[3][6][1) - fakaac * _cou[3][6][2) ) + fakac * _cou[9][2][2)
          _cou[12][7][1) = _wmp[2]*_cou[9][7][2) & 
               +    faka2 * (_cou[3][7][1) - fakaac * _cou[3][7][2) )
          _cou[12][8][1) = _wmp[2]*_cou[9][8][2) & 
               +    faka2 * (_cou[3][8][1) - fakaac * _cou[3][8][2) )
          _cou[12][9][1) = _wmp[2]*_cou[9][9][2) & 
               +faka2 * (_cou[3][9][1) - fakaac * _cou[3][9][2) ) + fakac2 * _cou[9][3][2)

     ndif


       if(abs(lmax1)+abs(lmax2).ge.4) then
          _cou[20][0][0) = _wmp[0]*_cou[10][0][1) & 
               +    faka3 * (_cou[7][0][0) - fakaac * _cou[7][0][1) ) 
          _cou[23][0][0) = _wmp[1]*_cou[10][0][1) 
          _cou[25][0][0) = _wmp[2]*_cou[10][0][1) 
          _cou[24][0][0) = _wmp[0]*_cou[11][0][1) 
          _cou[21][0][0) = _wmp[1]*_cou[11][0][1) & 
               +    faka3 * (_cou[8][0][0) - fakaac * _cou[8][0][1) ) 
          _cou[27][0][0) = _wmp[2]*_cou[11][0][1) 
          _cou[26][0][0) = _wmp[0]*_cou[12][0][1) 
          _cou[28][0][0) = _wmp[1]*_cou[12][0][1) 
          _cou[22][0][0) = _wmp[2]*_cou[12][0][1) & 
               +    faka3 * (_cou[9][0][0) - fakaac * _cou[9][0][1) ) 
          _cou[31][0][0) = _wmp[1]*_cou[13][0][1) & 
               +    faka * (_cou[7][0][0) - fakaac * _cou[7][0][1) ) 
          _cou[32][0][0) = _wmp[2]*_cou[13][0][1) 
          _cou[33][0][0) = _wmp[2]*_cou[14][0][1) 
          _cou[29][0][0) = _wmp[2]*_cou[15][0][1) & 
               +    faka * (_cou[7][0][0) - fakaac * _cou[7][0][1) ) 
          _cou[34][0][0) = _wmp[1]*_cou[16][0][1) 
          _cou[30][0][0) = _wmp[2]*_cou[17][0][1) & 
               +    faka * (_cou[8][0][0) - fakaac * _cou[8][0][1) ) 

          _cou[0][20][0) = _wmq[0]*_cou[0][10][1) & 
               +    fakc3 * (_cou[0][7][0) - fakaca * _cou[0][7][1) ) 
          _cou[0][23][0) = _wmq[1]*_cou[0][10][1) 
          _cou[0][25][0) = _wmq[2]*_cou[0][10][1) 
          _cou[0][24][0) = _wmq[0]*_cou[0][11][1) 
          _cou[0][21][0) = _wmq[1]*_cou[0][11][1) & 
               +    fakc3 * (_cou[0][8][0) - fakaca * _cou[0][8][1) ) 
          _cou[0][27][0) = _wmq[2]*_cou[0][11][1) 
          _cou[0][26][0) = _wmq[0]*_cou[0][12][1) 
          _cou[0][28][0) = _wmq[1]*_cou[0][12][1) 
          _cou[0][22][0) = _wmq[2]*_cou[0][12][1) & 
               +    fakc3 * (_cou[0][9][0) - fakaca * _cou[0][9][1) ) 
          _cou[0][31][0) = _wmq[1]*_cou[0][13][1) & 
               +    fakc * (_cou[0][7][0) - fakaca * _cou[0][7][1) ) 
          _cou[0][32][0) = _wmq[2]*_cou[0][13][1) 
          _cou[0][33][0) = _wmq[2]*_cou[0][14][1) 
          _cou[0][29][0) = _wmq[2]*_cou[0][15][1) & 
               +    fakc * (_cou[0][7][0) - fakaca * _cou[0][7][1) ) 
          _cou[0][34][0) = _wmq[1]*_cou[0][16][1) 
          _cou[0][30][0) = _wmq[2]*_cou[0][17][1) & 
               +    fakc * (_cou[0][8][0) - fakaca * _cou[0][8][1) ) 
       endif


       if(abs(lmax1)+abs(lmax2).ge.5) then

          _cou[20][0][1) = _wmp[0]*_cou[10][0][2) &
               +    faka3 * (_cou[7][0][1) - fakaac * _cou[7][0][2) ) 
          _cou[23][0][1) = _wmp[1]*_cou[10][0][2) 
          _cou[25][0][1) = _wmp[2]*_cou[10][0][2) 
          _cou[24][0][1) = _wmp[0]*_cou[11][0][2) 
          _cou[21][0][1) = _wmp[1]*_cou[11][0][2) & 
               +    faka3 * (_cou[8][0][1) - fakaac * _cou[8][0][2) ) 
          _cou[27][0][1) = _wmp[2]*_cou[11][0][2) 
          _cou[26][0][1) = _wmp[0]*_cou[12][0][2) 
          _cou[28][0][1) = _wmp[1]*_cou[12][0][2) 
          _cou[22][0][1) = _wmp[2]*_cou[12][0][2) &
               +    faka3 * (_cou[9][0][1) - fakaac * _cou[9][0][2) ) 
          _cou[31][0][1) = _wmp[1]*_cou[13][0][2)  & 
               +    faka * (_cou[7][0][1) - fakaac * _cou[7][0][2) ) 
          _cou[32][0][1) = _wmp[2]*_cou[13][0][2) 
          _cou[33][0][1) = _wmp[2]*_cou[14][0][2) 
          _cou[29][0][1) = _wmp[2]*_cou[15][0][2) & 
               +    faka * (_cou[7][0][1) - fakaac * _cou[7][0][2) ) 
          _cou[34][0][1) = _wmp[1]*_cou[16][0][2) 
          _cou[30][0][1) = _wmp[2]*_cou[17][0][2)  & 
               +    faka * (_cou[8][0][1) - fakaac * _cou[8][0][2) ) 

          _cou[0][20][1) = _wmq[0]*_cou[0][10][2)  & 
               +    fakc3 * (_cou[0][7][1) - fakaca * _cou[0][7][2) ) 
          _cou[0][23][1) = _wmq[1]*_cou[0][10][2) 
          _cou[0][25][1) = _wmq[2]*_cou[0][10][2) 
          _cou[0][24][1) = _wmq[0]*_cou[0][11][2) 
          _cou[0][21][1) = _wmq[1]*_cou[0][11][2)  & 
               +    fakc3 * (_cou[0][8][1) - fakaca * _cou[0][8][2) ) 
          _cou[0][27][1) = _wmq[2]*_cou[0][11][2) 
          _cou[0][26][1) = _wmq[0]*_cou[0][12][2) 
          _cou[0][28][1) = _wmq[1]*_cou[0][12][2) 
          _cou[0][22][1) = _wmq[2]*_cou[0][12][2)  & 
               +    fakc3 * (_cou[0][9][1) - fakaca * _cou[0][9][2) ) 
          _cou[0][31][1) = _wmq[1]*_cou[0][13][2)  & 
               +    fakc * (_cou[0][7][1) - fakaca * _cou[0][7][2) ) 
          _cou[0][32][1) = _wmq[2]*_cou[0][13][2) 
          _cou[0][33][1) = _wmq[2]*_cou[0][14][2) 
          _cou[0][29][1) = _wmq[2]*_cou[0][15][2)  & 
               +    fakc * (_cou[0][7][1) - fakaca * _cou[0][7][2) ) 
          _cou[0][34][1) = _wmq[1]*_cou[0][16][2)
          _cou[0][30][1) = _wmq[2]*_cou[0][17][2)  & 
               +    fakc * (_cou[0][8][1) - fakaca * _cou[0][8][2) ) 

          _cou[1][20][0) = _wmp[0]*_cou[0][20][1) +  fakac4 * _cou[0][10][1)
          _cou[2][20][0) = _wmp[1]*_cou[0][20][1) 
          _cou[3][20][0) = _wmp[2]*_cou[0][20][1) 
          _cou[1][21][0) = _wmp[0]*_cou[0][21][1) 
          _cou[2][21][0) = _wmp[1]*_cou[0][21][1) +  fakac4 * _cou[0][11][1)
          _cou[3][21][0) = _wmp[2]*_cou[0][21][1) 
          _cou[1][22][0) = _wmp[0]*_cou[0][22][1) 
          _cou[2][22][0) = _wmp[1]*_cou[0][22][1) 
          _cou[3][22][0) = _wmp[2]*_cou[0][22][1) +  fakac4 * _cou[0][12][1)
          _cou[1][23][0) = _wmp[0]*_cou[0][23][1) +  fakac3 * _cou[0][13][1)
          _cou[2][23][0) = _wmp[1]*_cou[0][23][1) +  fakac * _cou[0][10][1)
          _cou[3][23][0) = _wmp[2]*_cou[0][23][1) 
          _cou[1][24][0) = _wmp[0]*_cou[0][24][1) +  fakac * _cou[0][11][1)
          _cou[2][24][0) = _wmp[1]*_cou[0][24][1) +  fakac3 * _cou[0][14][1)
          _cou[3][24][0) = _wmp[2]*_cou[0][24][1) 
          _cou[1][25][0) = _wmp[0]*_cou[0][25][1) +  fakac3 * _cou[0][15][1)
          _cou[2][25][0) = _wmp[1]*_cou[0][25][1) 
          _cou[3][25][0) = _wmp[2]*_cou[0][25][1) +  fakac * _cou[0][10][1)
          _cou[1][26][0) = _wmp[0]*_cou[0][26][1) +  fakac * _cou[0][12][1)
          _cou[2][26][0) = _wmp[1]*_cou[0][26][1) 
          _cou[3][26][0) = _wmp[2]*_cou[0][26][1) +  fakac3 * _cou[0][16][1)
          _cou[1][27][0) = _wmp[0]*_cou[0][27][1) 
          _cou[2][27][0) = _wmp[1]*_cou[0][27][1) +  fakac3 * _cou[0][17][1)
          _cou[3][27][0) = _wmp[2]*_cou[0][27][1) +  fakac * _cou[0][11][1)
          _cou[1][28][0) = _wmp[0]*_cou[0][28][1) 
          _cou[2][28][0) = _wmp[1]*_cou[0][28][1) +  fakac * _cou[0][12][1)
          _cou[3][28][0) = _wmp[2]*_cou[0][28][1) +  fakac3 * _cou[0][18][1)
          _cou[1][29][0) = _wmp[0]*_cou[0][29][1) +  fakac2 * _cou[0][16][1)
          _cou[2][29][0) = _wmp[1]*_cou[0][29][1) 
          _cou[3][29][0) = _wmp[2]*_cou[0][29][1) +  fakac2 * _cou[0][15][1)
          _cou[1][30][0) = _wmp[0]*_cou[0][30][1) 
          _cou[2][30][0) = _wmp[1]*_cou[0][30][1) +  fakac2 * _cou[0][18][1)
          _cou[3][30][0) = _wmp[2]*_cou[0][30][1) +  fakac2 * _cou[0][17][1)
          _cou[1][31][0) = _wmp[0]*_cou[0][31][1) +  fakac2 * _cou[0][14][1)
          _cou[2][31][0) = _wmp[1]*_cou[0][31][1) +  fakac2 * _cou[0][13][1)
          _cou[3][31][0) = _wmp[2]*_cou[0][31][1) 
          _cou[1][32][0) = _wmp[0]*_cou[0][32][1) +  fakac2 * _cou[0][19][1)
          _cou[2][32][0) = _wmp[1]*_cou[0][32][1) +  fakac * _cou[0][15][1)
          _cou[3][32][0) = _wmp[2]*_cou[0][32][1) +  fakac * _cou[0][13][1)
          _cou[1][33][0) = _wmp[0]*_cou[0][33][1) +  fakac * _cou[0][17][1)
          _cou[2][33][0) = _wmp[1]*_cou[0][33][1) +  fakac2 * _cou[0][19][1)
          _cou[3][33][0) = _wmp[2]*_cou[0][33][1) +  fakac * _cou[0][14][1)
          _cou[1][34][0) = _wmp[0]*_cou[0][34][1) +  fakac * _cou[0][18][1)
          _cou[2][34][0) = _wmp[1]*_cou[0][34][1) +  fakac * _cou[0][16][1)
          _cou[3][34][0) = _wmp[2]*_cou[0][34][1) +  fakac2 * _cou[0][19][1)


          _cou[20][1][0) = _wmp[0]*_cou[10][1][1)  & 
               +  faka3 * (_cou[7][1][0) - fakaac * _cou[7][1][1) ) + fakac * _cou[10][0][1)
          _cou[23][1][0) = _wmp[1]*_cou[10][1][1) 
          _cou[25][1][0) = _wmp[2]*_cou[10][1][1) 
          _cou[20][2][0) = _wmp[0]*_cou[10][2][1)  & 
               +    faka3 * (_cou[7][2][0) - fakaac * _cou[7][2][1) ) 
          _cou[23][2][0) = _wmp[1]*_cou[10][2][1)  & 
               +    fakac * _cou[10][0][1)
          _cou[25][2][0) = _wmp[2]*_cou[10][2][1) 
          _cou[20][3][0) = _wmp[0]*_cou[10][3][1)  & 
               +    faka3 * (_cou[7][3][0) - fakaac * _cou[7][3][1) ) 
          _cou[23][3][0) = _wmp[1]*_cou[10][3][1) 
          _cou[25][3][0) = _wmp[2]*_cou[10][3][1) +  fakac * _cou[10][0][1)
          _cou[24][1][0) = _wmp[0]*_cou[11][1][1) +  fakac * _cou[11][0][1)
          _cou[21][1][0) = _wmp[1]*_cou[11][1][1)  & 
               +    faka3 * (_cou[8][1][0) - fakaac * _cou[8][1][1) ) 
          _cou[27][1][0) = _wmp[2]*_cou[11][1][1) 
          _cou[24][2][0) = _wmp[0]*_cou[11][2][1) 
          _cou[21][2][0) = _wmp[1]*_cou[11][2][1)  & 
               +  faka3 * (_cou[8][2][0) - fakaac * _cou[8][2][1) ) + fakac * _cou[11][0][1)
          _cou[27][2][0) = _wmp[2]*_cou[11][2][1) 
          _cou[24][3][0) = _wmp[0]*_cou[11][3][1) 
          _cou[21][3][0) = _wmp[1]*_cou[11][3][1)  & 
               +    faka3 * (_cou[8][3][0) - fakaac * _cou[8][3][1) ) 
          _cou[27][3][0) = _wmp[2]*_cou[11][3][1) +  fakac * _cou[11][0][1)
          _cou[26][1][0) = _wmp[0]*_cou[12][1][1) +  fakac * _cou[12][0][1)
          _cou[28][1][0) = _wmp[1]*_cou[12][1][1) 
          _cou[22][1][0) = _wmp[2]*_cou[12][1][1)  & 
               +    faka3 * (_cou[9][1][0) - fakaac * _cou[9][1][1) ) 
          _cou[26][2][0) = _wmp[0]*_cou[12][2][1) 
          _cou[28][2][0) = _wmp[1]*_cou[12][2][1) +  fakac * _cou[12][0][1)
          _cou[22][2][0) = _wmp[2]*_cou[12][2][1)  & 
               +    faka3 * (_cou[9][2][0) - fakaac * _cou[9][2][1) ) 
          _cou[26][3][0) = _wmp[0]*_cou[12][3][1) 
          _cou[28][3][0) = _wmp[1]*_cou[12][3][1) 
          _cou[22][3][0) = _wmp[2]*_cou[12][3][1)  & 
               +  faka3 * (_cou[9][3][0) - fakaac * _cou[9][3][1) ) + fakac * _cou[12][0][1)
          _cou[31][1][0) = _wmp[1]*_cou[13][1][1)  & 
               +    faka * (_cou[7][1][0) - fakaac * _cou[7][1][1) ) 
          _cou[32][1][0) = _wmp[2]*_cou[13][1][1) 
          _cou[31][2][0) = _wmp[1]*_cou[13][2][1)  & 
               +   faka * (_cou[7][2][0) - fakaac * _cou[7][2][1) ) + fakac * _cou[13][0][1)
          _cou[32][2][0) = _wmp[2]*_cou[13][2][1) 
          _cou[31][3][0) = _wmp[1]*_cou[13][3][1)  & 
               +    faka * (_cou[7][3][0) - fakaac * _cou[7][3][1) ) 
          _cou[32][3][0) = _wmp[2]*_cou[13][3][1) +  fakac * _cou[13][0][1)
          _cou[33][1][0) = _wmp[2]*_cou[14][1][1) 
          _cou[33][2][0) = _wmp[2]*_cou[14][2][1) 
          _cou[33][3][0) = _wmp[2]*_cou[14][3][1) +  fakac * _cou[14][0][1)
          _cou[29][1][0) = _wmp[2]*_cou[15][1][1)  & 
               +    faka * (_cou[7][1][0) - fakaac * _cou[7][1][1) ) 
          _cou[29][2][0) = _wmp[2]*_cou[15][2][1)  & 
               +    faka * (_cou[7][2][0) - fakaac * _cou[7][2][1) ) 
          _cou[29][3][0) = _wmp[2]*_cou[15][3][1)  & 
               +   faka * (_cou[7][3][0) - fakaac * _cou[7][3][1) ) + fakac * _cou[15][0][1)
          _cou[34][1][0) = _wmp[1]*_cou[16][1][1) 
          _cou[34][2][0) = _wmp[1]*_cou[16][2][1) +  fakac * _cou[16][0][1)
          _cou[34][3][0) = _wmp[1]*_cou[16][3][1) 
          _cou[30][1][0) = _wmp[2]*_cou[17][1][1)  & 
               +    faka * (_cou[8][1][0) - fakaac * _cou[8][1][1) ) 
          _cou[30][2][0) = _wmp[2]*_cou[17][2][1)  & 
               +    faka * (_cou[8][2][0) - fakaac * _cou[8][2][1) ) 
          _cou[30][3][0) = _wmp[2]*_cou[17][3][1)  & 
               +   faka * (_cou[8][3][0) - fakaac * _cou[8][3][1) ) + fakac * _cou[17][0][1)

       endif


       if(abs(lmax1)+abs(lmax2).ge.6) then
          _cou[20][0][2) = _wmp[0]*_cou[10][0][3)  & 
               +    faka3 * (_cou[7][0][2) - fakaac * _cou[7][0][3) ) 
          _cou[23][0][2) = _wmp[1]*_cou[10][0][3) 
          _cou[25][0][2) = _wmp[2]*_cou[10][0][3) 
          _cou[24][0][2) = _wmp[0]*_cou[11][0][3) 
          _cou[21][0][2) = _wmp[1]*_cou[11][0][3)  & 
               +    faka3 * (_cou[8][0][2) - fakaac * _cou[8][0][3) ) 
          _cou[27][0][2) = _wmp[2]*_cou[11][0][3) 
          _cou[26][0][2) = _wmp[0]*_cou[12][0][3) 
          _cou[28][0][2) = _wmp[1]*_cou[12][0][3) 
          _cou[22][0][2) = _wmp[2]*_cou[12][0][3)  & 
               +    faka3 * (_cou[9][0][2) - fakaac * _cou[9][0][3) ) 
          _cou[31][0][2) = _wmp[1]*_cou[13][0][3)  & 
               +    faka * (_cou[7][0][2) - fakaac * _cou[7][0][3) ) 
          _cou[32][0][2) = _wmp[2]*_cou[13][0][3) 
          _cou[33][0][2) = _wmp[2]*_cou[14][0][3) 
          _cou[29][0][2) = _wmp[2]*_cou[15][0][3)  & 
               +    faka * (_cou[7][0][2) - fakaac * _cou[7][0][3) ) 
          _cou[34][0][2) = _wmp[1]*_cou[16][0][3) 
          _cou[30][0][2) = _wmp[2]*_cou[17][0][3)  & 
               +    faka * (_cou[8][0][2) - fakaac * _cou[8][0][3) ) 

          _cou[0][20][2) = _wmq[0]*_cou[0][10][3)  & 
               +    fakc3 * (_cou[0][7][2) - fakaca * _cou[0][7][3) ) 
          _cou[0][23][2) = _wmq[1]*_cou[0][10][3) 
          _cou[0][25][2) = _wmq[2]*_cou[0][10][3) 
          _cou[0][24][2) = _wmq[0]*_cou[0][11][3) 
          _cou[0][21][2) = _wmq[1]*_cou[0][11][3)  & 
               +    fakc3 * (_cou[0][8][2) - fakaca * _cou[0][8][3) ) 
          _cou[0][27][2) = _wmq[2]*_cou[0][11][3) 
          _cou[0][26][2) = _wmq[0]*_cou[0][12][3) 
          _cou[0][28][2) = _wmq[1]*_cou[0][12][3) 
          _cou[0][22][2) = _wmq[2]*_cou[0][12][3)  & 
               +    fakc3 * (_cou[0][9][2) - fakaca * _cou[0][9][3) ) 
          _cou[0][31][2) = _wmq[1]*_cou[0][13][3)  & 
               +    fakc  * (_cou[0][7][2) - fakaca * _cou[0][7][3) ) 
          _cou[0][32][2) = _wmq[2]*_cou[0][13][3) 
          _cou[0][33][2) = _wmq[2]*_cou[0][14][3) 
          _cou[0][29][2) = _wmq[2]*_cou[0][15][3)  & 
               +    fakc  * (_cou[0][7][2) - fakaca * _cou[0][7][3) ) 
          _cou[0][34][2) = _wmq[1]*_cou[0][16][3) 
          _cou[0][30][2) = _wmq[2]*_cou[0][17][3)  & 
               +    fakc  * (_cou[0][8][2) - fakaca * _cou[0][8][3) ) 

          _cou[1][20][1) = _wmp[0]*_cou[0][20][2) +  fakac4 * _cou[0][10][2)
          _cou[2][20][1) = _wmp[1]*_cou[0][20][2) 
          _cou[3][20][1) = _wmp[2]*_cou[0][20][2) 
          _cou[1][21][1) = _wmp[0]*_cou[0][21][2) 
          _cou[2][21][1) = _wmp[1]*_cou[0][21][2) +  fakac4 * _cou[0][11][2)
          _cou[3][21][1) = _wmp[2]*_cou[0][21][2) 
          _cou[1][22][1) = _wmp[0]*_cou[0][22][2) 
          _cou[2][22][1) = _wmp[1]*_cou[0][22][2) 
          _cou[3][22][1) = _wmp[2]*_cou[0][22][2) +  fakac4 * _cou[0][12][2)
          _cou[1][23][1) = _wmp[0]*_cou[0][23][2) +  fakac3 * _cou[0][13][2)
          _cou[2][23][1) = _wmp[1]*_cou[0][23][2) +  fakac * _cou[0][10][2)
          _cou[3][23][1) = _wmp[2]*_cou[0][23][2) 
          _cou[1][24][1) = _wmp[0]*_cou[0][24][2) +  fakac * _cou[0][11][2)
          _cou[2][24][1) = _wmp[1]*_cou[0][24][2) +  fakac3 * _cou[0][14][2)
          _cou[3][24][1) = _wmp[2]*_cou[0][24][2) 
          _cou[1][25][1) = _wmp[0]*_cou[0][25][2) +  fakac3 * _cou[0][15][2)
          _cou[2][25][1) = _wmp[1]*_cou[0][25][2) 
          _cou[3][25][1) = _wmp[2]*_cou[0][25][2) +  fakac * _cou[0][10][2)
          _cou[1][26][1) = _wmp[0]*_cou[0][26][2) +  fakac * _cou[0][12][2)
          _cou[2][26][1) = _wmp[1]*_cou[0][26][2) 
          _cou[3][26][1) = _wmp[2]*_cou[0][26][2) +  fakac3 * _cou[0][16][2)
          _cou[1][27][1) = _wmp[0]*_cou[0][27][2) 
          _cou[2][27][1) = _wmp[1]*_cou[0][27][2) +  fakac3 * _cou[0][17][2)
          _cou[3][27][1) = _wmp[2]*_cou[0][27][2) +  fakac * _cou[0][11][2)
          _cou[1][28][1) = _wmp[0]*_cou[0][28][2) 
          _cou[2][28][1) = _wmp[1]*_cou[0][28][2) +  fakac * _cou[0][12][2)
          _cou[3][28][1) = _wmp[2]*_cou[0][28][2) +  fakac3 * _cou[0][18][2)
          _cou[1][29][1) = _wmp[0]*_cou[0][29][2) +  fakac2 * _cou[0][16][2)
          _cou[2][29][1) = _wmp[1]*_cou[0][29][2) 
          _cou[3][29][1) = _wmp[2]*_cou[0][29][2) +  fakac2 * _cou[0][15][2)
          _cou[1][30][1) = _wmp[0]*_cou[0][30][2) 
          _cou[2][30][1) = _wmp[1]*_cou[0][30][2) +  fakac2 * _cou[0][18][2)
          _cou[3][30][1) = _wmp[2]*_cou[0][30][2) +  fakac2 * _cou[0][17][2)
          _cou[1][31][1) = _wmp[0]*_cou[0][31][2) +  fakac2 * _cou[0][14][2)
          _cou[2][31][1) = _wmp[1]*_cou[0][31][2) +  fakac2 * _cou[0][13][2)
          _cou[3][31][1) = _wmp[2]*_cou[0][31][2) 
          _cou[1][32][1) = _wmp[0]*_cou[0][32][2) +  fakac2 * _cou[0][19][2)
          _cou[2][32][1) = _wmp[1]*_cou[0][32][2) +  fakac * _cou[0][15][2)
          _cou[3][32][1) = _wmp[2]*_cou[0][32][2) +  fakac * _cou[0][13][2)
          _cou[1][33][1) = _wmp[0]*_cou[0][33][2) +  fakac * _cou[0][17][2)
          _cou[2][33][1) = _wmp[1]*_cou[0][33][2) +  fakac2 * _cou[0][19][2)
          _cou[3][33][1) = _wmp[2]*_cou[0][33][2) +  fakac * _cou[0][14][2)
          _cou[1][34][1) = _wmp[0]*_cou[0][34][2) +  fakac * _cou[0][18][2)
          _cou[2][34][1) = _wmp[1]*_cou[0][34][2) +  fakac * _cou[0][16][2)
          _cou[3][34][1) = _wmp[2]*_cou[0][34][2) +  fakac2 * _cou[0][19][2)


          _cou[20][1][1) = _wmp[0]*_cou[10][1][2)  & 
               +  faka3 * (_cou[7][1][1) - fakaac * _cou[7][1][2) ) + fakac * _cou[10][0][2)
          _cou[23][1][1) = _wmp[1]*_cou[10][1][2) 
          _cou[25][1][1) = _wmp[2]*_cou[10][1][2) 
          _cou[20][2][1) = _wmp[0]*_cou[10][2][2)  & 
               +    faka3 * (_cou[7][2][1) - fakaac * _cou[7][2][2) ) 
          _cou[23][2][1) = _wmp[1]*_cou[10][2][2) +  fakac * _cou[10][0][2)
          _cou[25][2][1) = _wmp[2]*_cou[10][2][2) 
          _cou[20][3][1) = _wmp[0]*_cou[10][3][2)  & 
               +    faka3 * (_cou[7][3][1) - fakaac * _cou[7][3][2) ) 
          _cou[23][3][1) = _wmp[1]*_cou[10][3][2) 
          _cou[25][3][1) = _wmp[2]*_cou[10][3][2) +  fakac * _cou[10][0][2)
          _cou[24][1][1) = _wmp[0]*_cou[11][1][2) +  fakac * _cou[11][0][2)
          _cou[21][1][1) = _wmp[1]*_cou[11][1][2)  & 
               +    faka3 * (_cou[8][1][1) - fakaac * _cou[8][1][2) ) 
          _cou[27][1][1) = _wmp[2]*_cou[11][1][2) 
          _cou[24][2][1) = _wmp[0]*_cou[11][2][2) 
          _cou[21][2][1) = _wmp[1]*_cou[11][2][2)  & 
               +  faka3 * (_cou[8][2][1) - fakaac * _cou[8][2][2) ) + fakac * _cou[11][0][2)
          _cou[27][2][1) = _wmp[2]*_cou[11][2][2) 
          _cou[24][3][1) = _wmp[0]*_cou[11][3][2) 
          _cou[21][3][1) = _wmp[1]*_cou[11][3][2)  & 
               +    faka3 * (_cou[8][3][1) - fakaac * _cou[8][3][2) ) 
          _cou[27][3][1) = _wmp[2]*_cou[11][3][2) +  fakac * _cou[11][0][2)
          _cou[26][1][1) = _wmp[0]*_cou[12][1][2) +  fakac * _cou[12][0][2)
          _cou[28][1][1) = _wmp[1]*_cou[12][1][2) 
          _cou[22][1][1) = _wmp[2]*_cou[12][1][2)  & 
               +    faka3 * (_cou[9][1][1) - fakaac * _cou[9][1][2) ) 
          _cou[26][2][1) = _wmp[0]*_cou[12][2][2) 
          _cou[28][2][1) = _wmp[1]*_cou[12][2][2) +  fakac * _cou[12][0][2)
          _cou[22][2][1) = _wmp[2]*_cou[12][2][2)  & 
               +    faka3 * (_cou[9][2][1) - fakaac * _cou[9][2][2) ) 
          _cou[26][3][1) = _wmp[0]*_cou[12][3][2) 
          _cou[28][3][1) = _wmp[1]*_cou[12][3][2) 
          _cou[22][3][1) = _wmp[2]*_cou[12][3][2)  & 
               +  faka3 * (_cou[9][3][1) - fakaac * _cou[9][3][2) ) + fakac * _cou[12][0][2)
          _cou[31][1][1) = _wmp[1]*_cou[13][1][2)  & 
               +    faka * (_cou[7][1][1) - fakaac * _cou[7][1][2) ) 
          _cou[32][1][1) = _wmp[2]*_cou[13][1][2) 
          _cou[31][2][1) = _wmp[1]*_cou[13][2][2)  & 
               +   faka * (_cou[7][2][1) - fakaac * _cou[7][2][2) ) + fakac * _cou[13][0][2)
          _cou[32][2][1) = _wmp[2]*_cou[13][2][2) 
          _cou[31][3][1) = _wmp[1]*_cou[13][3][2)  & 
               +    faka * (_cou[7][3][1) - fakaac * _cou[7][3][2) ) 
          _cou[32][3][1) = _wmp[2]*_cou[13][3][2) +  fakac * _cou[13][0][2)
          _cou[33][1][1) = _wmp[2]*_cou[14][1][2) 
          _cou[33][2][1) = _wmp[2]*_cou[14][2][2) 
          _cou[33][3][1) = _wmp[2]*_cou[14][3][2) +  fakac * _cou[14][0][2)
          _cou[29][1][1) = _wmp[2]*_cou[15][1][2)  & 
               +    faka * (_cou[7][1][1) - fakaac * _cou[7][1][2) ) 
          _cou[29][2][1) = _wmp[2]*_cou[15][2][2)  & 
               +    faka * (_cou[7][2][1) - fakaac * _cou[7][2][2) ) 
          _cou[29][3][1) = _wmp[2]*_cou[15][3][2)  & 
               +   faka * (_cou[7][3][1) - fakaac * _cou[7][3][2) ) + fakac * _cou[15][0][2)
          _cou[34][1][1) = _wmp[1]*_cou[16][1][2) 
          _cou[34][2][1) = _wmp[1]*_cou[16][2][2) +  fakac * _cou[16][0][2)
          _cou[34][3][1) = _wmp[1]*_cou[16][3][2) 
          _cou[30][1][1) = _wmp[2]*_cou[17][1][2)  & 
               +    faka * (_cou[8][1][1) - fakaac * _cou[8][1][2) ) 
          _cou[30][2][1) = _wmp[2]*_cou[17][2][2)  & 
               +    faka * (_cou[8][2][1) - fakaac * _cou[8][2][2) ) 
          _cou[30][3][1) = _wmp[2]*_cou[17][3][2)  & 
               +   faka * (_cou[8][3][1) - fakaac * _cou[8][3][2) ) + fakac * _cou[17][0][2)


          _cou[7][20][0) = _wmp[0]*_cou[1][20][1)  & 
               +  faka * (_cou[0][20][0) - fakaac * _cou[0][20][1) ) + fakac4 * _cou[1][10][1)
          _cou[4][20][0) = _wmp[1]*_cou[1][20][1) 
          _cou[5][20][0) = _wmp[2]*_cou[1][20][1) 
          _cou[7][21][0) = _wmp[0]*_cou[1][21][1)  & 
               +    faka * (_cou[0][21][0) - fakaac * _cou[0][21][1) ) 
          _cou[4][21][0) = _wmp[1]*_cou[1][21][1) +  fakac4 * _cou[1][11][1)
          _cou[5][21][0) = _wmp[2]*_cou[1][21][1) 
          _cou[7][22][0) = _wmp[0]*_cou[1][22][1)  & 
               +    faka * (_cou[0][22][0) - fakaac * _cou[0][22][1) ) 
          _cou[4][22][0) = _wmp[1]*_cou[1][22][1) 
          _cou[5][22][0) = _wmp[2]*_cou[1][22][1) +  fakac4 * _cou[1][12][1)
          _cou[7][23][0) = _wmp[0]*_cou[1][23][1)  & 
               +  faka * (_cou[0][23][0) - fakaac * _cou[0][23][1) ) + fakac3 * _cou[1][13][1)
          _cou[4][23][0) = _wmp[1]*_cou[1][23][1) +  fakac * _cou[1][10][1)
          _cou[5][23][0) = _wmp[2]*_cou[1][23][1) 
          _cou[7][24][0) = _wmp[0]*_cou[1][24][1)  & 
               +   faka * (_cou[0][24][0) - fakaac * _cou[0][24][1) ) + fakac * _cou[1][11][1)
          _cou[4][24][0) = _wmp[1]*_cou[1][24][1) +  fakac3 * _cou[1][14][1)
          _cou[5][24][0) = _wmp[2]*_cou[1][24][1) 
          _cou[7][25][0) = _wmp[0]*_cou[1][25][1)  & 
               +  faka * (_cou[0][25][0) - fakaac * _cou[0][25][1) ) + fakac3 * _cou[1][15][1)
          _cou[4][25][0) = _wmp[1]*_cou[1][25][1) 
          _cou[5][25][0) = _wmp[2]*_cou[1][25][1) +  fakac * _cou[1][10][1)
          _cou[7][26][0) = _wmp[0]*_cou[1][26][1)  & 
               +   faka * (_cou[0][26][0) - fakaac * _cou[0][26][1) ) + fakac * _cou[1][12][1)
          _cou[4][26][0) = _wmp[1]*_cou[1][26][1) 
          _cou[5][26][0) = _wmp[2]*_cou[1][26][1) +  fakac3 * _cou[1][16][1)
          _cou[7][27][0) = _wmp[0]*_cou[1][27][1)  & 
               +    faka * (_cou[0][27][0) - fakaac * _cou[0][27][1) ) 
          _cou[4][27][0) = _wmp[1]*_cou[1][27][1) +  fakac3 * _cou[1][17][1)
          _cou[5][27][0) = _wmp[2]*_cou[1][27][1) +  fakac * _cou[1][11][1)
          _cou[7][28][0) = _wmp[0]*_cou[1][28][1)  & 
               +    faka * (_cou[0][28][0) - fakaac * _cou[0][28][1) ) 
          _cou[4][28][0) = _wmp[1]*_cou[1][28][1) +  fakac * _cou[1][12][1)
          _cou[5][28][0) = _wmp[2]*_cou[1][28][1) +  fakac3 * _cou[1][18][1)
          _cou[7][29][0) = _wmp[0]*_cou[1][29][1)  & 
               +  faka * (_cou[0][29][0) - fakaac * _cou[0][29][1) ) + fakac2 * _cou[1][16][1)
          _cou[4][29][0) = _wmp[1]*_cou[1][29][1) 
          _cou[5][29][0) = _wmp[2]*_cou[1][29][1) +  fakac2 * _cou[1][15][1)
          _cou[7][30][0) = _wmp[0]*_cou[1][30][1)  & 
               +    faka * (_cou[0][30][0) - fakaac * _cou[0][30][1) ) 
          _cou[4][30][0) = _wmp[1]*_cou[1][30][1) +  fakac2 * _cou[1][18][1)
          _cou[5][30][0) = _wmp[2]*_cou[1][30][1) +  fakac2 * _cou[1][17][1)
          _cou[7][31][0) = _wmp[0]*_cou[1][31][1)  & 
               +  faka * (_cou[0][31][0) - fakaac * _cou[0][31][1) ) + fakac2 * _cou[1][14][1)
          _cou[4][31][0) = _wmp[1]*_cou[1][31][1) +  fakac2 * _cou[1][13][1)
          _cou[5][31][0) = _wmp[2]*_cou[1][31][1) 
          _cou[7][32][0) = _wmp[0]*_cou[1][32][1)  & 
               +  faka * (_cou[0][32][0) - fakaac * _cou[0][32][1) ) + fakac2 * _cou[1][19][1)
          _cou[4][32][0) = _wmp[1]*_cou[1][32][1) +  fakac * _cou[1][15][1)
          _cou[5][32][0) = _wmp[2]*_cou[1][32][1) +  fakac * _cou[1][13][1)
          _cou[7][33][0) = _wmp[0]*_cou[1][33][1)  & 
               +  faka * (_cou[0][33][0) - fakaac * _cou[0][33][1) ) + fakac * _cou[1][17][1)
          _cou[4][33][0) = _wmp[1]*_cou[1][33][1) +  fakac2 * _cou[1][19][1)
          _cou[5][33][0) = _wmp[2]*_cou[1][33][1) +  fakac * _cou[1][14][1)
          _cou[7][34][0) = _wmp[0]*_cou[1][34][1)  & 
               +   faka * (_cou[0][34][0) - fakaac * _cou[0][34][1) ) + fakac * _cou[1][18][1)
          _cou[4][34][0) = _wmp[1]*_cou[1][34][1) +  fakac * _cou[1][16][1)
          _cou[5][34][0) = _wmp[2]*_cou[1][34][1) +  fakac2 * _cou[1][19][1)
          _cou[8][20][0) = _wmp[1]*_cou[2][20][1)  & 
               +    faka * (_cou[0][20][0) - fakaac * _cou[0][20][1) ) 
          _cou[6][20][0) = _wmp[2]*_cou[2][20][1) 
          _cou[8][21][0) = _wmp[1]*_cou[2][21][1)  & 
               +  faka * (_cou[0][21][0) - fakaac * _cou[0][21][1) ) + fakac4 * _cou[2][11][1)
          _cou[6][21][0) = _wmp[2]*_cou[2][21][1) 
          _cou[8][22][0) = _wmp[1]*_cou[2][22][1)  & 
               +    faka * (_cou[0][22][0) - fakaac * _cou[0][22][1) ) 
          _cou[6][22][0) = _wmp[2]*_cou[2][22][1) +  fakac4 * _cou[2][12][1)
          _cou[8][23][0) = _wmp[1]*_cou[2][23][1)  & 
               +   faka * (_cou[0][23][0) - fakaac * _cou[0][23][1) ) + fakac * _cou[2][10][1)
          _cou[6][23][0) = _wmp[2]*_cou[2][23][1) 
          _cou[8][24][0) = _wmp[1]*_cou[2][24][1)  & 
               +  faka * (_cou[0][24][0) - fakaac * _cou[0][24][1) ) + fakac3 * _cou[2][14][1)
          _cou[6][24][0) = _wmp[2]*_cou[2][24][1) 
          _cou[8][25][0) = _wmp[1]*_cou[2][25][1)  & 
               +    faka * (_cou[0][25][0) - fakaac * _cou[0][25][1) ) 
          _cou[6][25][0) = _wmp[2]*_cou[2][25][1) +  fakac * _cou[2][10][1)
          _cou[8][26][0) = _wmp[1]*_cou[2][26][1)  & 
               +    faka * (_cou[0][26][0) - fakaac * _cou[0][26][1) ) 
          _cou[6][26][0) = _wmp[2]*_cou[2][26][1) +  fakac3 * _cou[2][16][1)
          _cou[8][27][0) = _wmp[1]*_cou[2][27][1)  & 
               +  faka * (_cou[0][27][0) - fakaac * _cou[0][27][1) ) + fakac3 * _cou[2][17][1)
          _cou[6][27][0) = _wmp[2]*_cou[2][27][1) +  fakac * _cou[2][11][1)
          _cou[8][28][0) = _wmp[1]*_cou[2][28][1)  & 
               +   faka * (_cou[0][28][0) - fakaac * _cou[0][28][1) ) + fakac * _cou[2][12][1)
          _cou[6][28][0) = _wmp[2]*_cou[2][28][1) +  fakac3 * _cou[2][18][1)
          _cou[8][29][0) = _wmp[1]*_cou[2][29][1)  & 
               +    faka * (_cou[0][29][0) - fakaac * _cou[0][29][1) ) 
          _cou[6][29][0) = _wmp[2]*_cou[2][29][1) +  fakac2 * _cou[2][15][1)
          _cou[8][30][0) = _wmp[1]*_cou[2][30][1)  & 
               +  faka * (_cou[0][30][0) - fakaac * _cou[0][30][1) ) + fakac2 * _cou[2][18][1)
          _cou[6][30][0) = _wmp[2]*_cou[2][30][1) +  fakac2 * _cou[2][17][1)
          _cou[8][31][0) = _wmp[1]*_cou[2][31][1)  & 
               +  faka * (_cou[0][31][0) - fakaac * _cou[0][31][1) ) + fakac2 * _cou[2][13][1)
          _cou[6][31][0) = _wmp[2]*_cou[2][31][1) 
          _cou[8][32][0) = _wmp[1]*_cou[2][32][1)  & 
               +  faka * (_cou[0][32][0) - fakaac * _cou[0][32][1) ) + fakac * _cou[2][15][1)
          _cou[6][32][0) = _wmp[2]*_cou[2][32][1) +  fakac * _cou[2][13][1)
          _cou[8][33][0) = _wmp[1]*_cou[2][33][1)  & 
               +  faka * (_cou[0][33][0) - fakaac * _cou[0][33][1) ) + fakac2 * _cou[2][19][1)
          _cou[6][33][0) = _wmp[2]*_cou[2][33][1) +  fakac * _cou[2][14][1)
          _cou[8][34][0) = _wmp[1]*_cou[2][34][1)  & 
               +   faka * (_cou[0][34][0) - fakaac * _cou[0][34][1) ) + fakac * _cou[2][16][1)
          _cou[6][34][0) = _wmp[2]*_cou[2][34][1) +  fakac2 * _cou[2][19][1)
          _cou[9][20][0) = _wmp[2]*_cou[3][20][1)  & 
               +    faka * (_cou[0][20][0) - fakaac * _cou[0][20][1) ) 
          _cou[9][21][0) = _wmp[2]*_cou[3][21][1)  & 
               +    faka * (_cou[0][21][0) - fakaac * _cou[0][21][1) ) 
          _cou[9][22][0) = _wmp[2]*_cou[3][22][1)  & 
               +  faka * (_cou[0][22][0) - fakaac * _cou[0][22][1) ) + fakac4 * _cou[3][12][1)
          _cou[9][23][0) = _wmp[2]*_cou[3][23][1)  & 
               +    faka * (_cou[0][23][0) - fakaac * _cou[0][23][1) ) 
          _cou[9][24][0) = _wmp[2]*_cou[3][24][1)  & 
               +    faka * (_cou[0][24][0) - fakaac * _cou[0][24][1) ) 
          _cou[9][25][0) = _wmp[2]*_cou[3][25][1)  & 
               +   faka * (_cou[0][25][0) - fakaac * _cou[0][25][1) ) + fakac * _cou[3][10][1)
          _cou[9][26][0) = _wmp[2]*_cou[3][26][1)  & 
               +  faka * (_cou[0][26][0) - fakaac * _cou[0][26][1) ) + fakac3 * _cou[3][16][1)
          _cou[9][27][0) = _wmp[2]*_cou[3][27][1)  & 
               +  faka * (_cou[0][27][0) - fakaac * _cou[0][27][1) ) + fakac * _cou[3][11][1)
          _cou[9][28][0) = _wmp[2]*_cou[3][28][1)  & 
               +  faka * (_cou[0][28][0) - fakaac * _cou[0][28][1) ) + fakac3 * _cou[3][18][1)
          _cou[9][29][0) = _wmp[2]*_cou[3][29][1)  & 
               +  faka * (_cou[0][29][0) - fakaac * _cou[0][29][1) ) + fakac2 * _cou[3][15][1)
          _cou[9][30][0) = _wmp[2]*_cou[3][30][1)  & 
               +  faka * (_cou[0][30][0) - fakaac * _cou[0][30][1) ) + fakac2 * _cou[3][17][1)
          _cou[9][31][0) = _wmp[2]*_cou[3][31][1)  & 
               +    faka * (_cou[0][31][0) - fakaac * _cou[0][31][1) ) 
          _cou[9][32][0) = _wmp[2]*_cou[3][32][1)  & 
               +   faka * (_cou[0][32][0) - fakaac * _cou[0][32][1) ) + fakac * _cou[3][13][1)
          _cou[9][33][0) = _wmp[2]*_cou[3][33][1)  & 
               +   faka * (_cou[0][33][0) - fakaac * _cou[0][33][1) ) + fakac * _cou[3][14][1)
          _cou[9][34][0) = _wmp[2]*_cou[3][34][1)  & 
               +  faka * (_cou[0][34][0) - fakaac * _cou[0][34][1) ) + fakac2 * _cou[3][19][1)


          _cou[20][4][0) = _wmp[0]*_cou[10][4][1)  & 
               +  faka3 * (_cou[7][4][0) - fakaac * _cou[7][4][1) ) + fakac * _cou[10][2][1)
          _cou[23][4][0) = _wmp[1]*_cou[10][4][1) +  fakac * _cou[10][1][1)
          _cou[25][4][0) = _wmp[2]*_cou[10][4][1) 
          _cou[20][5][0) = _wmp[0]*_cou[10][5][1)  & 
               +  faka3 * (_cou[7][5][0) - fakaac * _cou[7][5][1) ) + fakac * _cou[10][3][1)
          _cou[23][5][0) = _wmp[1]*_cou[10][5][1) 
          _cou[25][5][0) = _wmp[2]*_cou[10][5][1) +  fakac * _cou[10][1][1)
          _cou[20][6][0) = _wmp[0]*_cou[10][6][1)  & 
               +    faka3 * (_cou[7][6][0) - fakaac * _cou[7][6][1) ) 
          _cou[23][6][0) = _wmp[1]*_cou[10][6][1) +  fakac * _cou[10][3][1)
          _cou[25][6][0) = _wmp[2]*_cou[10][6][1) +  fakac * _cou[10][2][1)
          _cou[20][7][0) = _wmp[0]*_cou[10][7][1)  & 
               + faka3 * (_cou[7][7][0) - fakaac * _cou[7][7][1) ) + fakac2 * _cou[10][1][1)
          _cou[23][7][0) = _wmp[1]*_cou[10][7][1) 
          _cou[25][7][0) = _wmp[2]*_cou[10][7][1) 
          _cou[20][8][0) = _wmp[0]*_cou[10][8][1)  & 
               +    faka3 * (_cou[7][8][0) - fakaac * _cou[7][8][1) ) 
          _cou[23][8][0) = _wmp[1]*_cou[10][8][1) +  fakac2 * _cou[10][2][1)
          _cou[25][8][0) = _wmp[2]*_cou[10][8][1) 
          _cou[20][9][0) = _wmp[0]*_cou[10][9][1)  & 
               +    faka3 * (_cou[7][9][0) - fakaac * _cou[7][9][1) ) 
          _cou[23][9][0) = _wmp[1]*_cou[10][9][1) 
          _cou[25][9][0) = _wmp[2]*_cou[10][9][1) +  fakac2 * _cou[10][3][1)
          _cou[24][4][0) = _wmp[0]*_cou[11][4][1) +  fakac * _cou[11][2][1)
          _cou[21][4][0) = _wmp[1]*_cou[11][4][1)  & 
               +  faka3 * (_cou[8][4][0) - fakaac * _cou[8][4][1) ) + fakac * _cou[11][1][1)
          _cou[27][4][0) = _wmp[2]*_cou[11][4][1) 
          _cou[24][5][0) = _wmp[0]*_cou[11][5][1) +  fakac * _cou[11][3][1)
          _cou[21][5][0) = _wmp[1]*_cou[11][5][1)  & 
               +    faka3 * (_cou[8][5][0) - fakaac * _cou[8][5][1) ) 
          _cou[27][5][0) = _wmp[2]*_cou[11][5][1) +  fakac * _cou[11][1][1)
          _cou[24][6][0) = _wmp[0]*_cou[11][6][1) 
          _cou[21][6][0) = _wmp[1]*_cou[11][6][1)  & 
               +  faka3 * (_cou[8][6][0) - fakaac * _cou[8][6][1) ) + fakac * _cou[11][3][1)
          _cou[27][6][0) = _wmp[2]*_cou[11][6][1) +  fakac * _cou[11][2][1)
          _cou[24][7][0) = _wmp[0]*_cou[11][7][1) +  fakac2 * _cou[11][1][1)
          _cou[21][7][0) = _wmp[1]*_cou[11][7][1)  & 
               +    faka3 * (_cou[8][7][0) - fakaac * _cou[8][7][1) ) 
          _cou[27][7][0) = _wmp[2]*_cou[11][7][1) 
          _cou[24][8][0) = _wmp[0]*_cou[11][8][1) 
          _cou[21][8][0) = _wmp[1]*_cou[11][8][1)  & 
               + faka3 * (_cou[8][8][0) - fakaac * _cou[8][8][1) ) + fakac2 * _cou[11][2][1)
          _cou[27][8][0) = _wmp[2]*_cou[11][8][1) 
          _cou[24][9][0) = _wmp[0]*_cou[11][9][1) 
          _cou[21][9][0) = _wmp[1]*_cou[11][9][1)  & 
               +    faka3 * (_cou[8][9][0) - fakaac * _cou[8][9][1) ) 
          _cou[27][9][0) = _wmp[2]*_cou[11][9][1) +  fakac2 * _cou[11][3][1)
          _cou[26][4][0) = _wmp[0]*_cou[12][4][1) +  fakac * _cou[12][2][1)
          _cou[28][4][0) = _wmp[1]*_cou[12][4][1) +  fakac * _cou[12][1][1)
          _cou[22][4][0) = _wmp[2]*_cou[12][4][1)  & 
               +    faka3 * (_cou[9][4][0) - fakaac * _cou[9][4][1) ) 
          _cou[26][5][0) = _wmp[0]*_cou[12][5][1) +  fakac * _cou[12][3][1)
          _cou[28][5][0) = _wmp[1]*_cou[12][5][1) 
          _cou[22][5][0) = _wmp[2]*_cou[12][5][1)  & 
               +  faka3 * (_cou[9][5][0) - fakaac * _cou[9][5][1) ) + fakac * _cou[12][1][1)
          _cou[26][6][0) = _wmp[0]*_cou[12][6][1) 
          _cou[28][6][0) = _wmp[1]*_cou[12][6][1) +  fakac * _cou[12][3][1)
          _cou[22][6][0) = _wmp[2]*_cou[12][6][1)  & 
               +  faka3 * (_cou[9][6][0) - fakaac * _cou[9][6][1) ) + fakac * _cou[12][2][1)
          _cou[26][7][0) = _wmp[0]*_cou[12][7][1) +  fakac2 * _cou[12][1][1)
          _cou[28][7][0) = _wmp[1]*_cou[12][7][1) 
          _cou[22][7][0) = _wmp[2]*_cou[12][7][1)  & 
               +    faka3 * (_cou[9][7][0) - fakaac * _cou[9][7][1) ) 
          _cou[26][8][0) = _wmp[0]*_cou[12][8][1) 
          _cou[28][8][0) = _wmp[1]*_cou[12][8][1) +  fakac2 * _cou[12][2][1)
          _cou[22][8][0) = _wmp[2]*_cou[12][8][1)  & 
               +    faka3 * (_cou[9][8][0) - fakaac * _cou[9][8][1) ) 
          _cou[26][9][0) = _wmp[0]*_cou[12][9][1) 
          _cou[28][9][0) = _wmp[1]*_cou[12][9][1) 
          _cou[22][9][0) = _wmp[2]*_cou[12][9][1)  & 
               + faka3 * (_cou[9][9][0) - fakaac * _cou[9][9][1) ) + fakac2 * _cou[12][3][1)
          _cou[31][4][0) = _wmp[1]*_cou[13][4][1)  & 
               +   faka * (_cou[7][4][0) - fakaac * _cou[7][4][1) ) + fakac * _cou[13][1][1)
          _cou[32][4][0) = _wmp[2]*_cou[13][4][1) 
          _cou[31][5][0) = _wmp[1]*_cou[13][5][1)  & 
               +    faka * (_cou[7][5][0) - fakaac * _cou[7][5][1) ) 
          _cou[32][5][0) = _wmp[2]*_cou[13][5][1) +  fakac * _cou[13][1][1)
          _cou[31][6][0) = _wmp[1]*_cou[13][6][1)  & 
               +   faka * (_cou[7][6][0) - fakaac * _cou[7][6][1) ) + fakac * _cou[13][3][1)
          _cou[32][6][0) = _wmp[2]*_cou[13][6][1) +  fakac * _cou[13][2][1)
          _cou[31][7][0) = _wmp[1]*_cou[13][7][1)  & 
               +    faka * (_cou[7][7][0) - fakaac * _cou[7][7][1) ) 
          _cou[32][7][0) = _wmp[2]*_cou[13][7][1) 
          _cou[31][8][0) = _wmp[1]*_cou[13][8][1)  & 
               +  faka * (_cou[7][8][0) - fakaac * _cou[7][8][1) ) + fakac2 * _cou[13][2][1)
          _cou[32][8][0) = _wmp[2]*_cou[13][8][1) 
          _cou[31][9][0) = _wmp[1]*_cou[13][9][1)  & 
               +    faka * (_cou[7][9][0) - fakaac * _cou[7][9][1) ) 
          _cou[32][9][0) = _wmp[2]*_cou[13][9][1) +  fakac2 * _cou[13][3][1)
          _cou[33][4][0) = _wmp[2]*_cou[14][4][1) 
          _cou[33][5][0) = _wmp[2]*_cou[14][5][1) +  fakac * _cou[14][1][1)
          _cou[33][6][0) = _wmp[2]*_cou[14][6][1) +  fakac * _cou[14][2][1)
          _cou[33][7][0) = _wmp[2]*_cou[14][7][1) 
          _cou[33][8][0) = _wmp[2]*_cou[14][8][1) 
          _cou[33][9][0) = _wmp[2]*_cou[14][9][1) +  fakac2 * _cou[14][3][1)
          _cou[29][4][0) = _wmp[2]*_cou[15][4][1)  & 
               +    faka * (_cou[7][4][0) - fakaac * _cou[7][4][1) ) 
          _cou[29][5][0) = _wmp[2]*_cou[15][5][1)  & 
               +   faka * (_cou[7][5][0) - fakaac * _cou[7][5][1) ) + fakac * _cou[15][1][1)
          _cou[29][6][0) = _wmp[2]*_cou[15][6][1)  & 
               +   faka * (_cou[7][6][0) - fakaac * _cou[7][6][1) ) + fakac * _cou[15][2][1)
          _cou[29][7][0) = _wmp[2]*_cou[15][7][1)  & 
               +    faka * (_cou[7][7][0) - fakaac * _cou[7][7][1) ) 
          _cou[29][8][0) = _wmp[2]*_cou[15][8][1)  & 
               +    faka * (_cou[7][8][0) - fakaac * _cou[7][8][1) ) 
          _cou[29][9][0) = _wmp[2]*_cou[15][9][1)  & 
               +  faka * (_cou[7][9][0) - fakaac * _cou[7][9][1) ) + fakac2 * _cou[15][3][1)
          _cou[34][4][0) = _wmp[1]*_cou[16][4][1) +  fakac * _cou[16][1][1)
          _cou[34][5][0) = _wmp[1]*_cou[16][5][1) 
          _cou[34][6][0) = _wmp[1]*_cou[16][6][1) +  fakac * _cou[16][3][1)
          _cou[34][7][0) = _wmp[1]*_cou[16][7][1) 
          _cou[34][8][0) = _wmp[1]*_cou[16][8][1) +  fakac2 * _cou[16][2][1)
          _cou[34][9][0) = _wmp[1]*_cou[16][9][1) 
          _cou[30][4][0) = _wmp[2]*_cou[17][4][1)  & 
               +    faka * (_cou[8][4][0) - fakaac * _cou[8][4][1) ) 
          _cou[30][5][0) = _wmp[2]*_cou[17][5][1)  & 
               +   faka * (_cou[8][5][0) - fakaac * _cou[8][5][1) ) + fakac * _cou[17][1][1)
          _cou[30][6][0) = _wmp[2]*_cou[17][6][1)  & 
               +   faka * (_cou[8][6][0) - fakaac * _cou[8][6][1) ) + fakac * _cou[17][2][1)
          _cou[30][7][0) = _wmp[2]*_cou[17][7][1)  & 
               +    faka * (_cou[8][7][0) - fakaac * _cou[8][7][1) ) 
          _cou[30][8][0) = _wmp[2]*_cou[17][8][1)  & 
               +    faka * (_cou[8][8][0) - fakaac * _cou[8][8][1) ) 
          _cou[30][9][0) = _wmp[2]*_cou[17][9][1)  & 
               +  faka * (_cou[8][9][0) - fakaac * _cou[8][9][1) ) + fakac2 * _cou[17][3][1)

       endif


       if(abs(lmax1)+abs(lmax2).ge.7) then

          _cou[1][0][6) = _wmp[0]*_cou[0][0][7) 
          _cou[0][1][6) = _wmq[0]*_cou[0][0][7) 
          _cou[2][0][6) = _wmp[1]*_cou[0][0][7) 
          _cou[0][2][6) = _wmq[1]*_cou[0][0][7) 
          _cou[3][0][6) = _wmp[2]*_cou[0][0][7) 
          _cou[0][3][6) = _wmq[2]*_cou[0][0][7) 

          _cou[1][1][5) = _wmp[0]*_cou[0][1][6) +  fakac * _cou[0][0][6)
          _cou[0][7][5) = _wmq[0]*_cou[0][1][6)  & 
               +    fakc * (_cou[0][0][5) - fakaca * _cou[0][0][6) ) 
          _cou[2][1][5) = _wmp[1]*_cou[0][1][6) 
          _cou[0][4][5) = _wmq[1]*_cou[0][1][6) 
          _cou[3][1][5) = _wmp[2]*_cou[0][1][6) 
          _cou[0][5][5) = _wmq[2]*_cou[0][1][6) 
          _cou[1][2][5) = _wmp[0]*_cou[0][2][6) 
          _cou[2][2][5) = _wmp[1]*_cou[0][2][6) +  fakac * _cou[0][0][6)
          _cou[0][8][5) = _wmq[1]*_cou[0][2][6)  & 
               +    fakc * (_cou[0][0][5) - fakaca * _cou[0][0][6) ) 
          _cou[3][2][5) = _wmp[2]*_cou[0][2][6) 
          _cou[0][6][5) = _wmq[2]*_cou[0][2][6) 
          _cou[1][3][5) = _wmp[0]*_cou[0][3][6) 
          _cou[2][3][5) = _wmp[1]*_cou[0][3][6) 
          _cou[3][3][5) = _wmp[2]*_cou[0][3][6) +  fakac * _cou[0][0][6)
          _cou[0][9][5) = _wmq[2]*_cou[0][3][6)  & 
               +    fakc * (_cou[0][0][5) - fakaca * _cou[0][0][6) ) 
          _cou[7][0][5) = _wmp[0]*_cou[1][0][6)  & 
               +    faka * (_cou[0][0][5) - fakaac * _cou[0][0][6) ) 
          _cou[4][0][5) = _wmp[1]*_cou[1][0][6) 
          _cou[5][0][5) = _wmp[2]*_cou[1][0][6) 
          _cou[8][0][5) = _wmp[1]*_cou[2][0][6)  & 
               +    faka * (_cou[0][0][5) - fakaac * _cou[0][0][6) ) 
          _cou[6][0][5) = _wmp[2]*_cou[2][0][6) 
          _cou[9][0][5) = _wmp[2]*_cou[3][0][6)  & 
               +    faka * (_cou[0][0][5) - fakaac * _cou[0][0][6) ) 
          _cou[1][4][4) = _wmp[0]*_cou[0][4][5) +  fakac * _cou[0][2][5)
          _cou[0][13][4) = _wmq[0]*_cou[0][4][5)  & 
               +    fakc * (_cou[0][2][4) - fakaca * _cou[0][2][5) ) 
          _cou[2][4][4) = _wmp[1]*_cou[0][4][5) +  fakac * _cou[0][1][5)
          _cou[0][14][4) = _wmq[1]*_cou[0][4][5)  & 
               +    fakc * (_cou[0][1][4) - fakaca * _cou[0][1][5) ) 
          _cou[3][4][4) = _wmp[2]*_cou[0][4][5) 
          _cou[0][19][4) = _wmq[2]*_cou[0][4][5) 
          _cou[1][5][4) = _wmp[0]*_cou[0][5][5) +  fakac * _cou[0][3][5)
          _cou[0][15][4) = _wmq[0]*_cou[0][5][5)  & 
               +    fakc * (_cou[0][3][4) - fakaca * _cou[0][3][5) ) 
          _cou[2][5][4) = _wmp[1]*_cou[0][5][5) 
          _cou[3][5][4) = _wmp[2]*_cou[0][5][5) +  fakac * _cou[0][1][5)
          _cou[0][16][4) = _wmq[2]*_cou[0][5][5)  & 
               +    fakc * (_cou[0][1][4) - fakaca * _cou[0][1][5) ) 
          _cou[1][6][4) = _wmp[0]*_cou[0][6][5) 
          _cou[2][6][4) = _wmp[1]*_cou[0][6][5) +  fakac * _cou[0][3][5)
          _cou[0][17][4) = _wmq[1]*_cou[0][6][5)  & 
               +    fakc * (_cou[0][3][4) - fakaca * _cou[0][3][5) ) 
          _cou[3][6][4) = _wmp[2]*_cou[0][6][5) +  fakac * _cou[0][2][5)
          _cou[0][18][4) = _wmq[2]*_cou[0][6][5)  & 
               +    fakc * (_cou[0][2][4) - fakaca * _cou[0][2][5) ) 
          _cou[1][7][4) = _wmp[0]*_cou[0][7][5) +  fakac2 * _cou[0][1][5)
          _cou[0][10][4) = _wmq[0]*_cou[0][7][5)  & 
               +    fakc2 * (_cou[0][1][4) - fakaca * _cou[0][1][5) ) 
          _cou[2][7][4) = _wmp[1]*_cou[0][7][5) 
          _cou[3][7][4) = _wmp[2]*_cou[0][7][5) 
          _cou[1][8][4) = _wmp[0]*_cou[0][8][5) 
          _cou[2][8][4) = _wmp[1]*_cou[0][8][5) +  fakac2 * _cou[0][2][5)
          _cou[0][11][4) = _wmq[1]*_cou[0][8][5)  & 
               +    fakc2 * (_cou[0][2][4) - fakaca * _cou[0][2][5) ) 
          _cou[3][8][4) = _wmp[2]*_cou[0][8][5) 
          _cou[1][9][4) = _wmp[0]*_cou[0][9][5) 
          _cou[2][9][4) = _wmp[1]*_cou[0][9][5) 
          _cou[3][9][4) = _wmp[2]*_cou[0][9][5) +  fakac2 * _cou[0][3][5)
          _cou[0][12][4) = _wmq[2]*_cou[0][9][5)  & 
               +    fakc2 * (_cou[0][3][4) - fakaca * _cou[0][3][5) ) 
          _cou[7][1][4) = _wmp[0]*_cou[1][1][5)  & 
               +   faka * (_cou[0][1][4) - fakaac * _cou[0][1][5) ) + fakac * _cou[1][0][5)
          _cou[4][1][4) = _wmp[1]*_cou[1][1][5) 
          _cou[5][1][4) = _wmp[2]*_cou[1][1][5) 
          _cou[7][2][4) = _wmp[0]*_cou[1][2][5)  & 
               +    faka * (_cou[0][2][4) - fakaac * _cou[0][2][5) ) 
          _cou[4][2][4) = _wmp[1]*_cou[1][2][5) +  fakac * _cou[1][0][5)
          _cou[5][2][4) = _wmp[2]*_cou[1][2][5) 
          _cou[7][3][4) = _wmp[0]*_cou[1][3][5)  & 
               +    faka * (_cou[0][3][4) - fakaac * _cou[0][3][5) ) 
          _cou[4][3][4) = _wmp[1]*_cou[1][3][5) 
          _cou[5][3][4) = _wmp[2]*_cou[1][3][5) +  fakac * _cou[1][0][5)
          _cou[8][1][4) = _wmp[1]*_cou[2][1][5)  & 
               +    faka * (_cou[0][1][4) - fakaac * _cou[0][1][5) ) 
          _cou[6][1][4) = _wmp[2]*_cou[2][1][5) 
          _cou[8][2][4) = _wmp[1]*_cou[2][2][5)  & 
               +   faka * (_cou[0][2][4) - fakaac * _cou[0][2][5) ) + fakac * _cou[2][0][5)
          _cou[6][2][4) = _wmp[2]*_cou[2][2][5) 
          _cou[8][3][4) = _wmp[1]*_cou[2][3][5)  & 
               +    faka * (_cou[0][3][4) - fakaac * _cou[0][3][5) ) 
          _cou[6][3][4) = _wmp[2]*_cou[2][3][5) +  fakac * _cou[2][0][5)
          _cou[9][1][4) = _wmp[2]*_cou[3][1][5)  & 
               +    faka * (_cou[0][1][4) - fakaac * _cou[0][1][5) ) 
          _cou[9][2][4) = _wmp[2]*_cou[3][2][5)  & 
               +    faka * (_cou[0][2][4) - fakaac * _cou[0][2][5) ) 
          _cou[9][3][4) = _wmp[2]*_cou[3][3][5)  & 
               +   faka * (_cou[0][3][4) - fakaac * _cou[0][3][5) ) + fakac * _cou[3][0][5)
          _cou[13][0][4) = _wmp[0]*_cou[4][0][5)  & 
               +    faka * (_cou[2][0][4) - fakaac * _cou[2][0][5) ) 
          _cou[14][0][4) = _wmp[1]*_cou[4][0][5)  & 
               +    faka * (_cou[1][0][4) - fakaac * _cou[1][0][5) ) 
          _cou[19][0][4) = _wmp[2]*_cou[4][0][5) 
          _cou[15][0][4) = _wmp[0]*_cou[5][0][5)  & 
               +    faka * (_cou[3][0][4) - fakaac * _cou[3][0][5) ) 
          _cou[16][0][4) = _wmp[2]*_cou[5][0][5)  & 
               +    faka * (_cou[1][0][4) - fakaac * _cou[1][0][5) ) 
          _cou[17][0][4) = _wmp[1]*_cou[6][0][5)  & 
               +    faka * (_cou[3][0][4) - fakaac * _cou[3][0][5) ) 
          _cou[18][0][4) = _wmp[2]*_cou[6][0][5)  & 
               +    faka * (_cou[2][0][4) - fakaac * _cou[2][0][5) ) 
          _cou[10][0][4) = _wmp[0]*_cou[7][0][5)  & 
               +  faka2 * (_cou[1][0][4) - fakaac * _cou[1][0][5) ) 
          _cou[11][0][4) = _wmp[1]*_cou[8][0][5)  & 
               +  faka2 * (_cou[2][0][4) - fakaac * _cou[2][0][5) ) 
          _cou[12][0][4) = _wmp[2]*_cou[9][0][5)  & 
               +  faka2 * (_cou[3][0][4) - fakaac * _cou[3][0][5) ) 
          _cou[1][10][3) = _wmp[0]*_cou[0][10][4) +  fakac3 * _cou[0][7][4)
          _cou[0][20][3) = _wmq[0]*_cou[0][10][4)  & 
               +    fakc3 * (_cou[0][7][3) - fakaca * _cou[0][7][4) ) 
          _cou[2][10][3) = _wmp[1]*_cou[0][10][4) 
          _cou[0][23][3) = _wmq[1]*_cou[0][10][4) 
          _cou[3][10][3) = _wmp[2]*_cou[0][10][4) 
          _cou[0][25][3) = _wmq[2]*_cou[0][10][4) 
          _cou[1][11][3) = _wmp[0]*_cou[0][11][4) 
          _cou[0][24][3) = _wmq[0]*_cou[0][11][4) 
          _cou[2][11][3) = _wmp[1]*_cou[0][11][4) +  fakac3 * _cou[0][8][4)
          _cou[0][21][3) = _wmq[1]*_cou[0][11][4)  & 
               +    fakc3 * (_cou[0][8][3) - fakaca * _cou[0][8][4) ) 
          _cou[3][11][3) = _wmp[2]*_cou[0][11][4) 
          _cou[0][27][3) = _wmq[2]*_cou[0][11][4) 
          _cou[1][12][3) = _wmp[0]*_cou[0][12][4) 
          _cou[0][26][3) = _wmq[0]*_cou[0][12][4) 
          _cou[2][12][3) = _wmp[1]*_cou[0][12][4) 
          _cou[0][28][3) = _wmq[1]*_cou[0][12][4) 
          _cou[3][12][3) = _wmp[2]*_cou[0][12][4) +  fakac3 * _cou[0][9][4)
          _cou[0][22][3) = _wmq[2]*_cou[0][12][4)  & 
               +    fakc3 * (_cou[0][9][3) - fakaca * _cou[0][9][4) ) 
          _cou[1][13][3) = _wmp[0]*_cou[0][13][4) +  fakac2 * _cou[0][4][4)
          _cou[2][13][3) = _wmp[1]*_cou[0][13][4) +  fakac * _cou[0][7][4)
          _cou[0][31][3) = _wmq[1]*_cou[0][13][4)  & 
               +    fakc * (_cou[0][7][3) - fakaca * _cou[0][7][4) ) 
          _cou[3][13][3) = _wmp[2]*_cou[0][13][4) 
          _cou[0][32][3) = _wmq[2]*_cou[0][13][4) 
          _cou[1][14][3) = _wmp[0]*_cou[0][14][4) +  fakac * _cou[0][8][4)
          _cou[2][14][3) = _wmp[1]*_cou[0][14][4) +  fakac2 * _cou[0][4][4)
          _cou[3][14][3) = _wmp[2]*_cou[0][14][4) 
          _cou[0][33][3) = _wmq[2]*_cou[0][14][4) 
          _cou[1][15][3) = _wmp[0]*_cou[0][15][4) +  fakac2 * _cou[0][5][4)
          _cou[2][15][3) = _wmp[1]*_cou[0][15][4) 
          _cou[3][15][3) = _wmp[2]*_cou[0][15][4) +  fakac * _cou[0][7][4)
          _cou[0][29][3) = _wmq[2]*_cou[0][15][4)  & 
               +    fakc * (_cou[0][7][3) - fakaca * _cou[0][7][4) ) 
          _cou[1][16][3) = _wmp[0]*_cou[0][16][4) +  fakac * _cou[0][9][4)
          _cou[2][16][3) = _wmp[1]*_cou[0][16][4) 
          _cou[0][34][3) = _wmq[1]*_cou[0][16][4) 
          _cou[3][16][3) = _wmp[2]*_cou[0][16][4) +  fakac2 * _cou[0][5][4)
          _cou[1][17][3) = _wmp[0]*_cou[0][17][4) 
          _cou[2][17][3) = _wmp[1]*_cou[0][17][4) +  fakac2 * _cou[0][6][4)
          _cou[3][17][3) = _wmp[2]*_cou[0][17][4) +  fakac * _cou[0][8][4)
          _cou[0][30][3) = _wmq[2]*_cou[0][17][4)  & 
               +    fakc * (_cou[0][8][3) - fakaca * _cou[0][8][4) ) 
          _cou[1][18][3) = _wmp[0]*_cou[0][18][4) 
          _cou[2][18][3) = _wmp[1]*_cou[0][18][4) +  fakac * _cou[0][9][4)
          _cou[3][18][3) = _wmp[2]*_cou[0][18][4) +  fakac2 * _cou[0][6][4)
          _cou[1][19][3) = _wmp[0]*_cou[0][19][4) +  fakac * _cou[0][6][4)
          _cou[2][19][3) = _wmp[1]*_cou[0][19][4) +  fakac * _cou[0][5][4)
          _cou[3][19][3) = _wmp[2]*_cou[0][19][4) +  fakac * _cou[0][4][4)
          _cou[7][4][3) = _wmp[0]*_cou[1][4][4)  & 
               +   faka * (_cou[0][4][3) - fakaac * _cou[0][4][4) ) + fakac * _cou[1][2][4)
          _cou[4][4][3) = _wmp[1]*_cou[1][4][4) +  fakac * _cou[1][1][4)
          _cou[5][4][3) = _wmp[2]*_cou[1][4][4) 
          _cou[7][5][3) = _wmp[0]*_cou[1][5][4)  & 
               +   faka * (_cou[0][5][3) - fakaac * _cou[0][5][4) ) + fakac * _cou[1][3][4)
          _cou[4][5][3) = _wmp[1]*_cou[1][5][4) 
          _cou[5][5][3) = _wmp[2]*_cou[1][5][4) +  fakac * _cou[1][1][4)
          _cou[7][6][3) = _wmp[0]*_cou[1][6][4)  & 
               +    faka * (_cou[0][6][3) - fakaac * _cou[0][6][4) ) 
          _cou[4][6][3) = _wmp[1]*_cou[1][6][4) +  fakac * _cou[1][3][4)
          _cou[5][6][3) = _wmp[2]*_cou[1][6][4) +  fakac * _cou[1][2][4)
          _cou[7][7][3) = _wmp[0]*_cou[1][7][4)  & 
               +  faka * (_cou[0][7][3) - fakaac * _cou[0][7][4) ) + fakac2 * _cou[1][1][4)
          _cou[4][7][3) = _wmp[1]*_cou[1][7][4) 
          _cou[5][7][3) = _wmp[2]*_cou[1][7][4) 
          _cou[7][8][3) = _wmp[0]*_cou[1][8][4)  & 
               +    faka * (_cou[0][8][3) - fakaac * _cou[0][8][4) ) 
          _cou[4][8][3) = _wmp[1]*_cou[1][8][4) +  fakac2 * _cou[1][2][4)
          _cou[5][8][3) = _wmp[2]*_cou[1][8][4) 
          _cou[7][9][3) = _wmp[0]*_cou[1][9][4)  & 
               +    faka * (_cou[0][9][3) - fakaac * _cou[0][9][4) ) 
          _cou[4][9][3) = _wmp[1]*_cou[1][9][4) 
          _cou[5][9][3) = _wmp[2]*_cou[1][9][4) +  fakac2 * _cou[1][3][4)
          _cou[8][4][3) = _wmp[1]*_cou[2][4][4)  & 
               +   faka * (_cou[0][4][3) - fakaac * _cou[0][4][4) ) + fakac * _cou[2][1][4)
          _cou[6][4][3) = _wmp[2]*_cou[2][4][4) 
          _cou[8][5][3) = _wmp[1]*_cou[2][5][4)  & 
               +    faka * (_cou[0][5][3) - fakaac * _cou[0][5][4) ) 
          _cou[6][5][3) = _wmp[2]*_cou[2][5][4) +  fakac * _cou[2][1][4)
          _cou[8][6][3) = _wmp[1]*_cou[2][6][4)  & 
               +   faka * (_cou[0][6][3) - fakaac * _cou[0][6][4) ) + fakac * _cou[2][3][4)
          _cou[6][6][3) = _wmp[2]*_cou[2][6][4) +  fakac * _cou[2][2][4)
          _cou[8][7][3) = _wmp[1]*_cou[2][7][4)  & 
               +    faka * (_cou[0][7][3) - fakaac * _cou[0][7][4) ) 
          _cou[6][7][3) = _wmp[2]*_cou[2][7][4) 
          _cou[8][8][3) = _wmp[1]*_cou[2][8][4)  & 
               +  faka * (_cou[0][8][3) - fakaac * _cou[0][8][4) ) + fakac2 * _cou[2][2][4)
          _cou[6][8][3) = _wmp[2]*_cou[2][8][4) 
          _cou[8][9][3) = _wmp[1]*_cou[2][9][4)  & 
               +    faka * (_cou[0][9][3) - fakaac * _cou[0][9][4) ) 
          _cou[6][9][3) = _wmp[2]*_cou[2][9][4) +  fakac2 * _cou[2][3][4)
          _cou[9][4][3) = _wmp[2]*_cou[3][4][4)  & 
               +    faka * (_cou[0][4][3) - fakaac * _cou[0][4][4) ) 
          _cou[9][5][3) = _wmp[2]*_cou[3][5][4)  & 
               +   faka * (_cou[0][5][3) - fakaac * _cou[0][5][4) ) + fakac * _cou[3][1][4)
          _cou[9][6][3) = _wmp[2]*_cou[3][6][4)  & 
               +   faka * (_cou[0][6][3) - fakaac * _cou[0][6][4) ) + fakac * _cou[3][2][4)
          _cou[9][7][3) = _wmp[2]*_cou[3][7][4)  & 
               +    faka * (_cou[0][7][3) - fakaac * _cou[0][7][4) ) 
          _cou[9][8][3) = _wmp[2]*_cou[3][8][4)  & 
               +    faka * (_cou[0][8][3) - fakaac * _cou[0][8][4) ) 
          _cou[9][9][3) = _wmp[2]*_cou[3][9][4)  & 
               +  faka * (_cou[0][9][3) - fakaac * _cou[0][9][4) ) + fakac2 * _cou[3][3][4)
          _cou[13][1][3) = _wmp[0]*_cou[4][1][4)  & 
               +   faka * (_cou[2][1][3) - fakaac * _cou[2][1][4) ) + fakac * _cou[4][0][4)
          _cou[14][1][3) = _wmp[1]*_cou[4][1][4)  & 
               +    faka * (_cou[1][1][3) - fakaac * _cou[1][1][4) ) 
          _cou[19][1][3) = _wmp[2]*_cou[4][1][4) 
          _cou[13][2][3) = _wmp[0]*_cou[4][2][4)  & 
               +    faka * (_cou[2][2][3) - fakaac * _cou[2][2][4) ) 
          _cou[14][2][3) = _wmp[1]*_cou[4][2][4)  & 
               +   faka * (_cou[1][2][3) - fakaac * _cou[1][2][4) ) + fakac * _cou[4][0][4)
          _cou[19][2][3) = _wmp[2]*_cou[4][2][4) 
          _cou[13][3][3) = _wmp[0]*_cou[4][3][4)  & 
               +    faka * (_cou[2][3][3) - fakaac * _cou[2][3][4) ) 
          _cou[14][3][3) = _wmp[1]*_cou[4][3][4)  & 
               +    faka * (_cou[1][3][3) - fakaac * _cou[1][3][4) ) 
          _cou[19][3][3) = _wmp[2]*_cou[4][3][4) +  fakac * _cou[4][0][4)
          _cou[15][1][3) = _wmp[0]*_cou[5][1][4)  & 
               +   faka * (_cou[3][1][3) - fakaac * _cou[3][1][4) ) + fakac * _cou[5][0][4)
          _cou[16][1][3) = _wmp[2]*_cou[5][1][4)  & 
               +    faka * (_cou[1][1][3) - fakaac * _cou[1][1][4) ) 
          _cou[15][2][3) = _wmp[0]*_cou[5][2][4)  & 
               +    faka * (_cou[3][2][3) - fakaac * _cou[3][2][4) ) 
          _cou[16][2][3) = _wmp[2]*_cou[5][2][4)  & 
               +    faka * (_cou[1][2][3) - fakaac * _cou[1][2][4) ) 
          _cou[15][3][3) = _wmp[0]*_cou[5][3][4)  & 
               +    faka * (_cou[3][3][3) - fakaac * _cou[3][3][4) ) 
          _cou[16][3][3) = _wmp[2]*_cou[5][3][4)  & 
               +   faka * (_cou[1][3][3) - fakaac * _cou[1][3][4) ) + fakac * _cou[5][0][4)
          _cou[17][1][3) = _wmp[1]*_cou[6][1][4)  & 
               +    faka * (_cou[3][1][3) - fakaac * _cou[3][1][4) ) 
          _cou[18][1][3) = _wmp[2]*_cou[6][1][4)  & 
               +    faka * (_cou[2][1][3) - fakaac * _cou[2][1][4) ) 
          _cou[17][2][3) = _wmp[1]*_cou[6][2][4)  & 
               +   faka * (_cou[3][2][3) - fakaac * _cou[3][2][4) ) + fakac * _cou[6][0][4)
          _cou[18][2][3) = _wmp[2]*_cou[6][2][4)  & 
               +    faka * (_cou[2][2][3) - fakaac * _cou[2][2][4) ) 
          _cou[17][3][3) = _wmp[1]*_cou[6][3][4)  & 
               +    faka * (_cou[3][3][3) - fakaac * _cou[3][3][4) ) 
          _cou[18][3][3) = _wmp[2]*_cou[6][3][4)  & 
               +   faka * (_cou[2][3][3) - fakaac * _cou[2][3][4) ) + fakac * _cou[6][0][4)
          _cou[10][1][3) = _wmp[0]*_cou[7][1][4)  & 
               +  faka2 * (_cou[1][1][3) - fakaac * _cou[1][1][4) ) + fakac * _cou[7][0][4)
          _cou[10][2][3) = _wmp[0]*_cou[7][2][4)  & 
               +  faka2 * (_cou[1][2][3) - fakaac * _cou[1][2][4) ) 
          _cou[10][3][3) = _wmp[0]*_cou[7][3][4)  & 
               +  faka2 * (_cou[1][3][3) - fakaac * _cou[1][3][4) ) 
          _cou[11][1][3) = _wmp[1]*_cou[8][1][4)  & 
               +  faka2 * (_cou[2][1][3) - fakaac * _cou[2][1][4) ) 
          _cou[11][2][3) = _wmp[1]*_cou[8][2][4)  & 
               +  faka2 * (_cou[2][2][3) - fakaac * _cou[2][2][4) ) + fakac * _cou[8][0][4)
          _cou[11][3][3) = _wmp[1]*_cou[8][3][4)  & 
               +  faka2 * (_cou[2][3][3) - fakaac * _cou[2][3][4) ) 
          _cou[12][1][3) = _wmp[2]*_cou[9][1][4)  & 
               +  faka2 * (_cou[3][1][3) - fakaac * _cou[3][1][4) ) 
          _cou[12][2][3) = _wmp[2]*_cou[9][2][4)  & 
               +  faka2 * (_cou[3][2][3) - fakaac * _cou[3][2][4) ) 
          _cou[12][3][3) = _wmp[2]*_cou[9][3][4)  & 
               +  faka2 * (_cou[3][3][3) - fakaac * _cou[3][3][4) ) + fakac * _cou[9][0][4)


          _cou[20][0][3) = _wmp[0]*_cou[10][0][4)  & 
               +    faka3 * (_cou[7][0][3) - fakaac * _cou[7][0][4) ) 
          _cou[23][0][3) = _wmp[1]*_cou[10][0][4) 
          _cou[25][0][3) = _wmp[2]*_cou[10][0][4) 
          _cou[24][0][3) = _wmp[0]*_cou[11][0][4) 
          _cou[21][0][3) = _wmp[1]*_cou[11][0][4)  & 
               +    faka3 * (_cou[8][0][3) - fakaac * _cou[8][0][4) ) 
          _cou[27][0][3) = _wmp[2]*_cou[11][0][4) 
          _cou[26][0][3) = _wmp[0]*_cou[12][0][4) 
          _cou[28][0][3) = _wmp[1]*_cou[12][0][4) 
          _cou[22][0][3) = _wmp[2]*_cou[12][0][4)  & 
               +    faka3 * (_cou[9][0][3) - fakaac * _cou[9][0][4) ) 
          _cou[31][0][3) = _wmp[1]*_cou[13][0][4)  & 
               +    faka * (_cou[7][0][3) - fakaac * _cou[7][0][4) ) 
          _cou[32][0][3) = _wmp[2]*_cou[13][0][4) 
          _cou[33][0][3) = _wmp[2]*_cou[14][0][4) 
          _cou[29][0][3) = _wmp[2]*_cou[15][0][4)  & 
               +    faka * (_cou[7][0][3) - fakaac * _cou[7][0][4) ) 
          _cou[34][0][3) = _wmp[1]*_cou[16][0][4) 
          _cou[30][0][3) = _wmp[2]*_cou[17][0][4)  & 
               +    faka * (_cou[8][0][3) - fakaac * _cou[8][0][4) ) 

          _cou[1][20][2) = _wmp[0]*_cou[0][20][3) +  fakac4 * _cou[0][10][3)
          _cou[2][20][2) = _wmp[1]*_cou[0][20][3) 
          _cou[3][20][2) = _wmp[2]*_cou[0][20][3) 
          _cou[1][21][2) = _wmp[0]*_cou[0][21][3) 
          _cou[2][21][2) = _wmp[1]*_cou[0][21][3) +  fakac4 * _cou[0][11][3)
          _cou[3][21][2) = _wmp[2]*_cou[0][21][3) 
          _cou[1][22][2) = _wmp[0]*_cou[0][22][3) 
          _cou[2][22][2) = _wmp[1]*_cou[0][22][3) 
          _cou[3][22][2) = _wmp[2]*_cou[0][22][3) +  fakac4 * _cou[0][12][3)
          _cou[1][23][2) = _wmp[0]*_cou[0][23][3) +  fakac3 * _cou[0][13][3)
          _cou[2][23][2) = _wmp[1]*_cou[0][23][3) +  fakac * _cou[0][10][3)
          _cou[3][23][2) = _wmp[2]*_cou[0][23][3) 
          _cou[1][24][2) = _wmp[0]*_cou[0][24][3) +  fakac * _cou[0][11][3)
          _cou[2][24][2) = _wmp[1]*_cou[0][24][3) +  fakac3 * _cou[0][14][3)
          _cou[3][24][2) = _wmp[2]*_cou[0][24][3) 
          _cou[1][25][2) = _wmp[0]*_cou[0][25][3) +  fakac3 * _cou[0][15][3)
          _cou[2][25][2) = _wmp[1]*_cou[0][25][3) 
          _cou[3][25][2) = _wmp[2]*_cou[0][25][3) +  fakac * _cou[0][10][3)
          _cou[1][26][2) = _wmp[0]*_cou[0][26][3) +  fakac * _cou[0][12][3)
          _cou[2][26][2) = _wmp[1]*_cou[0][26][3) 
          _cou[3][26][2) = _wmp[2]*_cou[0][26][3) +  fakac3 * _cou[0][16][3)
          _cou[1][27][2) = _wmp[0]*_cou[0][27][3) 
          _cou[2][27][2) = _wmp[1]*_cou[0][27][3) +  fakac3 * _cou[0][17][3)
          _cou[3][27][2) = _wmp[2]*_cou[0][27][3) +  fakac * _cou[0][11][3)
          _cou[1][28][2) = _wmp[0]*_cou[0][28][3) 
          _cou[2][28][2) = _wmp[1]*_cou[0][28][3) +  fakac * _cou[0][12][3)
          _cou[3][28][2) = _wmp[2]*_cou[0][28][3) +  fakac3 * _cou[0][18][3)
          _cou[1][29][2) = _wmp[0]*_cou[0][29][3) +  fakac2 * _cou[0][16][3)
          _cou[2][29][2) = _wmp[1]*_cou[0][29][3) 
          _cou[3][29][2) = _wmp[2]*_cou[0][29][3) +  fakac2 * _cou[0][15][3)
          _cou[1][30][2) = _wmp[0]*_cou[0][30][3) 
          _cou[2][30][2) = _wmp[1]*_cou[0][30][3) +  fakac2 * _cou[0][18][3)
          _cou[3][30][2) = _wmp[2]*_cou[0][30][3) +  fakac2 * _cou[0][17][3)
          _cou[1][31][2) = _wmp[0]*_cou[0][31][3) +  fakac2 * _cou[0][14][3)
          _cou[2][31][2) = _wmp[1]*_cou[0][31][3) +  fakac2 * _cou[0][13][3)
          _cou[3][31][2) = _wmp[2]*_cou[0][31][3) 
          _cou[1][32][2) = _wmp[0]*_cou[0][32][3) +  fakac2 * _cou[0][19][3)
          _cou[2][32][2) = _wmp[1]*_cou[0][32][3) +  fakac * _cou[0][15][3)
          _cou[3][32][2) = _wmp[2]*_cou[0][32][3) +  fakac * _cou[0][13][3)
          _cou[1][33][2) = _wmp[0]*_cou[0][33][3) +  fakac * _cou[0][17][3)
          _cou[2][33][2) = _wmp[1]*_cou[0][33][3) +  fakac2 * _cou[0][19][3)
          _cou[3][33][2) = _wmp[2]*_cou[0][33][3) +  fakac * _cou[0][14][3)
          _cou[1][34][2) = _wmp[0]*_cou[0][34][3) +  fakac * _cou[0][18][3)
          _cou[2][34][2) = _wmp[1]*_cou[0][34][3) +  fakac * _cou[0][16][3)
          _cou[3][34][2) = _wmp[2]*_cou[0][34][3) +  fakac2 * _cou[0][19][3)



          _cou[7][10][2) = _wmp[0]*_cou[1][10][3)  & 
               +  faka * (_cou[0][10][2) - fakaac * _cou[0][10][3) ) + fakac3 * _cou[1][7][3)
          _cou[4][10][2) = _wmp[1]*_cou[1][10][3) 
          _cou[5][10][2) = _wmp[2]*_cou[1][10][3) 
          _cou[7][11][2) = _wmp[0]*_cou[1][11][3)  & 
               +    faka * (_cou[0][11][2) - fakaac * _cou[0][11][3) ) 
          _cou[4][11][2) = _wmp[1]*_cou[1][11][3) +  fakac3 * _cou[1][8][3)
          _cou[5][11][2) = _wmp[2]*_cou[1][11][3) 
          _cou[7][12][2) = _wmp[0]*_cou[1][12][3)  & 
               +    faka * (_cou[0][12][2) - fakaac * _cou[0][12][3) ) 
          _cou[4][12][2) = _wmp[1]*_cou[1][12][3) 
          _cou[5][12][2) = _wmp[2]*_cou[1][12][3) +  fakac3 * _cou[1][9][3)
          _cou[7][13][2) = _wmp[0]*_cou[1][13][3)  & 
               +  faka * (_cou[0][13][2) - fakaac * _cou[0][13][3) ) + fakac2 * _cou[1][4][3)
          _cou[4][13][2) = _wmp[1]*_cou[1][13][3) +  fakac * _cou[1][7][3)
          _cou[5][13][2) = _wmp[2]*_cou[1][13][3) 
          _cou[7][14][2) = _wmp[0]*_cou[1][14][3)  & 
               +   faka * (_cou[0][14][2) - fakaac * _cou[0][14][3) ) + fakac * _cou[1][8][3)
          _cou[4][14][2) = _wmp[1]*_cou[1][14][3) +  fakac2 * _cou[1][4][3)
          _cou[5][14][2) = _wmp[2]*_cou[1][14][3) 
          _cou[7][15][2) = _wmp[0]*_cou[1][15][3)  & 
               +  faka * (_cou[0][15][2) - fakaac * _cou[0][15][3) ) + fakac2 * _cou[1][5][3)
          _cou[4][15][2) = _wmp[1]*_cou[1][15][3) 
          _cou[5][15][2) = _wmp[2]*_cou[1][15][3) +  fakac * _cou[1][7][3)
          _cou[7][16][2) = _wmp[0]*_cou[1][16][3)  & 
               +   faka * (_cou[0][16][2) - fakaac * _cou[0][16][3) ) + fakac * _cou[1][9][3)
          _cou[4][16][2) = _wmp[1]*_cou[1][16][3) 
          _cou[5][16][2) = _wmp[2]*_cou[1][16][3) +  fakac2 * _cou[1][5][3)
          _cou[7][17][2) = _wmp[0]*_cou[1][17][3)  & 
               +    faka * (_cou[0][17][2) - fakaac * _cou[0][17][3) ) 
          _cou[4][17][2) = _wmp[1]*_cou[1][17][3) +  fakac2 * _cou[1][6][3)
          _cou[5][17][2) = _wmp[2]*_cou[1][17][3) +  fakac * _cou[1][8][3)
          _cou[7][18][2) = _wmp[0]*_cou[1][18][3)  & 
               +    faka * (_cou[0][18][2) - fakaac * _cou[0][18][3) ) 
          _cou[4][18][2) = _wmp[1]*_cou[1][18][3) +  fakac * _cou[1][9][3)
          _cou[5][18][2) = _wmp[2]*_cou[1][18][3) +  fakac2 * _cou[1][6][3)
          _cou[7][19][2) = _wmp[0]*_cou[1][19][3)  & 
               +   faka * (_cou[0][19][2) - fakaac * _cou[0][19][3) ) + fakac * _cou[1][6][3)
          _cou[4][19][2) = _wmp[1]*_cou[1][19][3) +  fakac * _cou[1][5][3)
          _cou[5][19][2) = _wmp[2]*_cou[1][19][3) +  fakac * _cou[1][4][3)
          _cou[8][10][2) = _wmp[1]*_cou[2][10][3)  & 
               +    faka * (_cou[0][10][2) - fakaac * _cou[0][10][3) ) 
          _cou[6][10][2) = _wmp[2]*_cou[2][10][3) 
          _cou[8][11][2) = _wmp[1]*_cou[2][11][3)  & 
               +  faka * (_cou[0][11][2) - fakaac * _cou[0][11][3) ) + fakac3 * _cou[2][8][3)
          _cou[6][11][2) = _wmp[2]*_cou[2][11][3) 
          _cou[8][12][2) = _wmp[1]*_cou[2][12][3)  & 
               +    faka * (_cou[0][12][2) - fakaac * _cou[0][12][3) ) 
          _cou[6][12][2) = _wmp[2]*_cou[2][12][3) +  fakac3 * _cou[2][9][3)
          _cou[8][13][2) = _wmp[1]*_cou[2][13][3)  & 
               +   faka * (_cou[0][13][2) - fakaac * _cou[0][13][3) ) + fakac * _cou[2][7][3)
          _cou[6][13][2) = _wmp[2]*_cou[2][13][3) 
          _cou[8][14][2) = _wmp[1]*_cou[2][14][3)  & 
               +  faka * (_cou[0][14][2) - fakaac * _cou[0][14][3) ) + fakac2 * _cou[2][4][3)
          _cou[6][14][2) = _wmp[2]*_cou[2][14][3) 
          _cou[8][15][2) = _wmp[1]*_cou[2][15][3)  & 
               +    faka * (_cou[0][15][2) - fakaac * _cou[0][15][3) ) 
          _cou[6][15][2) = _wmp[2]*_cou[2][15][3) +  fakac * _cou[2][7][3)
          _cou[8][16][2) = _wmp[1]*_cou[2][16][3)  & 
               +    faka * (_cou[0][16][2) - fakaac * _cou[0][16][3) ) 
          _cou[6][16][2) = _wmp[2]*_cou[2][16][3) +  fakac2 * _cou[2][5][3)
          _cou[8][17][2) = _wmp[1]*_cou[2][17][3)  & 
               +  faka * (_cou[0][17][2) - fakaac * _cou[0][17][3) ) + fakac2 * _cou[2][6][3)
          _cou[6][17][2) = _wmp[2]*_cou[2][17][3) +  fakac * _cou[2][8][3)
          _cou[8][18][2) = _wmp[1]*_cou[2][18][3)  & 
               +   faka * (_cou[0][18][2) - fakaac * _cou[0][18][3) ) + fakac * _cou[2][9][3)
          _cou[6][18][2) = _wmp[2]*_cou[2][18][3) +  fakac2 * _cou[2][6][3)
          _cou[8][19][2) = _wmp[1]*_cou[2][19][3)  & 
               +   faka * (_cou[0][19][2) - fakaac * _cou[0][19][3) ) + fakac * _cou[2][5][3)
          _cou[6][19][2) = _wmp[2]*_cou[2][19][3) +  fakac * _cou[2][4][3)
          _cou[9][10][2) = _wmp[2]*_cou[3][10][3)  & 
               +    faka * (_cou[0][10][2) - fakaac * _cou[0][10][3) ) 
          _cou[9][11][2) = _wmp[2]*_cou[3][11][3)  & 
               +    faka * (_cou[0][11][2) - fakaac * _cou[0][11][3) ) 
          _cou[9][12][2) = _wmp[2]*_cou[3][12][3)  & 
               +  faka * (_cou[0][12][2) - fakaac * _cou[0][12][3) ) + fakac3 * _cou[3][9][3)
          _cou[9][13][2) = _wmp[2]*_cou[3][13][3)  & 
               +    faka * (_cou[0][13][2) - fakaac * _cou[0][13][3) ) 
          _cou[9][14][2) = _wmp[2]*_cou[3][14][3)  & 
               +    faka * (_cou[0][14][2) - fakaac * _cou[0][14][3) ) 
          _cou[9][15][2) = _wmp[2]*_cou[3][15][3)  & 
               +   faka * (_cou[0][15][2) - fakaac * _cou[0][15][3) ) + fakac * _cou[3][7][3)
          _cou[9][16][2) = _wmp[2]*_cou[3][16][3)  & 
               +  faka * (_cou[0][16][2) - fakaac * _cou[0][16][3) ) + fakac2 * _cou[3][5][3)
          _cou[9][17][2) = _wmp[2]*_cou[3][17][3)  & 
               +   faka * (_cou[0][17][2) - fakaac * _cou[0][17][3) ) + fakac * _cou[3][8][3)
          _cou[9][18][2) = _wmp[2]*_cou[3][18][3)  & 
               +  faka * (_cou[0][18][2) - fakaac * _cou[0][18][3) ) + fakac2 * _cou[3][6][3)
          _cou[9][19][2) = _wmp[2]*_cou[3][19][3)  & 
               +   faka * (_cou[0][19][2) - fakaac * _cou[0][19][3) ) + fakac * _cou[3][4][3)


          _cou[13][4][2) = _wmp[0]*_cou[4][4][3)  & 
               +  faka * (_cou[2][4][2) - fakaac * _cou[2][4][3) ) + fakac * _cou[4][2][3)
          _cou[14][4][2) = _wmp[1]*_cou[4][4][3)  & 
               +  faka * (_cou[1][4][2) - fakaac * _cou[1][4][3) ) + fakac * _cou[4][1][3)
          _cou[19][4][2) = _wmp[2]*_cou[4][4][3) 
          _cou[13][5][2) = _wmp[0]*_cou[4][5][3)  & 
               +  faka * (_cou[2][5][2) - fakaac * _cou[2][5][3) ) + fakac * _cou[4][3][3)
          _cou[14][5][2) = _wmp[1]*_cou[4][5][3)  & 
               +  faka * (_cou[1][5][2) - fakaac * _cou[1][5][3) ) 
          _cou[19][5][2) = _wmp[2]*_cou[4][5][3) +  fakac * _cou[4][1][3)
          _cou[13][6][2) = _wmp[0]*_cou[4][6][3)  & 
               +  faka * (_cou[2][6][2) - fakaac * _cou[2][6][3) ) 
          _cou[14][6][2) = _wmp[1]*_cou[4][6][3)  & 
               +  faka * (_cou[1][6][2) - fakaac * _cou[1][6][3) ) + fakac * _cou[4][3][3)
          _cou[19][6][2) = _wmp[2]*_cou[4][6][3) +  fakac * _cou[4][2][3)
          _cou[13][7][2) = _wmp[0]*_cou[4][7][3)  & 
               +  faka * (_cou[2][7][2) - fakaac * _cou[2][7][3) ) + fakac2 * _cou[4][1][3)
          _cou[14][7][2) = _wmp[1]*_cou[4][7][3)  & 
               +  faka * (_cou[1][7][2) - fakaac * _cou[1][7][3) ) 
          _cou[19][7][2) = _wmp[2]*_cou[4][7][3) 
          _cou[13][8][2) = _wmp[0]*_cou[4][8][3)  & 
               +  faka * (_cou[2][8][2) - fakaac * _cou[2][8][3) ) 
          _cou[14][8][2) = _wmp[1]*_cou[4][8][3)  & 
               +  faka * (_cou[1][8][2) - fakaac * _cou[1][8][3) ) + fakac2 * _cou[4][2][3)
          _cou[19][8][2) = _wmp[2]*_cou[4][8][3) 
          _cou[13][9][2) = _wmp[0]*_cou[4][9][3)  & 
               +  faka * (_cou[2][9][2) - fakaac * _cou[2][9][3) ) 
          _cou[14][9][2) = _wmp[1]*_cou[4][9][3)  & 
               +  faka * (_cou[1][9][2) - fakaac * _cou[1][9][3) ) 
          _cou[19][9][2) = _wmp[2]*_cou[4][9][3) +  fakac2 * _cou[4][3][3)
          _cou[15][4][2) = _wmp[0]*_cou[5][4][3)  & 
               +  faka * (_cou[3][4][2) - fakaac * _cou[3][4][3) ) + fakac * _cou[5][2][3)
          _cou[16][4][2) = _wmp[2]*_cou[5][4][3)  & 
               +  faka * (_cou[1][4][2) - fakaac * _cou[1][4][3) ) 
          _cou[15][5][2) = _wmp[0]*_cou[5][5][3)  & 
               +  faka * (_cou[3][5][2) - fakaac * _cou[3][5][3) ) + fakac * _cou[5][3][3)
          _cou[16][5][2) = _wmp[2]*_cou[5][5][3)  & 
               +  faka * (_cou[1][5][2) - fakaac * _cou[1][5][3) ) + fakac * _cou[5][1][3)
          _cou[15][6][2) = _wmp[0]*_cou[5][6][3)  & 
               +  faka * (_cou[3][6][2) - fakaac * _cou[3][6][3) ) 
          _cou[16][6][2) = _wmp[2]*_cou[5][6][3)  & 
               +  faka * (_cou[1][6][2) - fakaac * _cou[1][6][3) ) + fakac * _cou[5][2][3)
          _cou[15][7][2) = _wmp[0]*_cou[5][7][3)  & 
               +  faka * (_cou[3][7][2) - fakaac * _cou[3][7][3) ) + fakac2 * _cou[5][1][3)
          _cou[16][7][2) = _wmp[2]*_cou[5][7][3)  & 
               +  faka * (_cou[1][7][2) - fakaac * _cou[1][7][3) ) 
          _cou[15][8][2) = _wmp[0]*_cou[5][8][3)  & 
               +  faka * (_cou[3][8][2) - fakaac * _cou[3][8][3) ) 
          _cou[16][8][2) = _wmp[2]*_cou[5][8][3)  & 
               +  faka * (_cou[1][8][2) - fakaac * _cou[1][8][3) ) 
          _cou[15][9][2) = _wmp[0]*_cou[5][9][3)  & 
               +  faka * (_cou[3][9][2) - fakaac * _cou[3][9][3) ) 
          _cou[16][9][2) = _wmp[2]*_cou[5][9][3)  & 
               +  faka * (_cou[1][9][2) - fakaac * _cou[1][9][3) ) + fakac2 * _cou[5][3][3)
          _cou[17][4][2) = _wmp[1]*_cou[6][4][3)  & 
               +  faka * (_cou[3][4][2) - fakaac * _cou[3][4][3) ) + fakac * _cou[6][1][3)
          _cou[18][4][2) = _wmp[2]*_cou[6][4][3)  & 
               +  faka * (_cou[2][4][2) - fakaac * _cou[2][4][3) ) 
          _cou[17][5][2) = _wmp[1]*_cou[6][5][3)  & 
               +  faka * (_cou[3][5][2) - fakaac * _cou[3][5][3) ) 
          _cou[18][5][2) = _wmp[2]*_cou[6][5][3)  & 
               +  faka * (_cou[2][5][2) - fakaac * _cou[2][5][3) ) + fakac * _cou[6][1][3)
          _cou[17][6][2) = _wmp[1]*_cou[6][6][3)  & 
               +  faka * (_cou[3][6][2) - fakaac * _cou[3][6][3) ) + fakac * _cou[6][3][3)
          _cou[18][6][2) = _wmp[2]*_cou[6][6][3)  & 
               +  faka * (_cou[2][6][2) - fakaac * _cou[2][6][3) ) + fakac * _cou[6][2][3)
          _cou[17][7][2) = _wmp[1]*_cou[6][7][3)  & 
               +  faka * (_cou[3][7][2) - fakaac * _cou[3][7][3) ) 
          _cou[18][7][2) = _wmp[2]*_cou[6][7][3)  & 
               +  faka * (_cou[2][7][2) - fakaac * _cou[2][7][3) ) 
          _cou[17][8][2) = _wmp[1]*_cou[6][8][3)  & 
               +  faka * (_cou[3][8][2) - fakaac * _cou[3][8][3) ) + fakac2 * _cou[6][2][3)
          _cou[18][8][2) = _wmp[2]*_cou[6][8][3)  & 
               +  faka * (_cou[2][8][2) - fakaac * _cou[2][8][3) ) 
          _cou[17][9][2) = _wmp[1]*_cou[6][9][3)  & 
               +  faka * (_cou[3][9][2) - fakaac * _cou[3][9][3) ) 
          _cou[18][9][2) = _wmp[2]*_cou[6][9][3)  & 
               +  faka * (_cou[2][9][2) - fakaac * _cou[2][9][3) ) + fakac2 * _cou[6][3][3)
          _cou[10][4][2) = _wmp[0]*_cou[7][4][3)  & 
               +  faka2 * (_cou[1][4][2) - fakaac * _cou[1][4][3) ) + fakac * _cou[7][2][3)
          _cou[10][5][2) = _wmp[0]*_cou[7][5][3)  & 
               +  faka2 * (_cou[1][5][2) - fakaac * _cou[1][5][3) ) + fakac * _cou[7][3][3)
          _cou[10][6][2) = _wmp[0]*_cou[7][6][3)  & 
               +  faka2 * (_cou[1][6][2) - fakaac * _cou[1][6][3) ) 
          _cou[10][7][2) = _wmp[0]*_cou[7][7][3)  & 
               + faka2 * (_cou[1][7][2) - fakaac * _cou[1][7][3) ) + fakac2 * _cou[7][1][3)
          _cou[10][8][2) = _wmp[0]*_cou[7][8][3)  & 
               +  faka2 * (_cou[1][8][2) - fakaac * _cou[1][8][3) ) 
          _cou[10][9][2) = _wmp[0]*_cou[7][9][3)  & 
               +  faka2 * (_cou[1][9][2) - fakaac * _cou[1][9][3) ) 
          _cou[11][4][2) = _wmp[1]*_cou[8][4][3)  & 
               +  faka2 * (_cou[2][4][2) - fakaac * _cou[2][4][3) ) + fakac * _cou[8][1][3)
          _cou[11][5][2) = _wmp[1]*_cou[8][5][3)  & 
               +  faka2 * (_cou[2][5][2) - fakaac * _cou[2][5][3) ) 
          _cou[11][6][2) = _wmp[1]*_cou[8][6][3)  & 
               +  faka2 * (_cou[2][6][2) - fakaac * _cou[2][6][3) ) + fakac * _cou[8][3][3)
          _cou[11][7][2) = _wmp[1]*_cou[8][7][3)  & 
               +  faka2 * (_cou[2][7][2) - fakaac * _cou[2][7][3) ) 
          _cou[11][8][2) = _wmp[1]*_cou[8][8][3)  & 
               + faka2 * (_cou[2][8][2) - fakaac * _cou[2][8][3) ) + fakac2 * _cou[8][2][3)
          _cou[11][9][2) = _wmp[1]*_cou[8][9][3)  & 
               +  faka2 * (_cou[2][9][2) - fakaac * _cou[2][9][3) ) 
          _cou[12][4][2) = _wmp[2]*_cou[9][4][3)  & 
               +  faka2 * (_cou[3][4][2) - fakaac * _cou[3][4][3) ) 
          _cou[12][5][2) = _wmp[2]*_cou[9][5][3)  & 
               +  faka2 * (_cou[3][5][2) - fakaac * _cou[3][5][3) ) + fakac * _cou[9][1][3)
          _cou[12][6][2) = _wmp[2]*_cou[9][6][3)  & 
               +  faka2 * (_cou[3][6][2) - fakaac * _cou[3][6][3) ) + fakac * _cou[9][2][3)
          _cou[12][7][2) = _wmp[2]*_cou[9][7][3)  & 
               +  faka2 * (_cou[3][7][2) - fakaac * _cou[3][7][3) ) 
          _cou[12][8][2) = _wmp[2]*_cou[9][8][3)  & 
               +  faka2 * (_cou[3][8][2) - fakaac * _cou[3][8][3) ) 
          _cou[12][9][2) = _wmp[2]*_cou[9][9][3)  & 
               + faka2 * (_cou[3][9][2) - fakaac * _cou[3][9][3) ) + fakac2 * _cou[9][3][3)






          _cou[20][1][2) = _wmp[0]*_cou[10][1][3)  & 
               +  faka3 * (_cou[7][1][2) - fakaac * _cou[7][1][3) ) + fakac * _cou[10][0][3)
          _cou[23][1][2) = _wmp[1]*_cou[10][1][3) 
          _cou[25][1][2) = _wmp[2]*_cou[10][1][3) 
          _cou[20][2][2) = _wmp[0]*_cou[10][2][3)  & 
               +  faka3 * (_cou[7][2][2) - fakaac * _cou[7][2][3) ) 
          _cou[23][2][2) = _wmp[1]*_cou[10][2][3) +  fakac * _cou[10][0][3)
          _cou[25][2][2) = _wmp[2]*_cou[10][2][3) 
          _cou[20][3][2) = _wmp[0]*_cou[10][3][3)  & 
               +  faka3 * (_cou[7][3][2) - fakaac * _cou[7][3][3) ) 
          _cou[23][3][2) = _wmp[1]*_cou[10][3][3) 
          _cou[25][3][2) = _wmp[2]*_cou[10][3][3) +  fakac * _cou[10][0][3)
          _cou[24][1][2) = _wmp[0]*_cou[11][1][3) +  fakac * _cou[11][0][3)
          _cou[21][1][2) = _wmp[1]*_cou[11][1][3)  & 
               +  faka3 * (_cou[8][1][2) - fakaac * _cou[8][1][3) ) 
          _cou[27][1][2) = _wmp[2]*_cou[11][1][3) 
          _cou[24][2][2) = _wmp[0]*_cou[11][2][3) 
          _cou[21][2][2) = _wmp[1]*_cou[11][2][3)  & 
               +  faka3 * (_cou[8][2][2) - fakaac * _cou[8][2][3) ) + fakac * _cou[11][0][3)
          _cou[27][2][2) = _wmp[2]*_cou[11][2][3) 
          _cou[24][3][2) = _wmp[0]*_cou[11][3][3) 
          _cou[21][3][2) = _wmp[1]*_cou[11][3][3)  & 
               +  faka3 * (_cou[8][3][2) - fakaac * _cou[8][3][3) ) 
          _cou[27][3][2) = _wmp[2]*_cou[11][3][3) +  fakac * _cou[11][0][3)
          _cou[26][1][2) = _wmp[0]*_cou[12][1][3) +  fakac * _cou[12][0][3)
          _cou[28][1][2) = _wmp[1]*_cou[12][1][3) 
          _cou[22][1][2) = _wmp[2]*_cou[12][1][3)  & 
               +  faka3 * (_cou[9][1][2) - fakaac * _cou[9][1][3) ) 
          _cou[26][2][2) = _wmp[0]*_cou[12][2][3) 
          _cou[28][2][2) = _wmp[1]*_cou[12][2][3) +  fakac * _cou[12][0][3)
          _cou[22][2][2) = _wmp[2]*_cou[12][2][3)  & 
               +  faka3 * (_cou[9][2][2) - fakaac * _cou[9][2][3) ) 
          _cou[26][3][2) = _wmp[0]*_cou[12][3][3) 
          _cou[28][3][2) = _wmp[1]*_cou[12][3][3) 
          _cou[22][3][2) = _wmp[2]*_cou[12][3][3)  & 
               +  faka3 * (_cou[9][3][2) - fakaac * _cou[9][3][3) ) + fakac * _cou[12][0][3)
          _cou[31][1][2) = _wmp[1]*_cou[13][1][3)  & 
               +  faka * (_cou[7][1][2) - fakaac * _cou[7][1][3) ) 
          _cou[32][1][2) = _wmp[2]*_cou[13][1][3) 
          _cou[31][2][2) = _wmp[1]*_cou[13][2][3)  & 
               +  faka * (_cou[7][2][2) - fakaac * _cou[7][2][3) ) + fakac * _cou[13][0][3)
          _cou[32][2][2) = _wmp[2]*_cou[13][2][3) 
          _cou[31][3][2) = _wmp[1]*_cou[13][3][3)  & 
               +  faka * (_cou[7][3][2) - fakaac * _cou[7][3][3) ) 
          _cou[32][3][2) = _wmp[2]*_cou[13][3][3) +  fakac * _cou[13][0][3)
          _cou[33][1][2) = _wmp[2]*_cou[14][1][3) 
          _cou[33][2][2) = _wmp[2]*_cou[14][2][3) 
          _cou[33][3][2) = _wmp[2]*_cou[14][3][3) +  fakac * _cou[14][0][3)
          _cou[29][1][2) = _wmp[2]*_cou[15][1][3)  & 
               +  faka * (_cou[7][1][2) - fakaac * _cou[7][1][3) ) 
          _cou[29][2][2) = _wmp[2]*_cou[15][2][3)  & 
               +  faka * (_cou[7][2][2) - fakaac * _cou[7][2][3) ) 
          _cou[29][3][2) = _wmp[2]*_cou[15][3][3)  & 
               +  faka * (_cou[7][3][2) - fakaac * _cou[7][3][3) ) + fakac * _cou[15][0][3)
          _cou[34][1][2) = _wmp[1]*_cou[16][1][3) 
          _cou[34][2][2) = _wmp[1]*_cou[16][2][3) +  fakac * _cou[16][0][3)
          _cou[34][3][2) = _wmp[1]*_cou[16][3][3) 
          _cou[30][1][2) = _wmp[2]*_cou[17][1][3)  & 
               +  faka * (_cou[8][1][2) - fakaac * _cou[8][1][3) ) 
          _cou[30][2][2) = _wmp[2]*_cou[17][2][3)  & 
               +  faka * (_cou[8][2][2) - fakaac * _cou[8][2][3) ) 
          _cou[30][3][2) = _wmp[2]*_cou[17][3][3)  & 
               +  faka * (_cou[8][3][2) - fakaac * _cou[8][3][3) ) + fakac * _cou[17][0][3)






          _cou[7][20][1) = _wmp[0]*_cou[1][20][2)  & 
               +  faka * (_cou[0][20][1) - fakaac * _cou[0][20][2) ) + fakac4 * _cou[1][10][2)
          _cou[4][20][1) = _wmp[1]*_cou[1][20][2) 
          _cou[5][20][1) = _wmp[2]*_cou[1][20][2) 
          _cou[7][21][1) = _wmp[0]*_cou[1][21][2)  & 
               +  faka * (_cou[0][21][1) - fakaac * _cou[0][21][2) ) 
          _cou[4][21][1) = _wmp[1]*_cou[1][21][2) +  fakac4 * _cou[1][11][2)
          _cou[5][21][1) = _wmp[2]*_cou[1][21][2) 
          _cou[7][22][1) = _wmp[0]*_cou[1][22][2)  & 
               +  faka * (_cou[0][22][1) - fakaac * _cou[0][22][2) ) 
          _cou[4][22][1) = _wmp[1]*_cou[1][22][2) 
          _cou[5][22][1) = _wmp[2]*_cou[1][22][2) +  fakac4 * _cou[1][12][2)
          _cou[7][23][1) = _wmp[0]*_cou[1][23][2)  & 
               +  faka * (_cou[0][23][1) - fakaac * _cou[0][23][2) ) + fakac3 * _cou[1][13][2)
          _cou[4][23][1) = _wmp[1]*_cou[1][23][2) +  fakac * _cou[1][10][2)
          _cou[5][23][1) = _wmp[2]*_cou[1][23][2) 
          _cou[7][24][1) = _wmp[0]*_cou[1][24][2)  & 
               +  faka * (_cou[0][24][1) - fakaac * _cou[0][24][2) ) + fakac * _cou[1][11][2)
          _cou[4][24][1) = _wmp[1]*_cou[1][24][2) +  fakac3 * _cou[1][14][2)
          _cou[5][24][1) = _wmp[2]*_cou[1][24][2) 
          _cou[7][25][1) = _wmp[0]*_cou[1][25][2)  & 
               +  faka * (_cou[0][25][1) - fakaac * _cou[0][25][2) ) + fakac3 * _cou[1][15][2)
          _cou[4][25][1) = _wmp[1]*_cou[1][25][2) 
          _cou[5][25][1) = _wmp[2]*_cou[1][25][2) +  fakac * _cou[1][10][2)
          _cou[7][26][1) = _wmp[0]*_cou[1][26][2)  & 
               +  faka * (_cou[0][26][1) - fakaac * _cou[0][26][2) ) + fakac * _cou[1][12][2)
          _cou[4][26][1) = _wmp[1]*_cou[1][26][2) 
          _cou[5][26][1) = _wmp[2]*_cou[1][26][2) +  fakac3 * _cou[1][16][2)
          _cou[7][27][1) = _wmp[0]*_cou[1][27][2)  & 
               +  faka * (_cou[0][27][1) - fakaac * _cou[0][27][2) ) 
          _cou[4][27][1) = _wmp[1]*_cou[1][27][2) +  fakac3 * _cou[1][17][2)
          _cou[5][27][1) = _wmp[2]*_cou[1][27][2) +  fakac * _cou[1][11][2)
          _cou[7][28][1) = _wmp[0]*_cou[1][28][2)  & 
               +  faka * (_cou[0][28][1) - fakaac * _cou[0][28][2) ) 
          _cou[4][28][1) = _wmp[1]*_cou[1][28][2) +  fakac * _cou[1][12][2)
          _cou[5][28][1) = _wmp[2]*_cou[1][28][2) +  fakac3 * _cou[1][18][2)
          _cou[7][29][1) = _wmp[0]*_cou[1][29][2)  & 
               +  faka * (_cou[0][29][1) - fakaac * _cou[0][29][2) ) + fakac2 * _cou[1][16][2)
          _cou[4][29][1) = _wmp[1]*_cou[1][29][2) 
          _cou[5][29][1) = _wmp[2]*_cou[1][29][2) +  fakac2 * _cou[1][15][2)
          _cou[7][30][1) = _wmp[0]*_cou[1][30][2)  & 
               +  faka * (_cou[0][30][1) - fakaac * _cou[0][30][2) ) 
          _cou[4][30][1) = _wmp[1]*_cou[1][30][2) +  fakac2 * _cou[1][18][2)
          _cou[5][30][1) = _wmp[2]*_cou[1][30][2) +  fakac2 * _cou[1][17][2)
          _cou[7][31][1) = _wmp[0]*_cou[1][31][2)  & 
               +  faka * (_cou[0][31][1) - fakaac * _cou[0][31][2) ) + fakac2 * _cou[1][14][2)
          _cou[4][31][1) = _wmp[1]*_cou[1][31][2) +  fakac2 * _cou[1][13][2)
          _cou[5][31][1) = _wmp[2]*_cou[1][31][2) 
          _cou[7][32][1) = _wmp[0]*_cou[1][32][2)  & 
               +  faka * (_cou[0][32][1) - fakaac * _cou[0][32][2) ) + fakac2 * _cou[1][19][2)
          _cou[4][32][1) = _wmp[1]*_cou[1][32][2) +  fakac * _cou[1][15][2)
          _cou[5][32][1) = _wmp[2]*_cou[1][32][2) +  fakac * _cou[1][13][2)
          _cou[7][33][1) = _wmp[0]*_cou[1][33][2)  & 
               +  faka * (_cou[0][33][1) - fakaac * _cou[0][33][2) ) + fakac * _cou[1][17][2)
          _cou[4][33][1) = _wmp[1]*_cou[1][33][2) +  fakac2 * _cou[1][19][2)
          _cou[5][33][1) = _wmp[2]*_cou[1][33][2) +  fakac * _cou[1][14][2)
          _cou[7][34][1) = _wmp[0]*_cou[1][34][2)  & 
               +  faka * (_cou[0][34][1) - fakaac * _cou[0][34][2) ) + fakac * _cou[1][18][2)
          _cou[4][34][1) = _wmp[1]*_cou[1][34][2) +  fakac * _cou[1][16][2)
          _cou[5][34][1) = _wmp[2]*_cou[1][34][2) +  fakac2 * _cou[1][19][2)
          _cou[8][20][1) = _wmp[1]*_cou[2][20][2)  & 
               +  faka * (_cou[0][20][1) - fakaac * _cou[0][20][2) ) 
          _cou[6][20][1) = _wmp[2]*_cou[2][20][2) 
          _cou[8][21][1) = _wmp[1]*_cou[2][21][2)  & 
               +  faka * (_cou[0][21][1) - fakaac * _cou[0][21][2) ) + fakac4 * _cou[2][11][2)
          _cou[6][21][1) = _wmp[2]*_cou[2][21][2) 
          _cou[8][22][1) = _wmp[1]*_cou[2][22][2)  & 
               +  faka * (_cou[0][22][1) - fakaac * _cou[0][22][2) ) 
          _cou[6][22][1) = _wmp[2]*_cou[2][22][2) +  fakac4 * _cou[2][12][2)
          _cou[8][23][1) = _wmp[1]*_cou[2][23][2)  & 
               +  faka * (_cou[0][23][1) - fakaac * _cou[0][23][2) ) + fakac * _cou[2][10][2)
          _cou[6][23][1) = _wmp[2]*_cou[2][23][2) 
          _cou[8][24][1) = _wmp[1]*_cou[2][24][2)  & 
               +  faka * (_cou[0][24][1) - fakaac * _cou[0][24][2) ) + fakac3 * _cou[2][14][2)
          _cou[6][24][1) = _wmp[2]*_cou[2][24][2) 
          _cou[8][25][1) = _wmp[1]*_cou[2][25][2)  & 
               +  faka * (_cou[0][25][1) - fakaac * _cou[0][25][2) ) 
          _cou[6][25][1) = _wmp[2]*_cou[2][25][2) +  fakac * _cou[2][10][2)
          _cou[8][26][1) = _wmp[1]*_cou[2][26][2)  & 
               +  faka * (_cou[0][26][1) - fakaac * _cou[0][26][2) ) 
          _cou[6][26][1) = _wmp[2]*_cou[2][26][2) +  fakac3 * _cou[2][16][2)
          _cou[8][27][1) = _wmp[1]*_cou[2][27][2)  & 
               +  faka * (_cou[0][27][1) - fakaac * _cou[0][27][2) ) + fakac3 * _cou[2][17][2)
          _cou[6][27][1) = _wmp[2]*_cou[2][27][2) +  fakac * _cou[2][11][2)
          _cou[8][28][1) = _wmp[1]*_cou[2][28][2)  & 
               +  faka * (_cou[0][28][1) - fakaac * _cou[0][28][2) ) + fakac * _cou[2][12][2)
          _cou[6][28][1) = _wmp[2]*_cou[2][28][2) +  fakac3 * _cou[2][18][2)
          _cou[8][29][1) = _wmp[1]*_cou[2][29][2)  & 
               +  faka * (_cou[0][29][1) - fakaac * _cou[0][29][2) ) 
          _cou[6][29][1) = _wmp[2]*_cou[2][29][2) +  fakac2 * _cou[2][15][2)
          _cou[8][30][1) = _wmp[1]*_cou[2][30][2)  & 
               +  faka * (_cou[0][30][1) - fakaac * _cou[0][30][2) ) + fakac2 * _cou[2][18][2)
          _cou[6][30][1) = _wmp[2]*_cou[2][30][2) +  fakac2 * _cou[2][17][2)
          _cou[8][31][1) = _wmp[1]*_cou[2][31][2)  & 
               +  faka * (_cou[0][31][1) - fakaac * _cou[0][31][2) ) + fakac2 * _cou[2][13][2)
          _cou[6][31][1) = _wmp[2]*_cou[2][31][2) 
          _cou[8][32][1) = _wmp[1]*_cou[2][32][2)  & 
               +  faka * (_cou[0][32][1) - fakaac * _cou[0][32][2) ) + fakac * _cou[2][15][2)
          _cou[6][32][1) = _wmp[2]*_cou[2][32][2) +  fakac * _cou[2][13][2)
          _cou[8][33][1) = _wmp[1]*_cou[2][33][2)  & 
               +  faka * (_cou[0][33][1) - fakaac * _cou[0][33][2) ) + fakac2 * _cou[2][19][2)
          _cou[6][33][1) = _wmp[2]*_cou[2][33][2) +  fakac * _cou[2][14][2)
          _cou[8][34][1) = _wmp[1]*_cou[2][34][2)  & 
               +  faka * (_cou[0][34][1) - fakaac * _cou[0][34][2) ) + fakac * _cou[2][16][2)
          _cou[6][34][1) = _wmp[2]*_cou[2][34][2) +  fakac2 * _cou[2][19][2)
          _cou[9][20][1) = _wmp[2]*_cou[3][20][2)  & 
               +  faka * (_cou[0][20][1) - fakaac * _cou[0][20][2) ) 
          _cou[9][21][1) = _wmp[2]*_cou[3][21][2)  & 
               +  faka * (_cou[0][21][1) - fakaac * _cou[0][21][2) ) 
          _cou[9][22][1) = _wmp[2]*_cou[3][22][2)  & 
               +  faka * (_cou[0][22][1) - fakaac * _cou[0][22][2) ) + fakac4 * _cou[3][12][2)
          _cou[9][23][1) = _wmp[2]*_cou[3][23][2)  & 
               +  faka * (_cou[0][23][1) - fakaac * _cou[0][23][2) ) 
          _cou[9][24][1) = _wmp[2]*_cou[3][24][2)  & 
               +  faka * (_cou[0][24][1) - fakaac * _cou[0][24][2) ) 
          _cou[9][25][1) = _wmp[2]*_cou[3][25][2)  & 
               +  faka * (_cou[0][25][1) - fakaac * _cou[0][25][2) ) + fakac * _cou[3][10][2)
          _cou[9][26][1) = _wmp[2]*_cou[3][26][2)  & 
               +  faka * (_cou[0][26][1) - fakaac * _cou[0][26][2) ) + fakac3 * _cou[3][16][2)
          _cou[9][27][1) = _wmp[2]*_cou[3][27][2)  & 
               +  faka * (_cou[0][27][1) - fakaac * _cou[0][27][2) ) + fakac * _cou[3][11][2)
          _cou[9][28][1) = _wmp[2]*_cou[3][28][2)  & 
               +  faka * (_cou[0][28][1) - fakaac * _cou[0][28][2) ) + fakac3 * _cou[3][18][2)
          _cou[9][29][1) = _wmp[2]*_cou[3][29][2)  & 
               +  faka * (_cou[0][29][1) - fakaac * _cou[0][29][2) ) + fakac2 * _cou[3][15][2)
          _cou[9][30][1) = _wmp[2]*_cou[3][30][2)  & 
               +  faka * (_cou[0][30][1) - fakaac * _cou[0][30][2) ) + fakac2 * _cou[3][17][2)
          _cou[9][31][1) = _wmp[2]*_cou[3][31][2)  & 
               +  faka * (_cou[0][31][1) - fakaac * _cou[0][31][2) ) 
          _cou[9][32][1) = _wmp[2]*_cou[3][32][2)  & 
               +  faka * (_cou[0][32][1) - fakaac * _cou[0][32][2) ) + fakac * _cou[3][13][2)
          _cou[9][33][1) = _wmp[2]*_cou[3][33][2)  & 
               +  faka * (_cou[0][33][1) - fakaac * _cou[0][33][2) ) + fakac * _cou[3][14][2)
          _cou[9][34][1) = _wmp[2]*_cou[3][34][2)  & 
               +  faka * (_cou[0][34][1) - fakaac * _cou[0][34][2) ) + fakac2 * _cou[3][19][2)






          _cou[13][10][1) = _wmp[0]*_cou[4][10][2)  & 
               +  faka * (_cou[2][10][1) - fakaac * _cou[2][10][2) ) + fakac3 * _cou[4][7][2)
          _cou[14][10][1) = _wmp[1]*_cou[4][10][2)  & 
               +  faka * (_cou[1][10][1) - fakaac * _cou[1][10][2) ) 
          _cou[19][10][1) = _wmp[2]*_cou[4][10][2) 
          _cou[13][11][1) = _wmp[0]*_cou[4][11][2)  & 
               +  faka * (_cou[2][11][1) - fakaac * _cou[2][11][2) ) 
          _cou[14][11][1) = _wmp[1]*_cou[4][11][2)  & 
               +  faka * (_cou[1][11][1) - fakaac * _cou[1][11][2) ) + fakac3 * _cou[4][8][2)
          _cou[19][11][1) = _wmp[2]*_cou[4][11][2) 
          _cou[13][12][1) = _wmp[0]*_cou[4][12][2)  & 
               +  faka * (_cou[2][12][1) - fakaac * _cou[2][12][2) ) 
          _cou[14][12][1) = _wmp[1]*_cou[4][12][2)  & 
               +  faka * (_cou[1][12][1) - fakaac * _cou[1][12][2) ) 
          _cou[19][12][1) = _wmp[2]*_cou[4][12][2) +  fakac3 * _cou[4][9][2)
          _cou[13][13][1) = _wmp[0]*_cou[4][13][2)  & 
               +  faka * (_cou[2][13][1) - fakaac * _cou[2][13][2) ) + fakac2 * _cou[4][4][2)
          _cou[14][13][1) = _wmp[1]*_cou[4][13][2)  & 
               +  faka * (_cou[1][13][1) - fakaac * _cou[1][13][2) ) + fakac * _cou[4][7][2)
          _cou[19][13][1) = _wmp[2]*_cou[4][13][2) 
          _cou[13][14][1) = _wmp[0]*_cou[4][14][2)  & 
               +  faka * (_cou[2][14][1) - fakaac * _cou[2][14][2) ) + fakac * _cou[4][8][2)
          _cou[14][14][1) = _wmp[1]*_cou[4][14][2)  & 
               +  faka * (_cou[1][14][1) - fakaac * _cou[1][14][2) ) + fakac2 * _cou[4][4][2)
          _cou[19][14][1) = _wmp[2]*_cou[4][14][2) 
          _cou[13][15][1) = _wmp[0]*_cou[4][15][2)  & 
               +  faka * (_cou[2][15][1) - fakaac * _cou[2][15][2) ) + fakac2 * _cou[4][5][2)
          _cou[14][15][1) = _wmp[1]*_cou[4][15][2)  & 
               +  faka * (_cou[1][15][1) - fakaac * _cou[1][15][2) ) 
          _cou[19][15][1) = _wmp[2]*_cou[4][15][2) +  fakac * _cou[4][7][2)
          _cou[13][16][1) = _wmp[0]*_cou[4][16][2)  & 
               +  faka * (_cou[2][16][1) - fakaac * _cou[2][16][2) ) + fakac * _cou[4][9][2)
          _cou[14][16][1) = _wmp[1]*_cou[4][16][2)  & 
               +  faka * (_cou[1][16][1) - fakaac * _cou[1][16][2) ) 
          _cou[19][16][1) = _wmp[2]*_cou[4][16][2) +  fakac2 * _cou[4][5][2)
          _cou[13][17][1) = _wmp[0]*_cou[4][17][2)  & 
               +  faka * (_cou[2][17][1) - fakaac * _cou[2][17][2) ) 
          _cou[14][17][1) = _wmp[1]*_cou[4][17][2)  & 
               +  faka * (_cou[1][17][1) - fakaac * _cou[1][17][2) ) + fakac2 * _cou[4][6][2)
          _cou[19][17][1) = _wmp[2]*_cou[4][17][2) +  fakac * _cou[4][8][2)
          _cou[13][18][1) = _wmp[0]*_cou[4][18][2)  & 
               +  faka * (_cou[2][18][1) - fakaac * _cou[2][18][2) ) 
          _cou[14][18][1) = _wmp[1]*_cou[4][18][2)  & 
               +  faka * (_cou[1][18][1) - fakaac * _cou[1][18][2) ) + fakac * _cou[4][9][2)
          _cou[19][18][1) = _wmp[2]*_cou[4][18][2) +  fakac2 * _cou[4][6][2)
          _cou[13][19][1) = _wmp[0]*_cou[4][19][2)  & 
               +  faka * (_cou[2][19][1) - fakaac * _cou[2][19][2) ) + fakac * _cou[4][6][2)
          _cou[14][19][1) = _wmp[1]*_cou[4][19][2)  & 
               +  faka * (_cou[1][19][1) - fakaac * _cou[1][19][2) ) + fakac * _cou[4][5][2)
          _cou[19][19][1) = _wmp[2]*_cou[4][19][2) +  fakac * _cou[4][4][2)
          _cou[15][10][1) = _wmp[0]*_cou[5][10][2)  & 
               +  faka * (_cou[3][10][1) - fakaac * _cou[3][10][2) ) + fakac3 * _cou[5][7][2)
          _cou[16][10][1) = _wmp[2]*_cou[5][10][2)  & 
               +  faka * (_cou[1][10][1) - fakaac * _cou[1][10][2) ) 
          _cou[15][11][1) = _wmp[0]*_cou[5][11][2)  & 
               +  faka * (_cou[3][11][1) - fakaac * _cou[3][11][2) ) 
          _cou[16][11][1) = _wmp[2]*_cou[5][11][2)  & 
               +  faka * (_cou[1][11][1) - fakaac * _cou[1][11][2) ) 
          _cou[15][12][1) = _wmp[0]*_cou[5][12][2)  & 
               +  faka * (_cou[3][12][1) - fakaac * _cou[3][12][2) ) 
          _cou[16][12][1) = _wmp[2]*_cou[5][12][2)  & 
               +  faka * (_cou[1][12][1) - fakaac * _cou[1][12][2) ) + fakac3 * _cou[5][9][2)
          _cou[15][13][1) = _wmp[0]*_cou[5][13][2)  & 
               +  faka * (_cou[3][13][1) - fakaac * _cou[3][13][2) ) + fakac2 * _cou[5][4][2)
          _cou[16][13][1) = _wmp[2]*_cou[5][13][2)  & 
               +  faka * (_cou[1][13][1) - fakaac * _cou[1][13][2) ) 
          _cou[15][14][1) = _wmp[0]*_cou[5][14][2)  & 
               +  faka * (_cou[3][14][1) - fakaac * _cou[3][14][2) ) + fakac * _cou[5][8][2)
          _cou[16][14][1) = _wmp[2]*_cou[5][14][2)  & 
               +  faka * (_cou[1][14][1) - fakaac * _cou[1][14][2) ) 
          _cou[15][15][1) = _wmp[0]*_cou[5][15][2)  & 
               +  faka * (_cou[3][15][1) - fakaac * _cou[3][15][2) ) + fakac2 * _cou[5][5][2)
          _cou[16][15][1) = _wmp[2]*_cou[5][15][2)  & 
               +  faka * (_cou[1][15][1) - fakaac * _cou[1][15][2) ) + fakac * _cou[5][7][2)
          _cou[15][16][1) = _wmp[0]*_cou[5][16][2)  & 
               +  faka * (_cou[3][16][1) - fakaac * _cou[3][16][2) ) + fakac * _cou[5][9][2)
          _cou[16][16][1) = _wmp[2]*_cou[5][16][2)  & 
               +  faka * (_cou[1][16][1) - fakaac * _cou[1][16][2) ) + fakac2 * _cou[5][5][2)
          _cou[15][17][1) = _wmp[0]*_cou[5][17][2)  & 
               +  faka * (_cou[3][17][1) - fakaac * _cou[3][17][2) ) 
          _cou[16][17][1) = _wmp[2]*_cou[5][17][2)  & 
               +  faka * (_cou[1][17][1) - fakaac * _cou[1][17][2) ) + fakac * _cou[5][8][2)
          _cou[15][18][1) = _wmp[0]*_cou[5][18][2)  & 
               +  faka * (_cou[3][18][1) - fakaac * _cou[3][18][2) ) 
          _cou[16][18][1) = _wmp[2]*_cou[5][18][2)  & 
               +  faka * (_cou[1][18][1) - fakaac * _cou[1][18][2) ) + fakac2 * _cou[5][6][2)
          _cou[15][19][1) = _wmp[0]*_cou[5][19][2)  & 
               +  faka * (_cou[3][19][1) - fakaac * _cou[3][19][2) ) + fakac * _cou[5][6][2)
          _cou[16][19][1) = _wmp[2]*_cou[5][19][2)  & 
               +  faka * (_cou[1][19][1) - fakaac * _cou[1][19][2) ) + fakac * _cou[5][4][2)
          _cou[17][10][1) = _wmp[1]*_cou[6][10][2)  & 
               +  faka * (_cou[3][10][1) - fakaac * _cou[3][10][2) ) 
          _cou[18][10][1) = _wmp[2]*_cou[6][10][2)  & 
               +  faka * (_cou[2][10][1) - fakaac * _cou[2][10][2) ) 
          _cou[17][11][1) = _wmp[1]*_cou[6][11][2)  & 
               +  faka * (_cou[3][11][1) - fakaac * _cou[3][11][2) ) + fakac3 * _cou[6][8][2)
          _cou[18][11][1) = _wmp[2]*_cou[6][11][2)  & 
               +  faka * (_cou[2][11][1) - fakaac * _cou[2][11][2) ) 
          _cou[17][12][1) = _wmp[1]*_cou[6][12][2)  & 
               +  faka * (_cou[3][12][1) - fakaac * _cou[3][12][2) ) 
          _cou[18][12][1) = _wmp[2]*_cou[6][12][2)  & 
               +  faka * (_cou[2][12][1) - fakaac * _cou[2][12][2) ) + fakac3 * _cou[6][9][2)
          _cou[17][13][1) = _wmp[1]*_cou[6][13][2)  & 
               +  faka * (_cou[3][13][1) - fakaac * _cou[3][13][2) ) + fakac * _cou[6][7][2)
          _cou[18][13][1) = _wmp[2]*_cou[6][13][2)  & 
               +  faka * (_cou[2][13][1) - fakaac * _cou[2][13][2) ) 
          _cou[17][14][1) = _wmp[1]*_cou[6][14][2)  & 
               +  faka * (_cou[3][14][1) - fakaac * _cou[3][14][2) ) + fakac2 * _cou[6][4][2)
          _cou[18][14][1) = _wmp[2]*_cou[6][14][2)  & 
               +  faka * (_cou[2][14][1) - fakaac * _cou[2][14][2) ) 
          _cou[17][15][1) = _wmp[1]*_cou[6][15][2)  & 
               +  faka * (_cou[3][15][1) - fakaac * _cou[3][15][2) ) 
          _cou[18][15][1) = _wmp[2]*_cou[6][15][2)  & 
               +  faka * (_cou[2][15][1) - fakaac * _cou[2][15][2) ) + fakac * _cou[6][7][2)
          _cou[17][16][1) = _wmp[1]*_cou[6][16][2)  & 
               +  faka * (_cou[3][16][1) - fakaac * _cou[3][16][2) ) 
          _cou[18][16][1) = _wmp[2]*_cou[6][16][2)  & 
               +  faka * (_cou[2][16][1) - fakaac * _cou[2][16][2) ) + fakac2 * _cou[6][5][2)
          _cou[17][17][1) = _wmp[1]*_cou[6][17][2)  & 
               +  faka * (_cou[3][17][1) - fakaac * _cou[3][17][2) ) + fakac2 * _cou[6][6][2)
          _cou[18][17][1) = _wmp[2]*_cou[6][17][2)  & 
               +  faka * (_cou[2][17][1) - fakaac * _cou[2][17][2) ) + fakac * _cou[6][8][2)
          _cou[17][18][1) = _wmp[1]*_cou[6][18][2)  & 
               +  faka * (_cou[3][18][1) - fakaac * _cou[3][18][2) ) + fakac * _cou[6][9][2)
          _cou[18][18][1) = _wmp[2]*_cou[6][18][2)  & 
               +  faka * (_cou[2][18][1) - fakaac * _cou[2][18][2) ) + fakac2 * _cou[6][6][2)
          _cou[17][19][1) = _wmp[1]*_cou[6][19][2)  & 
               +  faka * (_cou[3][19][1) - fakaac * _cou[3][19][2) ) + fakac * _cou[6][5][2)
          _cou[18][19][1) = _wmp[2]*_cou[6][19][2)  & 
               +  faka * (_cou[2][19][1) - fakaac * _cou[2][19][2) ) + fakac * _cou[6][4][2)
          _cou[10][10][1) = _wmp[0]*_cou[7][10][2)  & 
               + faka2 * (_cou[1][10][1) - fakaac * _cou[1][10][2) ) + fakac3 * _cou[7][7][2)
          _cou[10][11][1) = _wmp[0]*_cou[7][11][2)  & 
               +  faka2 * (_cou[1][11][1) - fakaac * _cou[1][11][2) ) 
          _cou[10][12][1) = _wmp[0]*_cou[7][12][2)  & 
               +  faka2 * (_cou[1][12][1) - fakaac * _cou[1][12][2) ) 
          _cou[10][13][1) = _wmp[0]*_cou[7][13][2)  & 
               + faka2 * (_cou[1][13][1) - fakaac * _cou[1][13][2) ) + fakac2 * _cou[7][4][2)
          _cou[10][14][1) = _wmp[0]*_cou[7][14][2)  & 
               +  faka2 * (_cou[1][14][1) - fakaac * _cou[1][14][2) ) + fakac * _cou[7][8][2)
          _cou[10][15][1) = _wmp[0]*_cou[7][15][2)  & 
               + faka2 * (_cou[1][15][1) - fakaac * _cou[1][15][2) ) + fakac2 * _cou[7][5][2)
          _cou[10][16][1) = _wmp[0]*_cou[7][16][2)  & 
               +  faka2 * (_cou[1][16][1) - fakaac * _cou[1][16][2) ) + fakac * _cou[7][9][2)
          _cou[10][17][1) = _wmp[0]*_cou[7][17][2)  & 
               +  faka2 * (_cou[1][17][1) - fakaac * _cou[1][17][2) ) 
          _cou[10][18][1) = _wmp[0]*_cou[7][18][2)  & 
               +  faka2 * (_cou[1][18][1) - fakaac * _cou[1][18][2) ) 
          _cou[10][19][1) = _wmp[0]*_cou[7][19][2)  & 
               +  faka2 * (_cou[1][19][1) - fakaac * _cou[1][19][2) ) + fakac * _cou[7][6][2)
          _cou[11][10][1) = _wmp[1]*_cou[8][10][2)  & 
               +  faka2 * (_cou[2][10][1) - fakaac * _cou[2][10][2) ) 
          _cou[11][11][1) = _wmp[1]*_cou[8][11][2)  & 
               + faka2 * (_cou[2][11][1) - fakaac * _cou[2][11][2) ) + fakac3 * _cou[8][8][2)
          _cou[11][12][1) = _wmp[1]*_cou[8][12][2)  & 
               +  faka2 * (_cou[2][12][1) - fakaac * _cou[2][12][2) ) 
          _cou[11][13][1) = _wmp[1]*_cou[8][13][2)  & 
               +  faka2 * (_cou[2][13][1) - fakaac * _cou[2][13][2) ) + fakac * _cou[8][7][2)
          _cou[11][14][1) = _wmp[1]*_cou[8][14][2)  & 
               + faka2 * (_cou[2][14][1) - fakaac * _cou[2][14][2) ) + fakac2 * _cou[8][4][2)
          _cou[11][15][1) = _wmp[1]*_cou[8][15][2)  & 
               +  faka2 * (_cou[2][15][1) - fakaac * _cou[2][15][2) ) 
          _cou[11][16][1) = _wmp[1]*_cou[8][16][2)  & 
               +  faka2 * (_cou[2][16][1) - fakaac * _cou[2][16][2) ) 
          _cou[11][17][1) = _wmp[1]*_cou[8][17][2)  & 
               + faka2 * (_cou[2][17][1) - fakaac * _cou[2][17][2) ) + fakac2 * _cou[8][6][2)
          _cou[11][18][1) = _wmp[1]*_cou[8][18][2)  & 
               +  faka2 * (_cou[2][18][1) - fakaac * _cou[2][18][2) ) + fakac * _cou[8][9][2)
          _cou[11][19][1) = _wmp[1]*_cou[8][19][2)  & 
               +  faka2 * (_cou[2][19][1) - fakaac * _cou[2][19][2) ) + fakac * _cou[8][5][2)
          _cou[12][10][1) = _wmp[2]*_cou[9][10][2)  & 
               +  faka2 * (_cou[3][10][1) - fakaac * _cou[3][10][2) ) 
          _cou[12][11][1) = _wmp[2]*_cou[9][11][2)  & 
               +  faka2 * (_cou[3][11][1) - fakaac * _cou[3][11][2) ) 
          _cou[12][12][1) = _wmp[2]*_cou[9][12][2)  & 
               + faka2 * (_cou[3][12][1) - fakaac * _cou[3][12][2) ) + fakac3 * _cou[9][9][2)
          _cou[12][13][1) = _wmp[2]*_cou[9][13][2)  & 
               +  faka2 * (_cou[3][13][1) - fakaac * _cou[3][13][2) ) 
          _cou[12][14][1) = _wmp[2]*_cou[9][14][2)  & 
               +  faka2 * (_cou[3][14][1) - fakaac * _cou[3][14][2) ) 
          _cou[12][15][1) = _wmp[2]*_cou[9][15][2)  & 
               +  faka2 * (_cou[3][15][1) - fakaac * _cou[3][15][2) ) + fakac * _cou[9][7][2)
          _cou[12][16][1) = _wmp[2]*_cou[9][16][2)  & 
               + faka2 * (_cou[3][16][1) - fakaac * _cou[3][16][2) ) + fakac2 * _cou[9][5][2)
          _cou[12][17][1) = _wmp[2]*_cou[9][17][2)  & 
               +  faka2 * (_cou[3][17][1) - fakaac * _cou[3][17][2) ) + fakac * _cou[9][8][2)
          _cou[12][18][1) = _wmp[2]*_cou[9][18][2)  & 
               + faka2 * (_cou[3][18][1) - fakaac * _cou[3][18][2) ) + fakac2 * _cou[9][6][2)
          _cou[12][19][1) = _wmp[2]*_cou[9][19][2)  & 
               +  faka2 * (_cou[3][19][1) - fakaac * _cou[3][19][2) ) + fakac * _cou[9][4][2)





          _cou[20][4][1) = _wmp[0]*_cou[10][4][2)  & 
               +  faka3 * (_cou[7][4][1) - fakaac * _cou[7][4][2) ) + fakac * _cou[10][2][2)
          _cou[23][4][1) = _wmp[1]*_cou[10][4][2) +  fakac * _cou[10][1][2)
          _cou[25][4][1) = _wmp[2]*_cou[10][4][2) 
          _cou[20][5][1) = _wmp[0]*_cou[10][5][2)  & 
               +  faka3 * (_cou[7][5][1) - fakaac * _cou[7][5][2) ) + fakac * _cou[10][3][2)
          _cou[23][5][1) = _wmp[1]*_cou[10][5][2) 
          _cou[25][5][1) = _wmp[2]*_cou[10][5][2) +  fakac * _cou[10][1][2)
          _cou[20][6][1) = _wmp[0]*_cou[10][6][2)  & 
               +  faka3 * (_cou[7][6][1) - fakaac * _cou[7][6][2) ) 
          _cou[23][6][1) = _wmp[1]*_cou[10][6][2) +  fakac * _cou[10][3][2)
          _cou[25][6][1) = _wmp[2]*_cou[10][6][2) +  fakac * _cou[10][2][2)
          _cou[20][7][1) = _wmp[0]*_cou[10][7][2)  & 
               + faka3 * (_cou[7][7][1) - fakaac * _cou[7][7][2) ) + fakac2 * _cou[10][1][2)
          _cou[23][7][1) = _wmp[1]*_cou[10][7][2) 
          _cou[25][7][1) = _wmp[2]*_cou[10][7][2) 
          _cou[20][8][1) = _wmp[0]*_cou[10][8][2)  & 
               +  faka3 * (_cou[7][8][1) - fakaac * _cou[7][8][2) ) 
          _cou[23][8][1) = _wmp[1]*_cou[10][8][2) +  fakac2 * _cou[10][2][2)
          _cou[25][8][1) = _wmp[2]*_cou[10][8][2) 
          _cou[20][9][1) = _wmp[0]*_cou[10][9][2)  & 
               +  faka3 * (_cou[7][9][1) - fakaac * _cou[7][9][2) ) 
          _cou[23][9][1) = _wmp[1]*_cou[10][9][2) 
          _cou[25][9][1) = _wmp[2]*_cou[10][9][2) +  fakac2 * _cou[10][3][2)
          _cou[24][4][1) = _wmp[0]*_cou[11][4][2) +  fakac * _cou[11][2][2)
          _cou[21][4][1) = _wmp[1]*_cou[11][4][2)  & 
               +  faka3 * (_cou[8][4][1) - fakaac * _cou[8][4][2) ) + fakac * _cou[11][1][2)
          _cou[27][4][1) = _wmp[2]*_cou[11][4][2) 
          _cou[24][5][1) = _wmp[0]*_cou[11][5][2) +  fakac * _cou[11][3][2)
          _cou[21][5][1) = _wmp[1]*_cou[11][5][2)  & 
               +  faka3 * (_cou[8][5][1) - fakaac * _cou[8][5][2) ) 
          _cou[27][5][1) = _wmp[2]*_cou[11][5][2) +  fakac * _cou[11][1][2)
          _cou[24][6][1) = _wmp[0]*_cou[11][6][2) 
          _cou[21][6][1) = _wmp[1]*_cou[11][6][2)  & 
               +  faka3 * (_cou[8][6][1) - fakaac * _cou[8][6][2) ) + fakac * _cou[11][3][2)
          _cou[27][6][1) = _wmp[2]*_cou[11][6][2) +  fakac * _cou[11][2][2)
          _cou[24][7][1) = _wmp[0]*_cou[11][7][2) +  fakac2 * _cou[11][1][2)
          _cou[21][7][1) = _wmp[1]*_cou[11][7][2)  & 
               +  faka3 * (_cou[8][7][1) - fakaac * _cou[8][7][2) ) 
          _cou[27][7][1) = _wmp[2]*_cou[11][7][2) 
          _cou[24][8][1) = _wmp[0]*_cou[11][8][2) 
          _cou[21][8][1) = _wmp[1]*_cou[11][8][2)  & 
               + faka3 * (_cou[8][8][1) - fakaac * _cou[8][8][2) ) + fakac2 * _cou[11][2][2)
          _cou[27][8][1) = _wmp[2]*_cou[11][8][2) 
          _cou[24][9][1) = _wmp[0]*_cou[11][9][2) 
          _cou[21][9][1) = _wmp[1]*_cou[11][9][2)  & 
               +  faka3 * (_cou[8][9][1) - fakaac * _cou[8][9][2) ) 
          _cou[27][9][1) = _wmp[2]*_cou[11][9][2) +  fakac2 * _cou[11][3][2)
          _cou[26][4][1) = _wmp[0]*_cou[12][4][2) +  fakac * _cou[12][2][2)
          _cou[28][4][1) = _wmp[1]*_cou[12][4][2) +  fakac * _cou[12][1][2)
          _cou[22][4][1) = _wmp[2]*_cou[12][4][2)  & 
               +  faka3 * (_cou[9][4][1) - fakaac * _cou[9][4][2) ) 
          _cou[26][5][1) = _wmp[0]*_cou[12][5][2) +  fakac * _cou[12][3][2)
          _cou[28][5][1) = _wmp[1]*_cou[12][5][2) 
          _cou[22][5][1) = _wmp[2]*_cou[12][5][2)  & 
               +  faka3 * (_cou[9][5][1) - fakaac * _cou[9][5][2) ) + fakac * _cou[12][1][2)
          _cou[26][6][1) = _wmp[0]*_cou[12][6][2) 
          _cou[28][6][1) = _wmp[1]*_cou[12][6][2) +  fakac * _cou[12][3][2)
          _cou[22][6][1) = _wmp[2]*_cou[12][6][2)  & 
               +  faka3 * (_cou[9][6][1) - fakaac * _cou[9][6][2) ) + fakac * _cou[12][2][2)
          _cou[26][7][1) = _wmp[0]*_cou[12][7][2) +  fakac2 * _cou[12][1][2)
          _cou[28][7][1) = _wmp[1]*_cou[12][7][2) 
          _cou[22][7][1) = _wmp[2]*_cou[12][7][2)  & 
               +  faka3 * (_cou[9][7][1) - fakaac * _cou[9][7][2) ) 
          _cou[26][8][1) = _wmp[0]*_cou[12][8][2) 
          _cou[28][8][1) = _wmp[1]*_cou[12][8][2) +  fakac2 * _cou[12][2][2)
          _cou[22][8][1) = _wmp[2]*_cou[12][8][2)  & 
               +  faka3 * (_cou[9][8][1) - fakaac * _cou[9][8][2) ) 
          _cou[26][9][1) = _wmp[0]*_cou[12][9][2) 
          _cou[28][9][1) = _wmp[1]*_cou[12][9][2) 
          _cou[22][9][1) = _wmp[2]*_cou[12][9][2)  & 
               + faka3 * (_cou[9][9][1) - fakaac * _cou[9][9][2) ) + fakac2 * _cou[12][3][2)
          _cou[31][4][1) = _wmp[1]*_cou[13][4][2)  & 
               +  faka * (_cou[7][4][1) - fakaac * _cou[7][4][2) ) + fakac * _cou[13][1][2)
          _cou[32][4][1) = _wmp[2]*_cou[13][4][2) 
          _cou[31][5][1) = _wmp[1]*_cou[13][5][2)  & 
               +  faka * (_cou[7][5][1) - fakaac * _cou[7][5][2) ) 
          _cou[32][5][1) = _wmp[2]*_cou[13][5][2) +  fakac * _cou[13][1][2)
          _cou[31][6][1) = _wmp[1]*_cou[13][6][2)  & 
               +  faka * (_cou[7][6][1) - fakaac * _cou[7][6][2) ) + fakac * _cou[13][3][2)
          _cou[32][6][1) = _wmp[2]*_cou[13][6][2) +  fakac * _cou[13][2][2)
          _cou[31][7][1) = _wmp[1]*_cou[13][7][2)  & 
               +  faka * (_cou[7][7][1) - fakaac * _cou[7][7][2) ) 
          _cou[32][7][1) = _wmp[2]*_cou[13][7][2) 
          _cou[31][8][1) = _wmp[1]*_cou[13][8][2)  & 
               +  faka * (_cou[7][8][1) - fakaac * _cou[7][8][2) ) + fakac2 * _cou[13][2][2)
          _cou[32][8][1) = _wmp[2]*_cou[13][8][2) 
          _cou[31][9][1) = _wmp[1]*_cou[13][9][2)  & 
               +  faka * (_cou[7][9][1) - fakaac * _cou[7][9][2) ) 
          _cou[32][9][1) = _wmp[2]*_cou[13][9][2) +  fakac2 * _cou[13][3][2)
          _cou[33][4][1) = _wmp[2]*_cou[14][4][2) 
          _cou[33][5][1) = _wmp[2]*_cou[14][5][2) +  fakac * _cou[14][1][2)
          _cou[33][6][1) = _wmp[2]*_cou[14][6][2) +  fakac * _cou[14][2][2)
          _cou[33][7][1) = _wmp[2]*_cou[14][7][2) 
          _cou[33][8][1) = _wmp[2]*_cou[14][8][2) 
          _cou[33][9][1) = _wmp[2]*_cou[14][9][2) +  fakac2 * _cou[14][3][2)
          _cou[29][4][1) = _wmp[2]*_cou[15][4][2)  & 
               +  faka * (_cou[7][4][1) - fakaac * _cou[7][4][2) ) 
          _cou[29][5][1) = _wmp[2]*_cou[15][5][2)  & 
               +  faka * (_cou[7][5][1) - fakaac * _cou[7][5][2) ) + fakac * _cou[15][1][2)
          _cou[29][6][1) = _wmp[2]*_cou[15][6][2)  & 
               +  faka * (_cou[7][6][1) - fakaac * _cou[7][6][2) ) + fakac * _cou[15][2][2)
          _cou[29][7][1) = _wmp[2]*_cou[15][7][2)  & 
               +  faka * (_cou[7][7][1) - fakaac * _cou[7][7][2) ) 
          _cou[29][8][1) = _wmp[2]*_cou[15][8][2)  & 
               +  faka * (_cou[7][8][1) - fakaac * _cou[7][8][2) ) 
          _cou[29][9][1) = _wmp[2]*_cou[15][9][2)  & 
               +  faka * (_cou[7][9][1) - fakaac * _cou[7][9][2) ) + fakac2 * _cou[15][3][2)
          _cou[34][4][1) = _wmp[1]*_cou[16][4][2) +  fakac * _cou[16][1][2)
          _cou[34][5][1) = _wmp[1]*_cou[16][5][2) 
          _cou[34][6][1) = _wmp[1]*_cou[16][6][2) +  fakac * _cou[16][3][2)
          _cou[34][7][1) = _wmp[1]*_cou[16][7][2) 
          _cou[34][8][1) = _wmp[1]*_cou[16][8][2) +  fakac2 * _cou[16][2][2)
          _cou[34][9][1) = _wmp[1]*_cou[16][9][2) 
          _cou[30][4][1) = _wmp[2]*_cou[17][4][2)  & 
               +  faka * (_cou[8][4][1) - fakaac * _cou[8][4][2) ) 
          _cou[30][5][1) = _wmp[2]*_cou[17][5][2)  & 
               +  faka * (_cou[8][5][1) - fakaac * _cou[8][5][2) ) + fakac * _cou[17][1][2)
          _cou[30][6][1) = _wmp[2]*_cou[17][6][2)  & 
               +  faka * (_cou[8][6][1) - fakaac * _cou[8][6][2) ) + fakac * _cou[17][2][2)
          _cou[30][7][1) = _wmp[2]*_cou[17][7][2)  & 
               +  faka * (_cou[8][7][1) - fakaac * _cou[8][7][2) ) 
          _cou[30][8][1) = _wmp[2]*_cou[17][8][2)  & 
               +  faka * (_cou[8][8][1) - fakaac * _cou[8][8][2) ) 
          _cou[30][9][1) = _wmp[2]*_cou[17][9][2)  & 
               +  faka * (_cou[8][9][1) - fakaac * _cou[8][9][2) ) + fakac2 * _cou[17][3][2)




          _cou[13][20][0) = _wmp[0]*_cou[4][20][1)  & 
               +  faka * (_cou[2][20][0) - fakaac * _cou[2][20][1) ) + fakac4 * _cou[4][10][1)
          _cou[14][20][0) = _wmp[1]*_cou[4][20][1)  & 
               +  faka * (_cou[1][20][0) - fakaac * _cou[1][20][1) ) 
          _cou[19][20][0) = _wmp[2]*_cou[4][20][1) 
          _cou[13][21][0) = _wmp[0]*_cou[4][21][1)  & 
               +  faka * (_cou[2][21][0) - fakaac * _cou[2][21][1) ) 
          _cou[14][21][0) = _wmp[1]*_cou[4][21][1)  & 
               +  faka * (_cou[1][21][0) - fakaac * _cou[1][21][1) ) + fakac4 * _cou[4][11][1)
          _cou[19][21][0) = _wmp[2]*_cou[4][21][1) 
          _cou[13][22][0) = _wmp[0]*_cou[4][22][1)  & 
               +  faka * (_cou[2][22][0) - fakaac * _cou[2][22][1) ) 
          _cou[14][22][0) = _wmp[1]*_cou[4][22][1)  & 
               +  faka * (_cou[1][22][0) - fakaac * _cou[1][22][1) ) 
          _cou[19][22][0) = _wmp[2]*_cou[4][22][1) +  fakac4 * _cou[4][12][1)
          _cou[13][23][0) = _wmp[0]*_cou[4][23][1)  & 
               +  faka * (_cou[2][23][0) - fakaac * _cou[2][23][1) ) + fakac3 * _cou[4][13][1)
          _cou[14][23][0) = _wmp[1]*_cou[4][23][1)  & 
               +  faka * (_cou[1][23][0) - fakaac * _cou[1][23][1) ) + fakac * _cou[4][10][1)
          _cou[19][23][0) = _wmp[2]*_cou[4][23][1) 
          _cou[13][24][0) = _wmp[0]*_cou[4][24][1)  & 
               +  faka * (_cou[2][24][0) - fakaac * _cou[2][24][1) ) + fakac * _cou[4][11][1)
          _cou[14][24][0) = _wmp[1]*_cou[4][24][1)  & 
               +  faka * (_cou[1][24][0) - fakaac * _cou[1][24][1) ) + fakac3 * _cou[4][14][1)
          _cou[19][24][0) = _wmp[2]*_cou[4][24][1) 
          _cou[13][25][0) = _wmp[0]*_cou[4][25][1)  & 
               +  faka * (_cou[2][25][0) - fakaac * _cou[2][25][1) ) + fakac3 * _cou[4][15][1)
          _cou[14][25][0) = _wmp[1]*_cou[4][25][1)  & 
               +  faka * (_cou[1][25][0) - fakaac * _cou[1][25][1) ) 
          _cou[19][25][0) = _wmp[2]*_cou[4][25][1) +  fakac * _cou[4][10][1)
          _cou[13][26][0) = _wmp[0]*_cou[4][26][1)  & 
               +  faka * (_cou[2][26][0) - fakaac * _cou[2][26][1) ) + fakac * _cou[4][12][1)
          _cou[14][26][0) = _wmp[1]*_cou[4][26][1)  & 
               +  faka * (_cou[1][26][0) - fakaac * _cou[1][26][1) ) 
          _cou[19][26][0) = _wmp[2]*_cou[4][26][1) +  fakac3 * _cou[4][16][1)
          _cou[13][27][0) = _wmp[0]*_cou[4][27][1)  & 
               +  faka * (_cou[2][27][0) - fakaac * _cou[2][27][1) ) 
          _cou[14][27][0) = _wmp[1]*_cou[4][27][1)  & 
               +  faka * (_cou[1][27][0) - fakaac * _cou[1][27][1) ) + fakac3 * _cou[4][17][1)
          _cou[19][27][0) = _wmp[2]*_cou[4][27][1) +  fakac * _cou[4][11][1)
          _cou[13][28][0) = _wmp[0]*_cou[4][28][1)  & 
               +  faka * (_cou[2][28][0) - fakaac * _cou[2][28][1) ) 
          _cou[14][28][0) = _wmp[1]*_cou[4][28][1)  & 
               +  faka * (_cou[1][28][0) - fakaac * _cou[1][28][1) ) + fakac * _cou[4][12][1)
          _cou[19][28][0) = _wmp[2]*_cou[4][28][1) +  fakac3 * _cou[4][18][1)
          _cou[13][29][0) = _wmp[0]*_cou[4][29][1)  & 
               +  faka * (_cou[2][29][0) - fakaac * _cou[2][29][1) ) + fakac2 * _cou[4][16][1)
          _cou[14][29][0) = _wmp[1]*_cou[4][29][1)  & 
               +  faka * (_cou[1][29][0) - fakaac * _cou[1][29][1) ) 
          _cou[19][29][0) = _wmp[2]*_cou[4][29][1) +  fakac2 * _cou[4][15][1)
          _cou[13][30][0) = _wmp[0]*_cou[4][30][1)  & 
               +  faka * (_cou[2][30][0) - fakaac * _cou[2][30][1) ) 
          _cou[14][30][0) = _wmp[1]*_cou[4][30][1)  & 
               +  faka * (_cou[1][30][0) - fakaac * _cou[1][30][1) ) + fakac2 * _cou[4][18][1)
          _cou[19][30][0) = _wmp[2]*_cou[4][30][1) +  fakac2 * _cou[4][17][1)
          _cou[13][31][0) = _wmp[0]*_cou[4][31][1)  & 
               +  faka * (_cou[2][31][0) - fakaac * _cou[2][31][1) ) + fakac2 * _cou[4][14][1)
          _cou[14][31][0) = _wmp[1]*_cou[4][31][1)  & 
               +  faka * (_cou[1][31][0) - fakaac * _cou[1][31][1) ) + fakac2 * _cou[4][13][1)
          _cou[19][31][0) = _wmp[2]*_cou[4][31][1) 
          _cou[13][32][0) = _wmp[0]*_cou[4][32][1)  & 
               +  faka * (_cou[2][32][0) - fakaac * _cou[2][32][1) ) + fakac2 * _cou[4][19][1)
          _cou[14][32][0) = _wmp[1]*_cou[4][32][1)  & 
               +  faka * (_cou[1][32][0) - fakaac * _cou[1][32][1) ) + fakac * _cou[4][15][1)
          _cou[19][32][0) = _wmp[2]*_cou[4][32][1) +  fakac * _cou[4][13][1)
          _cou[13][33][0) = _wmp[0]*_cou[4][33][1)  & 
               +  faka * (_cou[2][33][0) - fakaac * _cou[2][33][1) ) + fakac * _cou[4][17][1)
          _cou[14][33][0) = _wmp[1]*_cou[4][33][1)  & 
               +  faka * (_cou[1][33][0) - fakaac * _cou[1][33][1) ) + fakac2 * _cou[4][19][1)
          _cou[19][33][0) = _wmp[2]*_cou[4][33][1) +  fakac * _cou[4][14][1)
          _cou[13][34][0) = _wmp[0]*_cou[4][34][1)  & 
               +  faka * (_cou[2][34][0) - fakaac * _cou[2][34][1) ) + fakac * _cou[4][18][1)
          _cou[14][34][0) = _wmp[1]*_cou[4][34][1)  & 
               +  faka * (_cou[1][34][0) - fakaac * _cou[1][34][1) ) + fakac * _cou[4][16][1)
          _cou[19][34][0) = _wmp[2]*_cou[4][34][1) +  fakac2 * _cou[4][19][1)
          _cou[15][20][0) = _wmp[0]*_cou[5][20][1)  & 
               +  faka * (_cou[3][20][0) - fakaac * _cou[3][20][1) ) + fakac4 * _cou[5][10][1)
          _cou[16][20][0) = _wmp[2]*_cou[5][20][1)  & 
               +  faka * (_cou[1][20][0) - fakaac * _cou[1][20][1) ) 
          _cou[15][21][0) = _wmp[0]*_cou[5][21][1)  & 
               +  faka * (_cou[3][21][0) - fakaac * _cou[3][21][1) ) 
          _cou[16][21][0) = _wmp[2]*_cou[5][21][1)  & 
               +  faka * (_cou[1][21][0) - fakaac * _cou[1][21][1) ) 
          _cou[15][22][0) = _wmp[0]*_cou[5][22][1)  & 
               +  faka * (_cou[3][22][0) - fakaac * _cou[3][22][1) ) 
          _cou[16][22][0) = _wmp[2]*_cou[5][22][1)  & 
               +  faka * (_cou[1][22][0) - fakaac * _cou[1][22][1) ) + fakac4 * _cou[5][12][1)
          _cou[15][23][0) = _wmp[0]*_cou[5][23][1)  & 
               +  faka * (_cou[3][23][0) - fakaac * _cou[3][23][1) ) + fakac3 * _cou[5][13][1)
          _cou[16][23][0) = _wmp[2]*_cou[5][23][1)  & 
               +  faka * (_cou[1][23][0) - fakaac * _cou[1][23][1) ) 
          _cou[15][24][0) = _wmp[0]*_cou[5][24][1)  & 
               +  faka * (_cou[3][24][0) - fakaac * _cou[3][24][1) ) + fakac * _cou[5][11][1)
          _cou[16][24][0) = _wmp[2]*_cou[5][24][1)  & 
               +  faka * (_cou[1][24][0) - fakaac * _cou[1][24][1) ) 
          _cou[15][25][0) = _wmp[0]*_cou[5][25][1)  & 
               +  faka * (_cou[3][25][0) - fakaac * _cou[3][25][1) ) + fakac3 * _cou[5][15][1)
          _cou[16][25][0) = _wmp[2]*_cou[5][25][1)  & 
               +  faka * (_cou[1][25][0) - fakaac * _cou[1][25][1) ) + fakac * _cou[5][10][1)
          _cou[15][26][0) = _wmp[0]*_cou[5][26][1)  & 
               +  faka * (_cou[3][26][0) - fakaac * _cou[3][26][1) ) + fakac * _cou[5][12][1)
          _cou[16][26][0) = _wmp[2]*_cou[5][26][1)  & 
               +  faka * (_cou[1][26][0) - fakaac * _cou[1][26][1) ) + fakac3 * _cou[5][16][1)
          _cou[15][27][0) = _wmp[0]*_cou[5][27][1)  & 
               +  faka * (_cou[3][27][0) - fakaac * _cou[3][27][1) ) 
          _cou[16][27][0) = _wmp[2]*_cou[5][27][1)  & 
               +  faka * (_cou[1][27][0) - fakaac * _cou[1][27][1) ) + fakac * _cou[5][11][1)
          _cou[15][28][0) = _wmp[0]*_cou[5][28][1)  & 
               +  faka * (_cou[3][28][0) - fakaac * _cou[3][28][1) ) 
          _cou[16][28][0) = _wmp[2]*_cou[5][28][1)  & 
               +  faka * (_cou[1][28][0) - fakaac * _cou[1][28][1) ) + fakac3 * _cou[5][18][1)
          _cou[15][29][0) = _wmp[0]*_cou[5][29][1)  & 
               +  faka * (_cou[3][29][0) - fakaac * _cou[3][29][1) ) + fakac2 * _cou[5][16][1)
          _cou[16][29][0) = _wmp[2]*_cou[5][29][1)  & 
               +  faka * (_cou[1][29][0) - fakaac * _cou[1][29][1) ) + fakac2 * _cou[5][15][1)
          _cou[15][30][0) = _wmp[0]*_cou[5][30][1)  & 
               +  faka * (_cou[3][30][0) - fakaac * _cou[3][30][1) ) 
          _cou[16][30][0) = _wmp[2]*_cou[5][30][1)  & 
               +  faka * (_cou[1][30][0) - fakaac * _cou[1][30][1) ) + fakac2 * _cou[5][17][1)
          _cou[15][31][0) = _wmp[0]*_cou[5][31][1)  & 
               +  faka * (_cou[3][31][0) - fakaac * _cou[3][31][1) ) + fakac2 * _cou[5][14][1)
          _cou[16][31][0) = _wmp[2]*_cou[5][31][1)  & 
               +  faka * (_cou[1][31][0) - fakaac * _cou[1][31][1) ) 
          _cou[15][32][0) = _wmp[0]*_cou[5][32][1)  & 
               +  faka * (_cou[3][32][0) - fakaac * _cou[3][32][1) ) + fakac2 * _cou[5][19][1)
          _cou[16][32][0) = _wmp[2]*_cou[5][32][1)  & 
               +  faka * (_cou[1][32][0) - fakaac * _cou[1][32][1) ) + fakac * _cou[5][13][1)
          _cou[15][33][0) = _wmp[0]*_cou[5][33][1)  & 
               +  faka * (_cou[3][33][0) - fakaac * _cou[3][33][1) ) + fakac * _cou[5][17][1)
          _cou[16][33][0) = _wmp[2]*_cou[5][33][1)  & 
               +  faka * (_cou[1][33][0) - fakaac * _cou[1][33][1) ) + fakac * _cou[5][14][1)
          _cou[15][34][0) = _wmp[0]*_cou[5][34][1)  & 
               +  faka * (_cou[3][34][0) - fakaac * _cou[3][34][1) ) + fakac * _cou[5][18][1)
          _cou[16][34][0) = _wmp[2]*_cou[5][34][1)  & 
               +  faka * (_cou[1][34][0) - fakaac * _cou[1][34][1) ) + fakac2 * _cou[5][19][1)
          _cou[17][20][0) = _wmp[1]*_cou[6][20][1)  & 
               +  faka * (_cou[3][20][0) - fakaac * _cou[3][20][1) ) 
          _cou[18][20][0) = _wmp[2]*_cou[6][20][1)  & 
               +  faka * (_cou[2][20][0) - fakaac * _cou[2][20][1) ) 
          _cou[17][21][0) = _wmp[1]*_cou[6][21][1)  & 
               +  faka * (_cou[3][21][0) - fakaac * _cou[3][21][1) ) + fakac4 * _cou[6][11][1)
          _cou[18][21][0) = _wmp[2]*_cou[6][21][1)  & 
               +  faka * (_cou[2][21][0) - fakaac * _cou[2][21][1) ) 
          _cou[17][22][0) = _wmp[1]*_cou[6][22][1)  & 
               +  faka * (_cou[3][22][0) - fakaac * _cou[3][22][1) ) 
          _cou[18][22][0) = _wmp[2]*_cou[6][22][1)  & 
               +  faka * (_cou[2][22][0) - fakaac * _cou[2][22][1) ) + fakac4 * _cou[6][12][1)
          _cou[17][23][0) = _wmp[1]*_cou[6][23][1)  & 
               +  faka * (_cou[3][23][0) - fakaac * _cou[3][23][1) ) + fakac * _cou[6][10][1)
          _cou[18][23][0) = _wmp[2]*_cou[6][23][1)  & 
               +  faka * (_cou[2][23][0) - fakaac * _cou[2][23][1) ) 
          _cou[17][24][0) = _wmp[1]*_cou[6][24][1)  & 
               +  faka * (_cou[3][24][0) - fakaac * _cou[3][24][1) ) + fakac3 * _cou[6][14][1)
          _cou[18][24][0) = _wmp[2]*_cou[6][24][1)  & 
               +  faka * (_cou[2][24][0) - fakaac * _cou[2][24][1) ) 
          _cou[17][25][0) = _wmp[1]*_cou[6][25][1)  & 
               +  faka * (_cou[3][25][0) - fakaac * _cou[3][25][1) ) 
          _cou[18][25][0) = _wmp[2]*_cou[6][25][1)  & 
               +  faka * (_cou[2][25][0) - fakaac * _cou[2][25][1) ) + fakac * _cou[6][10][1)
          _cou[17][26][0) = _wmp[1]*_cou[6][26][1)  & 
               +  faka * (_cou[3][26][0) - fakaac * _cou[3][26][1) ) 
          _cou[18][26][0) = _wmp[2]*_cou[6][26][1)  & 
               +  faka * (_cou[2][26][0) - fakaac * _cou[2][26][1) ) + fakac3 * _cou[6][16][1)
          _cou[17][27][0) = _wmp[1]*_cou[6][27][1)  & 
               +  faka * (_cou[3][27][0) - fakaac * _cou[3][27][1) ) + fakac3 * _cou[6][17][1)
          _cou[18][27][0) = _wmp[2]*_cou[6][27][1)  & 
               +  faka * (_cou[2][27][0) - fakaac * _cou[2][27][1) ) + fakac * _cou[6][11][1)
          _cou[17][28][0) = _wmp[1]*_cou[6][28][1)  & 
               +  faka * (_cou[3][28][0) - fakaac * _cou[3][28][1) ) + fakac * _cou[6][12][1)
          _cou[18][28][0) = _wmp[2]*_cou[6][28][1)  & 
               +  faka * (_cou[2][28][0) - fakaac * _cou[2][28][1) ) + fakac3 * _cou[6][18][1)
          _cou[17][29][0) = _wmp[1]*_cou[6][29][1)  & 
               +  faka * (_cou[3][29][0) - fakaac * _cou[3][29][1) ) 
          _cou[18][29][0) = _wmp[2]*_cou[6][29][1)  & 
               +  faka * (_cou[2][29][0) - fakaac * _cou[2][29][1) ) + fakac2 * _cou[6][15][1)
          _cou[17][30][0) = _wmp[1]*_cou[6][30][1)  & 
               +  faka * (_cou[3][30][0) - fakaac * _cou[3][30][1) ) + fakac2 * _cou[6][18][1)
          _cou[18][30][0) = _wmp[2]*_cou[6][30][1)  & 
               +  faka * (_cou[2][30][0) - fakaac * _cou[2][30][1) ) + fakac2 * _cou[6][17][1)
          _cou[17][31][0) = _wmp[1]*_cou[6][31][1)  & 
               +  faka * (_cou[3][31][0) - fakaac * _cou[3][31][1) ) + fakac2 * _cou[6][13][1)
          _cou[18][31][0) = _wmp[2]*_cou[6][31][1)  & 
               +  faka * (_cou[2][31][0) - fakaac * _cou[2][31][1) ) 
          _cou[17][32][0) = _wmp[1]*_cou[6][32][1)  & 
               +  faka * (_cou[3][32][0) - fakaac * _cou[3][32][1) ) + fakac * _cou[6][15][1)
          _cou[18][32][0) = _wmp[2]*_cou[6][32][1)  & 
               +  faka * (_cou[2][32][0) - fakaac * _cou[2][32][1) ) + fakac * _cou[6][13][1)
          _cou[17][33][0) = _wmp[1]*_cou[6][33][1)  & 
               +  faka * (_cou[3][33][0) - fakaac * _cou[3][33][1) ) + fakac2 * _cou[6][19][1)
          _cou[18][33][0) = _wmp[2]*_cou[6][33][1)  & 
               +  faka * (_cou[2][33][0) - fakaac * _cou[2][33][1) ) + fakac * _cou[6][14][1)
          _cou[17][34][0) = _wmp[1]*_cou[6][34][1)  & 
               +  faka * (_cou[3][34][0) - fakaac * _cou[3][34][1) ) + fakac * _cou[6][16][1)
          _cou[18][34][0) = _wmp[2]*_cou[6][34][1)  & 
               +  faka * (_cou[2][34][0) - fakaac * _cou[2][34][1) ) + fakac2 * _cou[6][19][1)
          _cou[10][20][0) = _wmp[0]*_cou[7][20][1)  & 
               + faka2 * (_cou[1][20][0) - fakaac * _cou[1][20][1) ) + fakac4 * _cou[7][10][1)
          _cou[10][21][0) = _wmp[0]*_cou[7][21][1)  & 
               +  faka2 * (_cou[1][21][0) - fakaac * _cou[1][21][1) ) 
          _cou[10][22][0) = _wmp[0]*_cou[7][22][1)  & 
               +  faka2 * (_cou[1][22][0) - fakaac * _cou[1][22][1) ) 
          _cou[10][23][0) = _wmp[0]*_cou[7][23][1)  & 
               + faka2 * (_cou[1][23][0) - fakaac * _cou[1][23][1) ) + fakac3 * _cou[7][13][1)
          _cou[10][24][0) = _wmp[0]*_cou[7][24][1)  & 
               +  faka2 * (_cou[1][24][0) - fakaac * _cou[1][24][1) ) + fakac * _cou[7][11][1)
          _cou[10][25][0) = _wmp[0]*_cou[7][25][1)  & 
               + faka2 * (_cou[1][25][0) - fakaac * _cou[1][25][1) ) + fakac3 * _cou[7][15][1)
          _cou[10][26][0) = _wmp[0]*_cou[7][26][1)  & 
               +  faka2 * (_cou[1][26][0) - fakaac * _cou[1][26][1) ) + fakac * _cou[7][12][1)
          _cou[10][27][0) = _wmp[0]*_cou[7][27][1)  & 
               +  faka2 * (_cou[1][27][0) - fakaac * _cou[1][27][1) ) 
          _cou[10][28][0) = _wmp[0]*_cou[7][28][1)  & 
               +  faka2 * (_cou[1][28][0) - fakaac * _cou[1][28][1) ) 
          _cou[10][29][0) = _wmp[0]*_cou[7][29][1)  & 
               + faka2 * (_cou[1][29][0) - fakaac * _cou[1][29][1) ) + fakac2 * _cou[7][16][1)
          _cou[10][30][0) = _wmp[0]*_cou[7][30][1)  & 
               +  faka2 * (_cou[1][30][0) - fakaac * _cou[1][30][1) ) 
          _cou[10][31][0) = _wmp[0]*_cou[7][31][1)  & 
               + faka2 * (_cou[1][31][0) - fakaac * _cou[1][31][1) ) + fakac2 * _cou[7][14][1)
          _cou[10][32][0) = _wmp[0]*_cou[7][32][1)  & 
               + faka2 * (_cou[1][32][0) - fakaac * _cou[1][32][1) ) + fakac2 * _cou[7][19][1)
          _cou[10][33][0) = _wmp[0]*_cou[7][33][1)  & 
               +  faka2 * (_cou[1][33][0) - fakaac * _cou[1][33][1) ) + fakac * _cou[7][17][1)
          _cou[10][34][0) = _wmp[0]*_cou[7][34][1)  & 
               +  faka2 * (_cou[1][34][0) - fakaac * _cou[1][34][1) ) + fakac * _cou[7][18][1)
          _cou[11][20][0) = _wmp[1]*_cou[8][20][1)  & 
               +  faka2 * (_cou[2][20][0) - fakaac * _cou[2][20][1) ) 
          _cou[11][21][0) = _wmp[1]*_cou[8][21][1)  & 
               + faka2 * (_cou[2][21][0) - fakaac * _cou[2][21][1) ) + fakac4 * _cou[8][11][1)
          _cou[11][22][0) = _wmp[1]*_cou[8][22][1)  & 
               +  faka2 * (_cou[2][22][0) - fakaac * _cou[2][22][1) ) 
          _cou[11][23][0) = _wmp[1]*_cou[8][23][1)  & 
               +  faka2 * (_cou[2][23][0) - fakaac * _cou[2][23][1) ) + fakac * _cou[8][10][1)
          _cou[11][24][0) = _wmp[1]*_cou[8][24][1)  & 
               + faka2 * (_cou[2][24][0) - fakaac * _cou[2][24][1) ) + fakac3 * _cou[8][14][1)
          _cou[11][25][0) = _wmp[1]*_cou[8][25][1)  & 
               +  faka2 * (_cou[2][25][0) - fakaac * _cou[2][25][1) ) 
          _cou[11][26][0) = _wmp[1]*_cou[8][26][1)  & 
               +  faka2 * (_cou[2][26][0) - fakaac * _cou[2][26][1) ) 
          _cou[11][27][0) = _wmp[1]*_cou[8][27][1)  & 
               + faka2 * (_cou[2][27][0) - fakaac * _cou[2][27][1) ) + fakac3 * _cou[8][17][1)
          _cou[11][28][0) = _wmp[1]*_cou[8][28][1)  & 
               +  faka2 * (_cou[2][28][0) - fakaac * _cou[2][28][1) ) + fakac * _cou[8][12][1)
          _cou[11][29][0) = _wmp[1]*_cou[8][29][1)  & 
               +  faka2 * (_cou[2][29][0) - fakaac * _cou[2][29][1) ) 
          _cou[11][30][0) = _wmp[1]*_cou[8][30][1)  & 
               + faka2 * (_cou[2][30][0) - fakaac * _cou[2][30][1) ) + fakac2 * _cou[8][18][1)
          _cou[11][31][0) = _wmp[1]*_cou[8][31][1)  & 
               + faka2 * (_cou[2][31][0) - fakaac * _cou[2][31][1) ) + fakac2 * _cou[8][13][1)
          _cou[11][32][0) = _wmp[1]*_cou[8][32][1)  & 
               +  faka2 * (_cou[2][32][0) - fakaac * _cou[2][32][1) ) + fakac * _cou[8][15][1)
          _cou[11][33][0) = _wmp[1]*_cou[8][33][1)  & 
               + faka2 * (_cou[2][33][0) - fakaac * _cou[2][33][1) ) + fakac2 * _cou[8][19][1)
          _cou[11][34][0) = _wmp[1]*_cou[8][34][1)  & 
               +  faka2 * (_cou[2][34][0) - fakaac * _cou[2][34][1) ) + fakac * _cou[8][16][1)
          _cou[12][20][0) = _wmp[2]*_cou[9][20][1)  & 
               +  faka2 * (_cou[3][20][0) - fakaac * _cou[3][20][1) ) 
          _cou[12][21][0) = _wmp[2]*_cou[9][21][1)  & 
               +  faka2 * (_cou[3][21][0) - fakaac * _cou[3][21][1) ) 
          _cou[12][22][0) = _wmp[2]*_cou[9][22][1)  & 
               + faka2 * (_cou[3][22][0) - fakaac * _cou[3][22][1) ) + fakac4 * _cou[9][12][1)
          _cou[12][23][0) = _wmp[2]*_cou[9][23][1)  & 
               +  faka2 * (_cou[3][23][0) - fakaac * _cou[3][23][1) ) 
          _cou[12][24][0) = _wmp[2]*_cou[9][24][1)  & 
               +  faka2 * (_cou[3][24][0) - fakaac * _cou[3][24][1) ) 
          _cou[12][25][0) = _wmp[2]*_cou[9][25][1)  & 
               +  faka2 * (_cou[3][25][0) - fakaac * _cou[3][25][1) ) + fakac * _cou[9][10][1)
          _cou[12][26][0) = _wmp[2]*_cou[9][26][1)  & 
               + faka2 * (_cou[3][26][0) - fakaac * _cou[3][26][1) ) + fakac3 * _cou[9][16][1)
          _cou[12][27][0) = _wmp[2]*_cou[9][27][1)  & 
               +  faka2 * (_cou[3][27][0) - fakaac * _cou[3][27][1) ) + fakac * _cou[9][11][1)
          _cou[12][28][0) = _wmp[2]*_cou[9][28][1)  & 
               + faka2 * (_cou[3][28][0) - fakaac * _cou[3][28][1) ) + fakac3 * _cou[9][18][1)
          _cou[12][29][0) = _wmp[2]*_cou[9][29][1)  & 
               + faka2 * (_cou[3][29][0) - fakaac * _cou[3][29][1) ) + fakac2 * _cou[9][15][1)
          _cou[12][30][0) = _wmp[2]*_cou[9][30][1)  & 
               + faka2 * (_cou[3][30][0) - fakaac * _cou[3][30][1) ) + fakac2 * _cou[9][17][1)
          _cou[12][31][0) = _wmp[2]*_cou[9][31][1)  & 
               +  faka2 * (_cou[3][31][0) - fakaac * _cou[3][31][1) ) 
          _cou[12][32][0) = _wmp[2]*_cou[9][32][1)  & 
               +  faka2 * (_cou[3][32][0) - fakaac * _cou[3][32][1) ) + fakac * _cou[9][13][1)
          _cou[12][33][0) = _wmp[2]*_cou[9][33][1)  & 
               +  faka2 * (_cou[3][33][0) - fakaac * _cou[3][33][1) ) + fakac * _cou[9][14][1)
          _cou[12][34][0) = _wmp[2]*_cou[9][34][1)  & 
               + faka2 * (_cou[3][34][0) - fakaac * _cou[3][34][1) ) + fakac2 * _cou[9][19][1)




          _cou[20][10][0) = _wmp[0]*_cou[10][10][1)  & 
               + faka3 * (_cou[7][10][0) - fakaac * _cou[7][10][1) ) + fakac3 * _cou[10][7][1)
          _cou[23][10][0) = _wmp[1]*_cou[10][10][1) 
          _cou[25][10][0) = _wmp[2]*_cou[10][10][1) 
          _cou[20][11][0) = _wmp[0]*_cou[10][11][1)  & 
               +  faka3 * (_cou[7][11][0) - fakaac * _cou[7][11][1) ) 
          _cou[23][11][0) = _wmp[1]*_cou[10][11][1) +  fakac3 * _cou[10][8][1)
          _cou[25][11][0) = _wmp[2]*_cou[10][11][1) 
          _cou[20][12][0) = _wmp[0]*_cou[10][12][1)  & 
               +  faka3 * (_cou[7][12][0) - fakaac * _cou[7][12][1) ) 
          _cou[23][12][0) = _wmp[1]*_cou[10][12][1) 
          _cou[25][12][0) = _wmp[2]*_cou[10][12][1) +  fakac3 * _cou[10][9][1)
          _cou[20][13][0) = _wmp[0]*_cou[10][13][1)  & 
               + faka3 * (_cou[7][13][0) - fakaac * _cou[7][13][1) ) + fakac2 * _cou[10][4][1)
          _cou[23][13][0) = _wmp[1]*_cou[10][13][1) +  fakac * _cou[10][7][1)
          _cou[25][13][0) = _wmp[2]*_cou[10][13][1) 
          _cou[20][14][0) = _wmp[0]*_cou[10][14][1)  & 
               +  faka3 * (_cou[7][14][0) - fakaac * _cou[7][14][1) ) + fakac * _cou[10][8][1)
          _cou[23][14][0) = _wmp[1]*_cou[10][14][1) +  fakac2 * _cou[10][4][1)
          _cou[25][14][0) = _wmp[2]*_cou[10][14][1) 
          _cou[20][15][0) = _wmp[0]*_cou[10][15][1)  & 
               + faka3 * (_cou[7][15][0) - fakaac * _cou[7][15][1) ) + fakac2 * _cou[10][5][1)
          _cou[23][15][0) = _wmp[1]*_cou[10][15][1) 
          _cou[25][15][0) = _wmp[2]*_cou[10][15][1) +  fakac * _cou[10][7][1)
          _cou[20][16][0) = _wmp[0]*_cou[10][16][1)  & 
               +  faka3 * (_cou[7][16][0) - fakaac * _cou[7][16][1) ) + fakac * _cou[10][9][1)
          _cou[23][16][0) = _wmp[1]*_cou[10][16][1) 
          _cou[25][16][0) = _wmp[2]*_cou[10][16][1) +  fakac2 * _cou[10][5][1)
          _cou[20][17][0) = _wmp[0]*_cou[10][17][1)  & 
               +  faka3 * (_cou[7][17][0) - fakaac * _cou[7][17][1) ) 
          _cou[23][17][0) = _wmp[1]*_cou[10][17][1) +  fakac2 * _cou[10][6][1)
          _cou[25][17][0) = _wmp[2]*_cou[10][17][1) +  fakac * _cou[10][8][1)
          _cou[20][18][0) = _wmp[0]*_cou[10][18][1)  & 
               +  faka3 * (_cou[7][18][0) - fakaac * _cou[7][18][1) ) 
          _cou[23][18][0) = _wmp[1]*_cou[10][18][1) +  fakac * _cou[10][9][1)
          _cou[25][18][0) = _wmp[2]*_cou[10][18][1) +  fakac2 * _cou[10][6][1)
          _cou[20][19][0) = _wmp[0]*_cou[10][19][1)  & 
               +  faka3 * (_cou[7][19][0) - fakaac * _cou[7][19][1) ) + fakac * _cou[10][6][1)
          _cou[23][19][0) = _wmp[1]*_cou[10][19][1) +  fakac * _cou[10][5][1)
          _cou[25][19][0) = _wmp[2]*_cou[10][19][1) +  fakac * _cou[10][4][1)
          _cou[24][10][0) = _wmp[0]*_cou[11][10][1) +  fakac3 * _cou[11][7][1)
          _cou[21][10][0) = _wmp[1]*_cou[11][10][1)  & 
               +  faka3 * (_cou[8][10][0) - fakaac * _cou[8][10][1) ) 
          _cou[27][10][0) = _wmp[2]*_cou[11][10][1) 
          _cou[24][11][0) = _wmp[0]*_cou[11][11][1) 
          _cou[21][11][0) = _wmp[1]*_cou[11][11][1)  & 
               + faka3 * (_cou[8][11][0) - fakaac * _cou[8][11][1) ) + fakac3 * _cou[11][8][1)
          _cou[27][11][0) = _wmp[2]*_cou[11][11][1) 
          _cou[24][12][0) = _wmp[0]*_cou[11][12][1) 
          _cou[21][12][0) = _wmp[1]*_cou[11][12][1)  & 
               +  faka3 * (_cou[8][12][0) - fakaac * _cou[8][12][1) ) 
          _cou[27][12][0) = _wmp[2]*_cou[11][12][1) +  fakac3 * _cou[11][9][1)
          _cou[24][13][0) = _wmp[0]*_cou[11][13][1) +  fakac2 * _cou[11][4][1)
          _cou[21][13][0) = _wmp[1]*_cou[11][13][1)  & 
               +  faka3 * (_cou[8][13][0) - fakaac * _cou[8][13][1) ) + fakac * _cou[11][7][1)
          _cou[27][13][0) = _wmp[2]*_cou[11][13][1) 
          _cou[24][14][0) = _wmp[0]*_cou[11][14][1) +  fakac * _cou[11][8][1)
          _cou[21][14][0) = _wmp[1]*_cou[11][14][1)  & 
               + faka3 * (_cou[8][14][0) - fakaac * _cou[8][14][1) ) + fakac2 * _cou[11][4][1)
          _cou[27][14][0) = _wmp[2]*_cou[11][14][1) 
          _cou[24][15][0) = _wmp[0]*_cou[11][15][1) +  fakac2 * _cou[11][5][1)
          _cou[21][15][0) = _wmp[1]*_cou[11][15][1)  & 
               +  faka3 * (_cou[8][15][0) - fakaac * _cou[8][15][1) ) 
          _cou[27][15][0) = _wmp[2]*_cou[11][15][1) +  fakac * _cou[11][7][1)
          _cou[24][16][0) = _wmp[0]*_cou[11][16][1) +  fakac * _cou[11][9][1)
          _cou[21][16][0) = _wmp[1]*_cou[11][16][1)  & 
               +  faka3 * (_cou[8][16][0) - fakaac * _cou[8][16][1) ) 
          _cou[27][16][0) = _wmp[2]*_cou[11][16][1) +  fakac2 * _cou[11][5][1)
          _cou[24][17][0) = _wmp[0]*_cou[11][17][1) 
          _cou[21][17][0) = _wmp[1]*_cou[11][17][1)  & 
               + faka3 * (_cou[8][17][0) - fakaac * _cou[8][17][1) ) + fakac2 * _cou[11][6][1)
          _cou[27][17][0) = _wmp[2]*_cou[11][17][1) +  fakac * _cou[11][8][1)
          _cou[24][18][0) = _wmp[0]*_cou[11][18][1) 
          _cou[21][18][0) = _wmp[1]*_cou[11][18][1)  & 
               +  faka3 * (_cou[8][18][0) - fakaac * _cou[8][18][1) ) + fakac * _cou[11][9][1)
          _cou[27][18][0) = _wmp[2]*_cou[11][18][1) +  fakac2 * _cou[11][6][1)
          _cou[24][19][0) = _wmp[0]*_cou[11][19][1) +  fakac * _cou[11][6][1)
          _cou[21][19][0) = _wmp[1]*_cou[11][19][1)  & 
               +  faka3 * (_cou[8][19][0) - fakaac * _cou[8][19][1) ) + fakac * _cou[11][5][1)
          _cou[27][19][0) = _wmp[2]*_cou[11][19][1) +  fakac * _cou[11][4][1)
          _cou[26][10][0) = _wmp[0]*_cou[12][10][1) +  fakac3 * _cou[12][7][1)
          _cou[28][10][0) = _wmp[1]*_cou[12][10][1) 
          _cou[22][10][0) = _wmp[2]*_cou[12][10][1)  & 
               +  faka3 * (_cou[9][10][0) - fakaac * _cou[9][10][1) ) 
          _cou[26][11][0) = _wmp[0]*_cou[12][11][1) 
          _cou[28][11][0) = _wmp[1]*_cou[12][11][1) +  fakac3 * _cou[12][8][1)
          _cou[22][11][0) = _wmp[2]*_cou[12][11][1)  & 
               +  faka3 * (_cou[9][11][0) - fakaac * _cou[9][11][1) ) 
          _cou[26][12][0) = _wmp[0]*_cou[12][12][1) 
          _cou[28][12][0) = _wmp[1]*_cou[12][12][1) 
          _cou[22][12][0) = _wmp[2]*_cou[12][12][1)  & 
               + faka3 * (_cou[9][12][0) - fakaac * _cou[9][12][1) ) + fakac3 * _cou[12][9][1)
          _cou[26][13][0) = _wmp[0]*_cou[12][13][1) +  fakac2 * _cou[12][4][1)
          _cou[28][13][0) = _wmp[1]*_cou[12][13][1) +  fakac * _cou[12][7][1)
          _cou[22][13][0) = _wmp[2]*_cou[12][13][1)  & 
               +  faka3 * (_cou[9][13][0) - fakaac * _cou[9][13][1) ) 
          _cou[26][14][0) = _wmp[0]*_cou[12][14][1) +  fakac * _cou[12][8][1)
          _cou[28][14][0) = _wmp[1]*_cou[12][14][1) +  fakac2 * _cou[12][4][1)
          _cou[22][14][0) = _wmp[2]*_cou[12][14][1)  & 
               +  faka3 * (_cou[9][14][0) - fakaac * _cou[9][14][1) ) 
          _cou[26][15][0) = _wmp[0]*_cou[12][15][1) +  fakac2 * _cou[12][5][1)
          _cou[28][15][0) = _wmp[1]*_cou[12][15][1) 
          _cou[22][15][0) = _wmp[2]*_cou[12][15][1)  & 
               +  faka3 * (_cou[9][15][0) - fakaac * _cou[9][15][1) ) + fakac * _cou[12][7][1)
          _cou[26][16][0) = _wmp[0]*_cou[12][16][1) +  fakac * _cou[12][9][1)
          _cou[28][16][0) = _wmp[1]*_cou[12][16][1) 
          _cou[22][16][0) = _wmp[2]*_cou[12][16][1)  & 
               + faka3 * (_cou[9][16][0) - fakaac * _cou[9][16][1) ) + fakac2 * _cou[12][5][1)
          _cou[26][17][0) = _wmp[0]*_cou[12][17][1) 
          _cou[28][17][0) = _wmp[1]*_cou[12][17][1) +  fakac2 * _cou[12][6][1)
          _cou[22][17][0) = _wmp[2]*_cou[12][17][1)  & 
               +  faka3 * (_cou[9][17][0) - fakaac * _cou[9][17][1) ) + fakac * _cou[12][8][1)
          _cou[26][18][0) = _wmp[0]*_cou[12][18][1) 
          _cou[28][18][0) = _wmp[1]*_cou[12][18][1) +  fakac * _cou[12][9][1)
          _cou[22][18][0) = _wmp[2]*_cou[12][18][1)  & 
               + faka3 * (_cou[9][18][0) - fakaac * _cou[9][18][1) ) + fakac2 * _cou[12][6][1)
          _cou[26][19][0) = _wmp[0]*_cou[12][19][1) +  fakac * _cou[12][6][1)
          _cou[28][19][0) = _wmp[1]*_cou[12][19][1) +  fakac * _cou[12][5][1)
          _cou[22][19][0) = _wmp[2]*_cou[12][19][1)  & 
               +  faka3 * (_cou[9][19][0) - fakaac * _cou[9][19][1) ) + fakac * _cou[12][4][1)
          _cou[31][10][0) = _wmp[1]*_cou[13][10][1)  & 
               +  faka * (_cou[7][10][0) - fakaac * _cou[7][10][1) ) 
          _cou[32][10][0) = _wmp[2]*_cou[13][10][1) 
          _cou[31][11][0) = _wmp[1]*_cou[13][11][1)  & 
               +  faka * (_cou[7][11][0) - fakaac * _cou[7][11][1) ) + fakac3 * _cou[13][8][1)
          _cou[32][11][0) = _wmp[2]*_cou[13][11][1) 
          _cou[31][12][0) = _wmp[1]*_cou[13][12][1)  & 
               +  faka * (_cou[7][12][0) - fakaac * _cou[7][12][1) ) 
          _cou[32][12][0) = _wmp[2]*_cou[13][12][1) +  fakac3 * _cou[13][9][1)
          _cou[31][13][0) = _wmp[1]*_cou[13][13][1)  & 
               +  faka * (_cou[7][13][0) - fakaac * _cou[7][13][1) ) + fakac * _cou[13][7][1)
          _cou[32][13][0) = _wmp[2]*_cou[13][13][1) 
          _cou[31][14][0) = _wmp[1]*_cou[13][14][1)  & 
               +  faka * (_cou[7][14][0) - fakaac * _cou[7][14][1) ) + fakac2 * _cou[13][4][1)
          _cou[32][14][0) = _wmp[2]*_cou[13][14][1) 
          _cou[31][15][0) = _wmp[1]*_cou[13][15][1)  & 
               +  faka * (_cou[7][15][0) - fakaac * _cou[7][15][1) ) 
          _cou[32][15][0) = _wmp[2]*_cou[13][15][1) +  fakac * _cou[13][7][1)
          _cou[31][16][0) = _wmp[1]*_cou[13][16][1)  & 
               +  faka * (_cou[7][16][0) - fakaac * _cou[7][16][1) ) 
          _cou[32][16][0) = _wmp[2]*_cou[13][16][1) +  fakac2 * _cou[13][5][1)
          _cou[31][17][0) = _wmp[1]*_cou[13][17][1)  & 
               +  faka * (_cou[7][17][0) - fakaac * _cou[7][17][1) ) + fakac2 * _cou[13][6][1)
          _cou[32][17][0) = _wmp[2]*_cou[13][17][1) +  fakac * _cou[13][8][1)
          _cou[31][18][0) = _wmp[1]*_cou[13][18][1)  & 
               +  faka * (_cou[7][18][0) - fakaac * _cou[7][18][1) ) + fakac * _cou[13][9][1)
          _cou[32][18][0) = _wmp[2]*_cou[13][18][1) +  fakac2 * _cou[13][6][1)
          _cou[31][19][0) = _wmp[1]*_cou[13][19][1)  & 
               +  faka * (_cou[7][19][0) - fakaac * _cou[7][19][1) ) + fakac * _cou[13][5][1)
          _cou[32][19][0) = _wmp[2]*_cou[13][19][1) +  fakac * _cou[13][4][1)
          _cou[33][10][0) = _wmp[2]*_cou[14][10][1) 
          _cou[33][11][0) = _wmp[2]*_cou[14][11][1) 
          _cou[33][12][0) = _wmp[2]*_cou[14][12][1) +  fakac3 * _cou[14][9][1)
          _cou[33][13][0) = _wmp[2]*_cou[14][13][1) 
          _cou[33][14][0) = _wmp[2]*_cou[14][14][1) 
          _cou[33][15][0) = _wmp[2]*_cou[14][15][1) +  fakac * _cou[14][7][1)
          _cou[33][16][0) = _wmp[2]*_cou[14][16][1) +  fakac2 * _cou[14][5][1)
          _cou[33][17][0) = _wmp[2]*_cou[14][17][1) +  fakac * _cou[14][8][1)
          _cou[33][18][0) = _wmp[2]*_cou[14][18][1) +  fakac2 * _cou[14][6][1)
          _cou[33][19][0) = _wmp[2]*_cou[14][19][1) +  fakac * _cou[14][4][1)
          _cou[29][10][0) = _wmp[2]*_cou[15][10][1)  & 
               +  faka * (_cou[7][10][0) - fakaac * _cou[7][10][1) ) 
          _cou[29][11][0) = _wmp[2]*_cou[15][11][1)  & 
               +  faka * (_cou[7][11][0) - fakaac * _cou[7][11][1) ) 
          _cou[29][12][0) = _wmp[2]*_cou[15][12][1)  & 
               +  faka * (_cou[7][12][0) - fakaac * _cou[7][12][1) ) + fakac3 * _cou[15][9][1)
          _cou[29][13][0) = _wmp[2]*_cou[15][13][1)  & 
               +  faka * (_cou[7][13][0) - fakaac * _cou[7][13][1) ) 
          _cou[29][14][0) = _wmp[2]*_cou[15][14][1)  & 
               +  faka * (_cou[7][14][0) - fakaac * _cou[7][14][1) ) 
          _cou[29][15][0) = _wmp[2]*_cou[15][15][1)  & 
               +  faka * (_cou[7][15][0) - fakaac * _cou[7][15][1) ) + fakac * _cou[15][7][1)
          _cou[29][16][0) = _wmp[2]*_cou[15][16][1)  & 
               +  faka * (_cou[7][16][0) - fakaac * _cou[7][16][1) ) + fakac2 * _cou[15][5][1)
          _cou[29][17][0) = _wmp[2]*_cou[15][17][1)  & 
               +  faka * (_cou[7][17][0) - fakaac * _cou[7][17][1) ) + fakac * _cou[15][8][1)
          _cou[29][18][0) = _wmp[2]*_cou[15][18][1)  & 
               +  faka * (_cou[7][18][0) - fakaac * _cou[7][18][1) ) + fakac2 * _cou[15][6][1)
          _cou[29][19][0) = _wmp[2]*_cou[15][19][1)  & 
               +  faka * (_cou[7][19][0) - fakaac * _cou[7][19][1) ) + fakac * _cou[15][4][1)
          _cou[34][10][0) = _wmp[1]*_cou[16][10][1) 
          _cou[34][11][0) = _wmp[1]*_cou[16][11][1) +  fakac3 * _cou[16][8][1)
          _cou[34][12][0) = _wmp[1]*_cou[16][12][1) 
          _cou[34][13][0) = _wmp[1]*_cou[16][13][1) +  fakac * _cou[16][7][1)
          _cou[34][14][0) = _wmp[1]*_cou[16][14][1) +  fakac2 * _cou[16][4][1)
          _cou[34][15][0) = _wmp[1]*_cou[16][15][1) 
          _cou[34][16][0) = _wmp[1]*_cou[16][16][1) 
          _cou[34][17][0) = _wmp[1]*_cou[16][17][1) +  fakac2 * _cou[16][6][1)
          _cou[34][18][0) = _wmp[1]*_cou[16][18][1) +  fakac * _cou[16][9][1)
          _cou[34][19][0) = _wmp[1]*_cou[16][19][1) +  fakac * _cou[16][5][1)
          _cou[30][10][0) = _wmp[2]*_cou[17][10][1)  & 
               +  faka * (_cou[8][10][0) - fakaac * _cou[8][10][1) ) 
          _cou[30][11][0) = _wmp[2]*_cou[17][11][1)  & 
               +  faka * (_cou[8][11][0) - fakaac * _cou[8][11][1) ) 
          _cou[30][12][0) = _wmp[2]*_cou[17][12][1)  & 
               +  faka * (_cou[8][12][0) - fakaac * _cou[8][12][1) ) + fakac3 * _cou[17][9][1)
          _cou[30][13][0) = _wmp[2]*_cou[17][13][1)  & 
               +  faka * (_cou[8][13][0) - fakaac * _cou[8][13][1) ) 
          _cou[30][14][0) = _wmp[2]*_cou[17][14][1)  & 
               +  faka * (_cou[8][14][0) - fakaac * _cou[8][14][1) ) 
          _cou[30][15][0) = _wmp[2]*_cou[17][15][1)  & 
               +  faka * (_cou[8][15][0) - fakaac * _cou[8][15][1) ) + fakac * _cou[17][7][1)
          _cou[30][16][0) = _wmp[2]*_cou[17][16][1)  & 
               +  faka * (_cou[8][16][0) - fakaac * _cou[8][16][1) ) + fakac2 * _cou[17][5][1)
          _cou[30][17][0) = _wmp[2]*_cou[17][17][1)  & 
               +  faka * (_cou[8][17][0) - fakaac * _cou[8][17][1) ) + fakac * _cou[17][8][1)
          _cou[30][18][0) = _wmp[2]*_cou[17][18][1)  & 
               +  faka * (_cou[8][18][0) - fakaac * _cou[8][18][1) ) + fakac2 * _cou[17][6][1)
          _cou[30][19][0) = _wmp[2]*_cou[17][19][1)  & 
               +  faka * (_cou[8][19][0) - fakaac * _cou[8][19][1) ) + fakac * _cou[17][4][1)

       endif


       if(abs(lmax1)+abs(lmax2).ge.8) then


          _cou[1][0][7) = _wmp[0]*_cou[0][0][8) 
          _cou[0][1][7) = _wmq[0]*_cou[0][0][8) 
          _cou[2][0][7) = _wmp[1]*_cou[0][0][8) 
          _cou[0][2][7) = _wmq[1]*_cou[0][0][8) 
          _cou[3][0][7) = _wmp[2]*_cou[0][0][8) 
          _cou[0][3][7) = _wmq[2]*_cou[0][0][8) 
          _cou[1][1][6) = _wmp[0]*_cou[0][1][7) +  fakac * _cou[0][0][7)
          _cou[0][7][6) = _wmq[0]*_cou[0][1][7)  & 
               +  fakc * (_cou[0][0][6) - fakaca * _cou[0][0][7) ) 
          _cou[2][1][6) = _wmp[1]*_cou[0][1][7) 
          _cou[0][4][6) = _wmq[1]*_cou[0][1][7) 
          _cou[3][1][6) = _wmp[2]*_cou[0][1][7) 
          _cou[0][5][6) = _wmq[2]*_cou[0][1][7) 
          _cou[1][2][6) = _wmp[0]*_cou[0][2][7) 
          _cou[2][2][6) = _wmp[1]*_cou[0][2][7) +  fakac * _cou[0][0][7)
          _cou[0][8][6) = _wmq[1]*_cou[0][2][7)  & 
               +  fakc * (_cou[0][0][6) - fakaca * _cou[0][0][7) ) 
          _cou[3][2][6) = _wmp[2]*_cou[0][2][7) 
          _cou[0][6][6) = _wmq[2]*_cou[0][2][7) 
          _cou[1][3][6) = _wmp[0]*_cou[0][3][7) 
          _cou[2][3][6) = _wmp[1]*_cou[0][3][7) 
          _cou[3][3][6) = _wmp[2]*_cou[0][3][7) +  fakac * _cou[0][0][7)
          _cou[0][9][6) = _wmq[2]*_cou[0][3][7)  & 
               +  fakc * (_cou[0][0][6) - fakaca * _cou[0][0][7) ) 
          _cou[7][0][6) = _wmp[0]*_cou[1][0][7)  & 
               +  faka * (_cou[0][0][6) - fakaac * _cou[0][0][7) ) 
          _cou[4][0][6) = _wmp[1]*_cou[1][0][7) 
          _cou[5][0][6) = _wmp[2]*_cou[1][0][7) 
          _cou[8][0][6) = _wmp[1]*_cou[2][0][7)  & 
               +  faka * (_cou[0][0][6) - fakaac * _cou[0][0][7) ) 
          _cou[6][0][6) = _wmp[2]*_cou[2][0][7) 
          _cou[9][0][6) = _wmp[2]*_cou[3][0][7)  & 
               +  faka * (_cou[0][0][6) - fakaac * _cou[0][0][7) ) 
          _cou[1][4][5) = _wmp[0]*_cou[0][4][6) +  fakac * _cou[0][2][6)
          _cou[0][13][5) = _wmq[0]*_cou[0][4][6)  & 
               +  fakc * (_cou[0][2][5) - fakaca * _cou[0][2][6) ) 
          _cou[2][4][5) = _wmp[1]*_cou[0][4][6) +  fakac * _cou[0][1][6)
          _cou[0][14][5) = _wmq[1]*_cou[0][4][6)  & 
               +  fakc * (_cou[0][1][5) - fakaca * _cou[0][1][6) ) 
          _cou[3][4][5) = _wmp[2]*_cou[0][4][6) 
          _cou[0][19][5) = _wmq[2]*_cou[0][4][6) 
          _cou[1][5][5) = _wmp[0]*_cou[0][5][6) +  fakac * _cou[0][3][6)
          _cou[0][15][5) = _wmq[0]*_cou[0][5][6)  & 
               +  fakc * (_cou[0][3][5) - fakaca * _cou[0][3][6) ) 
          _cou[2][5][5) = _wmp[1]*_cou[0][5][6) 
          _cou[3][5][5) = _wmp[2]*_cou[0][5][6) +  fakac * _cou[0][1][6)
          _cou[0][16][5) = _wmq[2]*_cou[0][5][6)  & 
               +  fakc * (_cou[0][1][5) - fakaca * _cou[0][1][6) ) 
          _cou[1][6][5) = _wmp[0]*_cou[0][6][6) 
          _cou[2][6][5) = _wmp[1]*_cou[0][6][6) +  fakac * _cou[0][3][6)
          _cou[0][17][5) = _wmq[1]*_cou[0][6][6)  & 
               +  fakc * (_cou[0][3][5) - fakaca * _cou[0][3][6) ) 
          _cou[3][6][5) = _wmp[2]*_cou[0][6][6) +  fakac * _cou[0][2][6)
          _cou[0][18][5) = _wmq[2]*_cou[0][6][6)  & 
               +  fakc * (_cou[0][2][5) - fakaca * _cou[0][2][6) ) 
          _cou[1][7][5) = _wmp[0]*_cou[0][7][6) +  fakac2 * _cou[0][1][6)
          _cou[0][10][5) = _wmq[0]*_cou[0][7][6)  & 
               +  fakc2 * (_cou[0][1][5) - fakaca * _cou[0][1][6) ) 
          _cou[2][7][5) = _wmp[1]*_cou[0][7][6) 
          _cou[3][7][5) = _wmp[2]*_cou[0][7][6) 
          _cou[1][8][5) = _wmp[0]*_cou[0][8][6) 
          _cou[2][8][5) = _wmp[1]*_cou[0][8][6) +  fakac2 * _cou[0][2][6)
          _cou[0][11][5) = _wmq[1]*_cou[0][8][6)  & 
               +  fakc2 * (_cou[0][2][5) - fakaca * _cou[0][2][6) ) 
          _cou[3][8][5) = _wmp[2]*_cou[0][8][6) 
          _cou[1][9][5) = _wmp[0]*_cou[0][9][6) 
          _cou[2][9][5) = _wmp[1]*_cou[0][9][6) 
          _cou[3][9][5) = _wmp[2]*_cou[0][9][6) +  fakac2 * _cou[0][3][6)
          _cou[0][12][5) = _wmq[2]*_cou[0][9][6)  & 
               +  fakc2 * (_cou[0][3][5) - fakaca * _cou[0][3][6) ) 
          _cou[7][1][5) = _wmp[0]*_cou[1][1][6)  & 
               +  faka * (_cou[0][1][5) - fakaac * _cou[0][1][6) ) + fakac * _cou[1][0][6)
          _cou[4][1][5) = _wmp[1]*_cou[1][1][6) 
          _cou[5][1][5) = _wmp[2]*_cou[1][1][6) 
          _cou[7][2][5) = _wmp[0]*_cou[1][2][6)  & 
               +  faka * (_cou[0][2][5) - fakaac * _cou[0][2][6) ) 
          _cou[4][2][5) = _wmp[1]*_cou[1][2][6) +  fakac * _cou[1][0][6)
          _cou[5][2][5) = _wmp[2]*_cou[1][2][6) 
          _cou[7][3][5) = _wmp[0]*_cou[1][3][6)  & 
               +  faka * (_cou[0][3][5) - fakaac * _cou[0][3][6) ) 
          _cou[4][3][5) = _wmp[1]*_cou[1][3][6) 
          _cou[5][3][5) = _wmp[2]*_cou[1][3][6) +  fakac * _cou[1][0][6)
          _cou[8][1][5) = _wmp[1]*_cou[2][1][6)  & 
               +  faka * (_cou[0][1][5) - fakaac * _cou[0][1][6) ) 
          _cou[6][1][5) = _wmp[2]*_cou[2][1][6) 
          _cou[8][2][5) = _wmp[1]*_cou[2][2][6)  & 
               +  faka * (_cou[0][2][5) - fakaac * _cou[0][2][6) ) + fakac * _cou[2][0][6)
          _cou[6][2][5) = _wmp[2]*_cou[2][2][6) 
          _cou[8][3][5) = _wmp[1]*_cou[2][3][6)  & 
               +  faka * (_cou[0][3][5) - fakaac * _cou[0][3][6) ) 
          _cou[6][3][5) = _wmp[2]*_cou[2][3][6) +  fakac * _cou[2][0][6)
          _cou[9][1][5) = _wmp[2]*_cou[3][1][6)  & 
               +  faka * (_cou[0][1][5) - fakaac * _cou[0][1][6) ) 
          _cou[9][2][5) = _wmp[2]*_cou[3][2][6)  & 
               +  faka * (_cou[0][2][5) - fakaac * _cou[0][2][6) ) 
          _cou[9][3][5) = _wmp[2]*_cou[3][3][6)  & 
               +  faka * (_cou[0][3][5) - fakaac * _cou[0][3][6) ) + fakac * _cou[3][0][6)
          _cou[13][0][5) = _wmp[0]*_cou[4][0][6)  & 
               +  faka * (_cou[2][0][5) - fakaac * _cou[2][0][6) ) 
          _cou[14][0][5) = _wmp[1]*_cou[4][0][6)  & 
               +  faka * (_cou[1][0][5) - fakaac * _cou[1][0][6) ) 
          _cou[19][0][5) = _wmp[2]*_cou[4][0][6) 
          _cou[15][0][5) = _wmp[0]*_cou[5][0][6)  & 
               +  faka * (_cou[3][0][5) - fakaac * _cou[3][0][6) ) 
          _cou[16][0][5) = _wmp[2]*_cou[5][0][6)  & 
               +  faka * (_cou[1][0][5) - fakaac * _cou[1][0][6) ) 
          _cou[17][0][5) = _wmp[1]*_cou[6][0][6)  & 
               +  faka * (_cou[3][0][5) - fakaac * _cou[3][0][6) ) 
          _cou[18][0][5) = _wmp[2]*_cou[6][0][6)  & 
               +  faka * (_cou[2][0][5) - fakaac * _cou[2][0][6) ) 
          _cou[10][0][5) = _wmp[0]*_cou[7][0][6)  & 
               +  faka2 * (_cou[1][0][5) - fakaac * _cou[1][0][6) ) 
          _cou[11][0][5) = _wmp[1]*_cou[8][0][6)  & 
               +  faka2 * (_cou[2][0][5) - fakaac * _cou[2][0][6) ) 
          _cou[12][0][5) = _wmp[2]*_cou[9][0][6)  & 
               +  faka2 * (_cou[3][0][5) - fakaac * _cou[3][0][6) ) 
          _cou[1][10][4) = _wmp[0]*_cou[0][10][5) +  fakac3 * _cou[0][7][5)
          _cou[0][20][4) = _wmq[0]*_cou[0][10][5)  & 
               +  fakc3 * (_cou[0][7][4) - fakaca * _cou[0][7][5) ) 
          _cou[2][10][4) = _wmp[1]*_cou[0][10][5) 
          _cou[0][23][4) = _wmq[1]*_cou[0][10][5) 
          _cou[3][10][4) = _wmp[2]*_cou[0][10][5) 
          _cou[0][25][4) = _wmq[2]*_cou[0][10][5) 
          _cou[1][11][4) = _wmp[0]*_cou[0][11][5) 
          _cou[0][24][4) = _wmq[0]*_cou[0][11][5) 
          _cou[2][11][4) = _wmp[1]*_cou[0][11][5) +  fakac3 * _cou[0][8][5)
          _cou[0][21][4) = _wmq[1]*_cou[0][11][5)  & 
               +  fakc3 * (_cou[0][8][4) - fakaca * _cou[0][8][5) ) 
          _cou[3][11][4) = _wmp[2]*_cou[0][11][5) 
          _cou[0][27][4) = _wmq[2]*_cou[0][11][5) 
          _cou[1][12][4) = _wmp[0]*_cou[0][12][5) 
          _cou[0][26][4) = _wmq[0]*_cou[0][12][5) 
          _cou[2][12][4) = _wmp[1]*_cou[0][12][5) 
          _cou[0][28][4) = _wmq[1]*_cou[0][12][5) 
          _cou[3][12][4) = _wmp[2]*_cou[0][12][5) +  fakac3 * _cou[0][9][5)
          _cou[0][22][4) = _wmq[2]*_cou[0][12][5)  & 
               +  fakc3 * (_cou[0][9][4) - fakaca * _cou[0][9][5) ) 
          _cou[1][13][4) = _wmp[0]*_cou[0][13][5) +  fakac2 * _cou[0][4][5)
          _cou[2][13][4) = _wmp[1]*_cou[0][13][5) +  fakac * _cou[0][7][5)
          _cou[0][31][4) = _wmq[1]*_cou[0][13][5)  & 
               +  fakc * (_cou[0][7][4) - fakaca * _cou[0][7][5) ) 
          _cou[3][13][4) = _wmp[2]*_cou[0][13][5) 
          _cou[0][32][4) = _wmq[2]*_cou[0][13][5) 
          _cou[1][14][4) = _wmp[0]*_cou[0][14][5) +  fakac * _cou[0][8][5)
          _cou[2][14][4) = _wmp[1]*_cou[0][14][5) +  fakac2 * _cou[0][4][5)
          _cou[3][14][4) = _wmp[2]*_cou[0][14][5) 
          _cou[0][33][4) = _wmq[2]*_cou[0][14][5) 
          _cou[1][15][4) = _wmp[0]*_cou[0][15][5) +  fakac2 * _cou[0][5][5)
          _cou[2][15][4) = _wmp[1]*_cou[0][15][5) 
          _cou[3][15][4) = _wmp[2]*_cou[0][15][5) +  fakac * _cou[0][7][5)
          _cou[0][29][4) = _wmq[2]*_cou[0][15][5)  & 
               +  fakc * (_cou[0][7][4) - fakaca * _cou[0][7][5) ) 
          _cou[1][16][4) = _wmp[0]*_cou[0][16][5) +  fakac * _cou[0][9][5)
          _cou[2][16][4) = _wmp[1]*_cou[0][16][5) 
          _cou[0][34][4) = _wmq[1]*_cou[0][16][5) 
          _cou[3][16][4) = _wmp[2]*_cou[0][16][5) +  fakac2 * _cou[0][5][5)
          _cou[1][17][4) = _wmp[0]*_cou[0][17][5) 
          _cou[2][17][4) = _wmp[1]*_cou[0][17][5) +  fakac2 * _cou[0][6][5)
          _cou[3][17][4) = _wmp[2]*_cou[0][17][5) +  fakac * _cou[0][8][5)
          _cou[0][30][4) = _wmq[2]*_cou[0][17][5)  & 
               +  fakc * (_cou[0][8][4) - fakaca * _cou[0][8][5) ) 
          _cou[1][18][4) = _wmp[0]*_cou[0][18][5) 
          _cou[2][18][4) = _wmp[1]*_cou[0][18][5) +  fakac * _cou[0][9][5)
          _cou[3][18][4) = _wmp[2]*_cou[0][18][5) +  fakac2 * _cou[0][6][5)
          _cou[1][19][4) = _wmp[0]*_cou[0][19][5) +  fakac * _cou[0][6][5)
          _cou[2][19][4) = _wmp[1]*_cou[0][19][5) +  fakac * _cou[0][5][5)
          _cou[3][19][4) = _wmp[2]*_cou[0][19][5) +  fakac * _cou[0][4][5)
          _cou[7][4][4) = _wmp[0]*_cou[1][4][5)  & 
               +  faka * (_cou[0][4][4) - fakaac * _cou[0][4][5) ) + fakac * _cou[1][2][5)
          _cou[4][4][4) = _wmp[1]*_cou[1][4][5) +  fakac * _cou[1][1][5)
          _cou[5][4][4) = _wmp[2]*_cou[1][4][5) 
          _cou[7][5][4) = _wmp[0]*_cou[1][5][5)  & 
               +  faka * (_cou[0][5][4) - fakaac * _cou[0][5][5) ) + fakac * _cou[1][3][5)
          _cou[4][5][4) = _wmp[1]*_cou[1][5][5) 
          _cou[5][5][4) = _wmp[2]*_cou[1][5][5) +  fakac * _cou[1][1][5)
          _cou[7][6][4) = _wmp[0]*_cou[1][6][5)  & 
               +  faka * (_cou[0][6][4) - fakaac * _cou[0][6][5) ) 
          _cou[4][6][4) = _wmp[1]*_cou[1][6][5) +  fakac * _cou[1][3][5)
          _cou[5][6][4) = _wmp[2]*_cou[1][6][5) +  fakac * _cou[1][2][5)
          _cou[7][7][4) = _wmp[0]*_cou[1][7][5)  & 
               +  faka * (_cou[0][7][4) - fakaac * _cou[0][7][5) ) + fakac2 * _cou[1][1][5)
          _cou[4][7][4) = _wmp[1]*_cou[1][7][5) 
          _cou[5][7][4) = _wmp[2]*_cou[1][7][5) 
          _cou[7][8][4) = _wmp[0]*_cou[1][8][5)  & 
               +  faka * (_cou[0][8][4) - fakaac * _cou[0][8][5) ) 
          _cou[4][8][4) = _wmp[1]*_cou[1][8][5) +  fakac2 * _cou[1][2][5)
          _cou[5][8][4) = _wmp[2]*_cou[1][8][5) 
          _cou[7][9][4) = _wmp[0]*_cou[1][9][5)  & 
               +  faka * (_cou[0][9][4) - fakaac * _cou[0][9][5) ) 
          _cou[4][9][4) = _wmp[1]*_cou[1][9][5) 
          _cou[5][9][4) = _wmp[2]*_cou[1][9][5) +  fakac2 * _cou[1][3][5)
          _cou[8][4][4) = _wmp[1]*_cou[2][4][5)  & 
               +  faka * (_cou[0][4][4) - fakaac * _cou[0][4][5) ) + fakac * _cou[2][1][5)
          _cou[6][4][4) = _wmp[2]*_cou[2][4][5) 
          _cou[8][5][4) = _wmp[1]*_cou[2][5][5)  & 
               +  faka * (_cou[0][5][4) - fakaac * _cou[0][5][5) ) 
          _cou[6][5][4) = _wmp[2]*_cou[2][5][5) +  fakac * _cou[2][1][5)
          _cou[8][6][4) = _wmp[1]*_cou[2][6][5)  & 
               +  faka * (_cou[0][6][4) - fakaac * _cou[0][6][5) ) + fakac * _cou[2][3][5)
          _cou[6][6][4) = _wmp[2]*_cou[2][6][5) +  fakac * _cou[2][2][5)
          _cou[8][7][4) = _wmp[1]*_cou[2][7][5)  & 
               +  faka * (_cou[0][7][4) - fakaac * _cou[0][7][5) ) 
          _cou[6][7][4) = _wmp[2]*_cou[2][7][5) 
          _cou[8][8][4) = _wmp[1]*_cou[2][8][5)  & 
               +  faka * (_cou[0][8][4) - fakaac * _cou[0][8][5) ) + fakac2 * _cou[2][2][5)
          _cou[6][8][4) = _wmp[2]*_cou[2][8][5) 
          _cou[8][9][4) = _wmp[1]*_cou[2][9][5)  & 
               +  faka * (_cou[0][9][4) - fakaac * _cou[0][9][5) ) 
          _cou[6][9][4) = _wmp[2]*_cou[2][9][5) +  fakac2 * _cou[2][3][5)
          _cou[9][4][4) = _wmp[2]*_cou[3][4][5)  & 
               +  faka * (_cou[0][4][4) - fakaac * _cou[0][4][5) ) 
          _cou[9][5][4) = _wmp[2]*_cou[3][5][5)  & 
               +  faka * (_cou[0][5][4) - fakaac * _cou[0][5][5) ) + fakac * _cou[3][1][5)
          _cou[9][6][4) = _wmp[2]*_cou[3][6][5)  & 
               +  faka * (_cou[0][6][4) - fakaac * _cou[0][6][5) ) + fakac * _cou[3][2][5)
          _cou[9][7][4) = _wmp[2]*_cou[3][7][5)  & 
               +  faka * (_cou[0][7][4) - fakaac * _cou[0][7][5) ) 
          _cou[9][8][4) = _wmp[2]*_cou[3][8][5)  & 
               +  faka * (_cou[0][8][4) - fakaac * _cou[0][8][5) ) 
          _cou[9][9][4) = _wmp[2]*_cou[3][9][5)  & 
               +  faka * (_cou[0][9][4) - fakaac * _cou[0][9][5) ) + fakac2 * _cou[3][3][5)
          _cou[13][1][4) = _wmp[0]*_cou[4][1][5)  & 
               +  faka * (_cou[2][1][4) - fakaac * _cou[2][1][5) ) + fakac * _cou[4][0][5)
          _cou[14][1][4) = _wmp[1]*_cou[4][1][5)  & 
               +  faka * (_cou[1][1][4) - fakaac * _cou[1][1][5) ) 
          _cou[19][1][4) = _wmp[2]*_cou[4][1][5) 
          _cou[13][2][4) = _wmp[0]*_cou[4][2][5)  & 
               +  faka * (_cou[2][2][4) - fakaac * _cou[2][2][5) ) 
          _cou[14][2][4) = _wmp[1]*_cou[4][2][5)  & 
               +  faka * (_cou[1][2][4) - fakaac * _cou[1][2][5) ) + fakac * _cou[4][0][5)
          _cou[19][2][4) = _wmp[2]*_cou[4][2][5) 
          _cou[13][3][4) = _wmp[0]*_cou[4][3][5)  & 
               +  faka * (_cou[2][3][4) - fakaac * _cou[2][3][5) ) 
          _cou[14][3][4) = _wmp[1]*_cou[4][3][5)  & 
               +  faka * (_cou[1][3][4) - fakaac * _cou[1][3][5) ) 
          _cou[19][3][4) = _wmp[2]*_cou[4][3][5) +  fakac * _cou[4][0][5)
          _cou[15][1][4) = _wmp[0]*_cou[5][1][5)  & 
               +  faka * (_cou[3][1][4) - fakaac * _cou[3][1][5) ) + fakac * _cou[5][0][5)
          _cou[16][1][4) = _wmp[2]*_cou[5][1][5)  & 
               +  faka * (_cou[1][1][4) - fakaac * _cou[1][1][5) ) 
          _cou[15][2][4) = _wmp[0]*_cou[5][2][5)  & 
               +  faka * (_cou[3][2][4) - fakaac * _cou[3][2][5) ) 
          _cou[16][2][4) = _wmp[2]*_cou[5][2][5)  & 
               +  faka * (_cou[1][2][4) - fakaac * _cou[1][2][5) ) 
          _cou[15][3][4) = _wmp[0]*_cou[5][3][5)  & 
               +  faka * (_cou[3][3][4) - fakaac * _cou[3][3][5) ) 
          _cou[16][3][4) = _wmp[2]*_cou[5][3][5)  & 
               +  faka * (_cou[1][3][4) - fakaac * _cou[1][3][5) ) + fakac * _cou[5][0][5)
          _cou[17][1][4) = _wmp[1]*_cou[6][1][5)  & 
               +  faka * (_cou[3][1][4) - fakaac * _cou[3][1][5) ) 
          _cou[18][1][4) = _wmp[2]*_cou[6][1][5)  & 
               +  faka * (_cou[2][1][4) - fakaac * _cou[2][1][5) ) 
          _cou[17][2][4) = _wmp[1]*_cou[6][2][5)  & 
               +  faka * (_cou[3][2][4) - fakaac * _cou[3][2][5) ) + fakac * _cou[6][0][5)
          _cou[18][2][4) = _wmp[2]*_cou[6][2][5)  & 
               +  faka * (_cou[2][2][4) - fakaac * _cou[2][2][5) ) 
          _cou[17][3][4) = _wmp[1]*_cou[6][3][5)  & 
               +  faka * (_cou[3][3][4) - fakaac * _cou[3][3][5) ) 
          _cou[18][3][4) = _wmp[2]*_cou[6][3][5)  & 
               +  faka * (_cou[2][3][4) - fakaac * _cou[2][3][5) ) + fakac * _cou[6][0][5)
          _cou[10][1][4) = _wmp[0]*_cou[7][1][5)  & 
               +  faka2 * (_cou[1][1][4) - fakaac * _cou[1][1][5) ) + fakac * _cou[7][0][5)
          _cou[10][2][4) = _wmp[0]*_cou[7][2][5)  & 
               +  faka2 * (_cou[1][2][4) - fakaac * _cou[1][2][5) ) 
          _cou[10][3][4) = _wmp[0]*_cou[7][3][5)  & 
               +  faka2 * (_cou[1][3][4) - fakaac * _cou[1][3][5) ) 
          _cou[11][1][4) = _wmp[1]*_cou[8][1][5)  & 
               +  faka2 * (_cou[2][1][4) - fakaac * _cou[2][1][5) ) 
          _cou[11][2][4) = _wmp[1]*_cou[8][2][5)  & 
               +  faka2 * (_cou[2][2][4) - fakaac * _cou[2][2][5) ) + fakac * _cou[8][0][5)
          _cou[11][3][4) = _wmp[1]*_cou[8][3][5)  & 
               +  faka2 * (_cou[2][3][4) - fakaac * _cou[2][3][5) ) 
          _cou[12][1][4) = _wmp[2]*_cou[9][1][5)  & 
               +  faka2 * (_cou[3][1][4) - fakaac * _cou[3][1][5) ) 
          _cou[12][2][4) = _wmp[2]*_cou[9][2][5)  & 
               +  faka2 * (_cou[3][2][4) - fakaac * _cou[3][2][5) ) 
          _cou[12][3][4) = _wmp[2]*_cou[9][3][5)  & 
               +  faka2 * (_cou[3][3][4) - fakaac * _cou[3][3][5) ) + fakac * _cou[9][0][5)


          _cou[20][0][4) = _wmp[0]*_cou[10][0][5)  & 
               +  faka3 * (_cou[7][0][4) - fakaac * _cou[7][0][5) ) 
          _cou[23][0][4) = _wmp[1]*_cou[10][0][5) 
          _cou[25][0][4) = _wmp[2]*_cou[10][0][5) 
          _cou[24][0][4) = _wmp[0]*_cou[11][0][5) 
          _cou[21][0][4) = _wmp[1]*_cou[11][0][5)  & 
               +  faka3 * (_cou[8][0][4) - fakaac * _cou[8][0][5) ) 
          _cou[27][0][4) = _wmp[2]*_cou[11][0][5) 
          _cou[26][0][4) = _wmp[0]*_cou[12][0][5) 
          _cou[28][0][4) = _wmp[1]*_cou[12][0][5) 
          _cou[22][0][4) = _wmp[2]*_cou[12][0][5)  & 
               +  faka3 * (_cou[9][0][4) - fakaac * _cou[9][0][5) ) 
          _cou[31][0][4) = _wmp[1]*_cou[13][0][5)  & 
               +  faka * (_cou[7][0][4) - fakaac * _cou[7][0][5) ) 
          _cou[32][0][4) = _wmp[2]*_cou[13][0][5) 
          _cou[33][0][4) = _wmp[2]*_cou[14][0][5) 
          _cou[29][0][4) = _wmp[2]*_cou[15][0][5)  & 
               +  faka * (_cou[7][0][4) - fakaac * _cou[7][0][5) ) 
          _cou[34][0][4) = _wmp[1]*_cou[16][0][5) 
          _cou[30][0][4) = _wmp[2]*_cou[17][0][5)  & 
               +  faka * (_cou[8][0][4) - fakaac * _cou[8][0][5) ) 

          _cou[1][20][3) = _wmp[0]*_cou[0][20][4) +  fakac4 * _cou[0][10][4)
          _cou[2][20][3) = _wmp[1]*_cou[0][20][4) 
          _cou[3][20][3) = _wmp[2]*_cou[0][20][4) 
          _cou[1][21][3) = _wmp[0]*_cou[0][21][4) 
          _cou[2][21][3) = _wmp[1]*_cou[0][21][4) +  fakac4 * _cou[0][11][4)
          _cou[3][21][3) = _wmp[2]*_cou[0][21][4) 
          _cou[1][22][3) = _wmp[0]*_cou[0][22][4) 
          _cou[2][22][3) = _wmp[1]*_cou[0][22][4) 
          _cou[3][22][3) = _wmp[2]*_cou[0][22][4) +  fakac4 * _cou[0][12][4)
          _cou[1][23][3) = _wmp[0]*_cou[0][23][4) +  fakac3 * _cou[0][13][4)
          _cou[2][23][3) = _wmp[1]*_cou[0][23][4) +  fakac * _cou[0][10][4)
          _cou[3][23][3) = _wmp[2]*_cou[0][23][4) 
          _cou[1][24][3) = _wmp[0]*_cou[0][24][4) +  fakac * _cou[0][11][4)
          _cou[2][24][3) = _wmp[1]*_cou[0][24][4) +  fakac3 * _cou[0][14][4)
          _cou[3][24][3) = _wmp[2]*_cou[0][24][4) 
          _cou[1][25][3) = _wmp[0]*_cou[0][25][4) +  fakac3 * _cou[0][15][4)
          _cou[2][25][3) = _wmp[1]*_cou[0][25][4) 
          _cou[3][25][3) = _wmp[2]*_cou[0][25][4) +  fakac * _cou[0][10][4)
          _cou[1][26][3) = _wmp[0]*_cou[0][26][4) +  fakac * _cou[0][12][4)
          _cou[2][26][3) = _wmp[1]*_cou[0][26][4) 
          _cou[3][26][3) = _wmp[2]*_cou[0][26][4) +  fakac3 * _cou[0][16][4)
          _cou[1][27][3) = _wmp[0]*_cou[0][27][4) 
          _cou[2][27][3) = _wmp[1]*_cou[0][27][4) +  fakac3 * _cou[0][17][4)
          _cou[3][27][3) = _wmp[2]*_cou[0][27][4) +  fakac * _cou[0][11][4)
          _cou[1][28][3) = _wmp[0]*_cou[0][28][4) 
          _cou[2][28][3) = _wmp[1]*_cou[0][28][4) +  fakac * _cou[0][12][4)
          _cou[3][28][3) = _wmp[2]*_cou[0][28][4) +  fakac3 * _cou[0][18][4)
          _cou[1][29][3) = _wmp[0]*_cou[0][29][4) +  fakac2 * _cou[0][16][4)
          _cou[2][29][3) = _wmp[1]*_cou[0][29][4) 
          _cou[3][29][3) = _wmp[2]*_cou[0][29][4) +  fakac2 * _cou[0][15][4)
          _cou[1][30][3) = _wmp[0]*_cou[0][30][4) 
          _cou[2][30][3) = _wmp[1]*_cou[0][30][4) +  fakac2 * _cou[0][18][4)
          _cou[3][30][3) = _wmp[2]*_cou[0][30][4) +  fakac2 * _cou[0][17][4)
          _cou[1][31][3) = _wmp[0]*_cou[0][31][4) +  fakac2 * _cou[0][14][4)
          _cou[2][31][3) = _wmp[1]*_cou[0][31][4) +  fakac2 * _cou[0][13][4)
          _cou[3][31][3) = _wmp[2]*_cou[0][31][4) 
          _cou[1][32][3) = _wmp[0]*_cou[0][32][4) +  fakac2 * _cou[0][19][4)
          _cou[2][32][3) = _wmp[1]*_cou[0][32][4) +  fakac * _cou[0][15][4)
          _cou[3][32][3) = _wmp[2]*_cou[0][32][4) +  fakac * _cou[0][13][4)
          _cou[1][33][3) = _wmp[0]*_cou[0][33][4) +  fakac * _cou[0][17][4)
          _cou[2][33][3) = _wmp[1]*_cou[0][33][4) +  fakac2 * _cou[0][19][4)
          _cou[3][33][3) = _wmp[2]*_cou[0][33][4) +  fakac * _cou[0][14][4)
          _cou[1][34][3) = _wmp[0]*_cou[0][34][4) +  fakac * _cou[0][18][4)
          _cou[2][34][3) = _wmp[1]*_cou[0][34][4) +  fakac * _cou[0][16][4)
          _cou[3][34][3) = _wmp[2]*_cou[0][34][4) +  fakac2 * _cou[0][19][4)


          _cou[7][10][3) = _wmp[0]*_cou[1][10][4)  & 
               +  faka * (_cou[0][10][3) - fakaac * _cou[0][10][4) ) + fakac3 * _cou[1][7][4)
          _cou[4][10][3) = _wmp[1]*_cou[1][10][4) 
          _cou[5][10][3) = _wmp[2]*_cou[1][10][4) 
          _cou[7][11][3) = _wmp[0]*_cou[1][11][4)  & 
               +  faka * (_cou[0][11][3) - fakaac * _cou[0][11][4) ) 
          _cou[4][11][3) = _wmp[1]*_cou[1][11][4) +  fakac3 * _cou[1][8][4)
          _cou[5][11][3) = _wmp[2]*_cou[1][11][4) 
          _cou[7][12][3) = _wmp[0]*_cou[1][12][4)  & 
               +  faka * (_cou[0][12][3) - fakaac * _cou[0][12][4) ) 
          _cou[4][12][3) = _wmp[1]*_cou[1][12][4) 
          _cou[5][12][3) = _wmp[2]*_cou[1][12][4) +  fakac3 * _cou[1][9][4)
          _cou[7][13][3) = _wmp[0]*_cou[1][13][4)  & 
               +  faka * (_cou[0][13][3) - fakaac * _cou[0][13][4) ) + fakac2 * _cou[1][4][4)
          _cou[4][13][3) = _wmp[1]*_cou[1][13][4) +  fakac * _cou[1][7][4)
          _cou[5][13][3) = _wmp[2]*_cou[1][13][4) 
          _cou[7][14][3) = _wmp[0]*_cou[1][14][4)  & 
               +  faka * (_cou[0][14][3) - fakaac * _cou[0][14][4) ) + fakac * _cou[1][8][4)
          _cou[4][14][3) = _wmp[1]*_cou[1][14][4) +  fakac2 * _cou[1][4][4)
          _cou[5][14][3) = _wmp[2]*_cou[1][14][4) 
          _cou[7][15][3) = _wmp[0]*_cou[1][15][4)  & 
               +  faka * (_cou[0][15][3) - fakaac * _cou[0][15][4) ) + fakac2 * _cou[1][5][4)
          _cou[4][15][3) = _wmp[1]*_cou[1][15][4) 
          _cou[5][15][3) = _wmp[2]*_cou[1][15][4) +  fakac * _cou[1][7][4)
          _cou[7][16][3) = _wmp[0]*_cou[1][16][4)  & 
               +  faka * (_cou[0][16][3) - fakaac * _cou[0][16][4) ) + fakac * _cou[1][9][4)
          _cou[4][16][3) = _wmp[1]*_cou[1][16][4) 
          _cou[5][16][3) = _wmp[2]*_cou[1][16][4) +  fakac2 * _cou[1][5][4)
          _cou[7][17][3) = _wmp[0]*_cou[1][17][4)  & 
               +  faka * (_cou[0][17][3) - fakaac * _cou[0][17][4) ) 
          _cou[4][17][3) = _wmp[1]*_cou[1][17][4) +  fakac2 * _cou[1][6][4)
          _cou[5][17][3) = _wmp[2]*_cou[1][17][4) +  fakac * _cou[1][8][4)
          _cou[7][18][3) = _wmp[0]*_cou[1][18][4)  & 
               +  faka * (_cou[0][18][3) - fakaac * _cou[0][18][4) ) 
          _cou[4][18][3) = _wmp[1]*_cou[1][18][4) +  fakac * _cou[1][9][4)
          _cou[5][18][3) = _wmp[2]*_cou[1][18][4) +  fakac2 * _cou[1][6][4)
          _cou[7][19][3) = _wmp[0]*_cou[1][19][4)  & 
               +  faka * (_cou[0][19][3) - fakaac * _cou[0][19][4) ) + fakac * _cou[1][6][4)
          _cou[4][19][3) = _wmp[1]*_cou[1][19][4) +  fakac * _cou[1][5][4)
          _cou[5][19][3) = _wmp[2]*_cou[1][19][4) +  fakac * _cou[1][4][4)
          _cou[8][10][3) = _wmp[1]*_cou[2][10][4)  & 
               +  faka * (_cou[0][10][3) - fakaac * _cou[0][10][4) ) 
          _cou[6][10][3) = _wmp[2]*_cou[2][10][4) 
          _cou[8][11][3) = _wmp[1]*_cou[2][11][4)  & 
               +  faka * (_cou[0][11][3) - fakaac * _cou[0][11][4) ) + fakac3 * _cou[2][8][4)
          _cou[6][11][3) = _wmp[2]*_cou[2][11][4) 
          _cou[8][12][3) = _wmp[1]*_cou[2][12][4)  & 
               +  faka * (_cou[0][12][3) - fakaac * _cou[0][12][4) ) 
          _cou[6][12][3) = _wmp[2]*_cou[2][12][4) +  fakac3 * _cou[2][9][4)
          _cou[8][13][3) = _wmp[1]*_cou[2][13][4)  & 
               +  faka * (_cou[0][13][3) - fakaac * _cou[0][13][4) ) + fakac * _cou[2][7][4)
          _cou[6][13][3) = _wmp[2]*_cou[2][13][4) 
          _cou[8][14][3) = _wmp[1]*_cou[2][14][4)  & 
               +  faka * (_cou[0][14][3) - fakaac * _cou[0][14][4) ) + fakac2 * _cou[2][4][4)
          _cou[6][14][3) = _wmp[2]*_cou[2][14][4) 
          _cou[8][15][3) = _wmp[1]*_cou[2][15][4)  & 
               +  faka * (_cou[0][15][3) - fakaac * _cou[0][15][4) ) 
          _cou[6][15][3) = _wmp[2]*_cou[2][15][4) +  fakac * _cou[2][7][4)
          _cou[8][16][3) = _wmp[1]*_cou[2][16][4)  & 
               +  faka * (_cou[0][16][3) - fakaac * _cou[0][16][4) ) 
          _cou[6][16][3) = _wmp[2]*_cou[2][16][4) +  fakac2 * _cou[2][5][4)
          _cou[8][17][3) = _wmp[1]*_cou[2][17][4)  & 
               +  faka * (_cou[0][17][3) - fakaac * _cou[0][17][4) ) + fakac2 * _cou[2][6][4)
          _cou[6][17][3) = _wmp[2]*_cou[2][17][4) +  fakac * _cou[2][8][4)
          _cou[8][18][3) = _wmp[1]*_cou[2][18][4)  & 
               +  faka * (_cou[0][18][3) - fakaac * _cou[0][18][4) ) + fakac * _cou[2][9][4)
          _cou[6][18][3) = _wmp[2]*_cou[2][18][4) +  fakac2 * _cou[2][6][4)
          _cou[8][19][3) = _wmp[1]*_cou[2][19][4)  & 
               +  faka * (_cou[0][19][3) - fakaac * _cou[0][19][4) ) + fakac * _cou[2][5][4)
          _cou[6][19][3) = _wmp[2]*_cou[2][19][4) +  fakac * _cou[2][4][4)
          _cou[9][10][3) = _wmp[2]*_cou[3][10][4)  & 
               +  faka * (_cou[0][10][3) - fakaac * _cou[0][10][4) ) 
          _cou[9][11][3) = _wmp[2]*_cou[3][11][4)  & 
               +  faka * (_cou[0][11][3) - fakaac * _cou[0][11][4) ) 
          _cou[9][12][3) = _wmp[2]*_cou[3][12][4)  & 
               +  faka * (_cou[0][12][3) - fakaac * _cou[0][12][4) ) + fakac3 * _cou[3][9][4)
          _cou[9][13][3) = _wmp[2]*_cou[3][13][4)  & 
               +  faka * (_cou[0][13][3) - fakaac * _cou[0][13][4) ) 
          _cou[9][14][3) = _wmp[2]*_cou[3][14][4)  & 
               +  faka * (_cou[0][14][3) - fakaac * _cou[0][14][4) ) 
          _cou[9][15][3) = _wmp[2]*_cou[3][15][4)  & 
               +  faka * (_cou[0][15][3) - fakaac * _cou[0][15][4) ) + fakac * _cou[3][7][4)
          _cou[9][16][3) = _wmp[2]*_cou[3][16][4)  & 
               +  faka * (_cou[0][16][3) - fakaac * _cou[0][16][4) ) + fakac2 * _cou[3][5][4)
          _cou[9][17][3) = _wmp[2]*_cou[3][17][4)  & 
               +  faka * (_cou[0][17][3) - fakaac * _cou[0][17][4) ) + fakac * _cou[3][8][4)
          _cou[9][18][3) = _wmp[2]*_cou[3][18][4)  & 
               +  faka * (_cou[0][18][3) - fakaac * _cou[0][18][4) ) + fakac2 * _cou[3][6][4)
          _cou[9][19][3) = _wmp[2]*_cou[3][19][4)  & 
               +  faka * (_cou[0][19][3) - fakaac * _cou[0][19][4) ) + fakac * _cou[3][4][4)


          _cou[13][4][3) = _wmp[0]*_cou[4][4][4)  & 
               +  faka * (_cou[2][4][3) - fakaac * _cou[2][4][4) ) + fakac * _cou[4][2][4)
          _cou[14][4][3) = _wmp[1]*_cou[4][4][4)  & 
               +  faka * (_cou[1][4][3) - fakaac * _cou[1][4][4) ) + fakac * _cou[4][1][4)
          _cou[19][4][3) = _wmp[2]*_cou[4][4][4) 
          _cou[13][5][3) = _wmp[0]*_cou[4][5][4)  & 
               +  faka * (_cou[2][5][3) - fakaac * _cou[2][5][4) ) + fakac * _cou[4][3][4)
          _cou[14][5][3) = _wmp[1]*_cou[4][5][4)  & 
               +  faka * (_cou[1][5][3) - fakaac * _cou[1][5][4) ) 
          _cou[19][5][3) = _wmp[2]*_cou[4][5][4) +  fakac * _cou[4][1][4)
          _cou[13][6][3) = _wmp[0]*_cou[4][6][4)  & 
               +  faka * (_cou[2][6][3) - fakaac * _cou[2][6][4) ) 
          _cou[14][6][3) = _wmp[1]*_cou[4][6][4)  & 
               +  faka * (_cou[1][6][3) - fakaac * _cou[1][6][4) ) + fakac * _cou[4][3][4)
          _cou[19][6][3) = _wmp[2]*_cou[4][6][4) +  fakac * _cou[4][2][4)
          _cou[13][7][3) = _wmp[0]*_cou[4][7][4)  & 
               +  faka * (_cou[2][7][3) - fakaac * _cou[2][7][4) ) + fakac2 * _cou[4][1][4)
          _cou[14][7][3) = _wmp[1]*_cou[4][7][4)  & 
               +  faka * (_cou[1][7][3) - fakaac * _cou[1][7][4) ) 
          _cou[19][7][3) = _wmp[2]*_cou[4][7][4) 
          _cou[13][8][3) = _wmp[0]*_cou[4][8][4)  & 
               +  faka * (_cou[2][8][3) - fakaac * _cou[2][8][4) ) 
          _cou[14][8][3) = _wmp[1]*_cou[4][8][4)  & 
               +  faka * (_cou[1][8][3) - fakaac * _cou[1][8][4) ) + fakac2 * _cou[4][2][4)
          _cou[19][8][3) = _wmp[2]*_cou[4][8][4) 
          _cou[13][9][3) = _wmp[0]*_cou[4][9][4)  & 
               +  faka * (_cou[2][9][3) - fakaac * _cou[2][9][4) ) 
          _cou[14][9][3) = _wmp[1]*_cou[4][9][4)  & 
               +  faka * (_cou[1][9][3) - fakaac * _cou[1][9][4) ) 
          _cou[19][9][3) = _wmp[2]*_cou[4][9][4) +  fakac2 * _cou[4][3][4)
          _cou[15][4][3) = _wmp[0]*_cou[5][4][4)  & 
               +  faka * (_cou[3][4][3) - fakaac * _cou[3][4][4) ) + fakac * _cou[5][2][4)
          _cou[16][4][3) = _wmp[2]*_cou[5][4][4)  & 
               +  faka * (_cou[1][4][3) - fakaac * _cou[1][4][4) ) 
          _cou[15][5][3) = _wmp[0]*_cou[5][5][4)  & 
               +  faka * (_cou[3][5][3) - fakaac * _cou[3][5][4) ) + fakac * _cou[5][3][4)
          _cou[16][5][3) = _wmp[2]*_cou[5][5][4)  & 
               +  faka * (_cou[1][5][3) - fakaac * _cou[1][5][4) ) + fakac * _cou[5][1][4)
          _cou[15][6][3) = _wmp[0]*_cou[5][6][4)  & 
               +  faka * (_cou[3][6][3) - fakaac * _cou[3][6][4) ) 
          _cou[16][6][3) = _wmp[2]*_cou[5][6][4)  & 
               +  faka * (_cou[1][6][3) - fakaac * _cou[1][6][4) ) + fakac * _cou[5][2][4)
          _cou[15][7][3) = _wmp[0]*_cou[5][7][4)  & 
               +  faka * (_cou[3][7][3) - fakaac * _cou[3][7][4) ) + fakac2 * _cou[5][1][4)
          _cou[16][7][3) = _wmp[2]*_cou[5][7][4)  & 
               +  faka * (_cou[1][7][3) - fakaac * _cou[1][7][4) ) 
          _cou[15][8][3) = _wmp[0]*_cou[5][8][4)  & 
               +  faka * (_cou[3][8][3) - fakaac * _cou[3][8][4) ) 
          _cou[16][8][3) = _wmp[2]*_cou[5][8][4)  & 
               +  faka * (_cou[1][8][3) - fakaac * _cou[1][8][4) ) 
          _cou[15][9][3) = _wmp[0]*_cou[5][9][4)  & 
               +  faka * (_cou[3][9][3) - fakaac * _cou[3][9][4) ) 
          _cou[16][9][3) = _wmp[2]*_cou[5][9][4)  & 
               +  faka * (_cou[1][9][3) - fakaac * _cou[1][9][4) ) + fakac2 * _cou[5][3][4)
          _cou[17][4][3) = _wmp[1]*_cou[6][4][4)  & 
               +  faka * (_cou[3][4][3) - fakaac * _cou[3][4][4) ) + fakac * _cou[6][1][4)
          _cou[18][4][3) = _wmp[2]*_cou[6][4][4)  & 
               +  faka * (_cou[2][4][3) - fakaac * _cou[2][4][4) ) 
          _cou[17][5][3) = _wmp[1]*_cou[6][5][4)  & 
               +  faka * (_cou[3][5][3) - fakaac * _cou[3][5][4) ) 
          _cou[18][5][3) = _wmp[2]*_cou[6][5][4)  & 
               +  faka * (_cou[2][5][3) - fakaac * _cou[2][5][4) ) + fakac * _cou[6][1][4)
          _cou[17][6][3) = _wmp[1]*_cou[6][6][4)  & 
               +  faka * (_cou[3][6][3) - fakaac * _cou[3][6][4) ) + fakac * _cou[6][3][4)
          _cou[18][6][3) = _wmp[2]*_cou[6][6][4)  & 
               +  faka * (_cou[2][6][3) - fakaac * _cou[2][6][4) ) + fakac * _cou[6][2][4)
          _cou[17][7][3) = _wmp[1]*_cou[6][7][4)  & 
               +  faka * (_cou[3][7][3) - fakaac * _cou[3][7][4) ) 
          _cou[18][7][3) = _wmp[2]*_cou[6][7][4)  & 
               +  faka * (_cou[2][7][3) - fakaac * _cou[2][7][4) ) 
          _cou[17][8][3) = _wmp[1]*_cou[6][8][4)  & 
               +  faka * (_cou[3][8][3) - fakaac * _cou[3][8][4) ) + fakac2 * _cou[6][2][4)
          _cou[18][8][3) = _wmp[2]*_cou[6][8][4)  & 
               +  faka * (_cou[2][8][3) - fakaac * _cou[2][8][4) ) 
          _cou[17][9][3) = _wmp[1]*_cou[6][9][4)  & 
               +  faka * (_cou[3][9][3) - fakaac * _cou[3][9][4) ) 
          _cou[18][9][3) = _wmp[2]*_cou[6][9][4)  & 
               +  faka * (_cou[2][9][3) - fakaac * _cou[2][9][4) ) + fakac2 * _cou[6][3][4)
          _cou[10][4][3) = _wmp[0]*_cou[7][4][4)  & 
               +  faka2 * (_cou[1][4][3) - fakaac * _cou[1][4][4) ) + fakac * _cou[7][2][4)
          _cou[10][5][3) = _wmp[0]*_cou[7][5][4)  & 
               +  faka2 * (_cou[1][5][3) - fakaac * _cou[1][5][4) ) + fakac * _cou[7][3][4)
          _cou[10][6][3) = _wmp[0]*_cou[7][6][4)  & 
               +  faka2 * (_cou[1][6][3) - fakaac * _cou[1][6][4) ) 
          _cou[10][7][3) = _wmp[0]*_cou[7][7][4)  & 
               + faka2 * (_cou[1][7][3) - fakaac * _cou[1][7][4) ) + fakac2 * _cou[7][1][4)
          _cou[10][8][3) = _wmp[0]*_cou[7][8][4)  & 
               +  faka2 * (_cou[1][8][3) - fakaac * _cou[1][8][4) ) 
          _cou[10][9][3) = _wmp[0]*_cou[7][9][4)  & 
               +  faka2 * (_cou[1][9][3) - fakaac * _cou[1][9][4) ) 
          _cou[11][4][3) = _wmp[1]*_cou[8][4][4)  & 
               +  faka2 * (_cou[2][4][3) - fakaac * _cou[2][4][4) ) + fakac * _cou[8][1][4)
          _cou[11][5][3) = _wmp[1]*_cou[8][5][4)  & 
               +  faka2 * (_cou[2][5][3) - fakaac * _cou[2][5][4) ) 
          _cou[11][6][3) = _wmp[1]*_cou[8][6][4)  & 
               +  faka2 * (_cou[2][6][3) - fakaac * _cou[2][6][4) ) + fakac * _cou[8][3][4)
          _cou[11][7][3) = _wmp[1]*_cou[8][7][4)  & 
               +  faka2 * (_cou[2][7][3) - fakaac * _cou[2][7][4) ) 
          _cou[11][8][3) = _wmp[1]*_cou[8][8][4)  & 
               + faka2 * (_cou[2][8][3) - fakaac * _cou[2][8][4) ) + fakac2 * _cou[8][2][4)
          _cou[11][9][3) = _wmp[1]*_cou[8][9][4)  & 
               +  faka2 * (_cou[2][9][3) - fakaac * _cou[2][9][4) ) 
          _cou[12][4][3) = _wmp[2]*_cou[9][4][4)  & 
               +  faka2 * (_cou[3][4][3) - fakaac * _cou[3][4][4) ) 
          _cou[12][5][3) = _wmp[2]*_cou[9][5][4)  & 
               +  faka2 * (_cou[3][5][3) - fakaac * _cou[3][5][4) ) + fakac * _cou[9][1][4)
          _cou[12][6][3) = _wmp[2]*_cou[9][6][4)  & 
               +  faka2 * (_cou[3][6][3) - fakaac * _cou[3][6][4) ) + fakac * _cou[9][2][4)
          _cou[12][7][3) = _wmp[2]*_cou[9][7][4)  & 
               +  faka2 * (_cou[3][7][3) - fakaac * _cou[3][7][4) ) 
          _cou[12][8][3) = _wmp[2]*_cou[9][8][4)  & 
               +  faka2 * (_cou[3][8][3) - fakaac * _cou[3][8][4) ) 
          _cou[12][9][3) = _wmp[2]*_cou[9][9][4)  & 
               + faka2 * (_cou[3][9][3) - fakaac * _cou[3][9][4) ) + fakac2 * _cou[9][3][4)

          _cou[20][1][3) = _wmp[0]*_cou[10][1][4)  & 
               +  faka3 * (_cou[7][1][3) - fakaac * _cou[7][1][4) ) + fakac * _cou[10][0][4)
          _cou[23][1][3) = _wmp[1]*_cou[10][1][4) 
          _cou[25][1][3) = _wmp[2]*_cou[10][1][4) 
          _cou[20][2][3) = _wmp[0]*_cou[10][2][4)  & 
               +  faka3 * (_cou[7][2][3) - fakaac * _cou[7][2][4) ) 
          _cou[23][2][3) = _wmp[1]*_cou[10][2][4) +  fakac * _cou[10][0][4)
          _cou[25][2][3) = _wmp[2]*_cou[10][2][4) 
          _cou[20][3][3) = _wmp[0]*_cou[10][3][4)  & 
               +  faka3 * (_cou[7][3][3) - fakaac * _cou[7][3][4) ) 
          _cou[23][3][3) = _wmp[1]*_cou[10][3][4) 
          _cou[25][3][3) = _wmp[2]*_cou[10][3][4) +  fakac * _cou[10][0][4)
          _cou[24][1][3) = _wmp[0]*_cou[11][1][4) +  fakac * _cou[11][0][4)
          _cou[21][1][3) = _wmp[1]*_cou[11][1][4)  & 
               +  faka3 * (_cou[8][1][3) - fakaac * _cou[8][1][4) ) 
          _cou[27][1][3) = _wmp[2]*_cou[11][1][4) 
          _cou[24][2][3) = _wmp[0]*_cou[11][2][4) 
          _cou[21][2][3) = _wmp[1]*_cou[11][2][4)  & 
               +  faka3 * (_cou[8][2][3) - fakaac * _cou[8][2][4) ) + fakac * _cou[11][0][4)
          _cou[27][2][3) = _wmp[2]*_cou[11][2][4) 
          _cou[24][3][3) = _wmp[0]*_cou[11][3][4) 
          _cou[21][3][3) = _wmp[1]*_cou[11][3][4)  & 
               +  faka3 * (_cou[8][3][3) - fakaac * _cou[8][3][4) ) 
          _cou[27][3][3) = _wmp[2]*_cou[11][3][4) +  fakac * _cou[11][0][4)
          _cou[26][1][3) = _wmp[0]*_cou[12][1][4) +  fakac * _cou[12][0][4)
          _cou[28][1][3) = _wmp[1]*_cou[12][1][4) 
          _cou[22][1][3) = _wmp[2]*_cou[12][1][4)  & 
               +  faka3 * (_cou[9][1][3) - fakaac * _cou[9][1][4) ) 
          _cou[26][2][3) = _wmp[0]*_cou[12][2][4) 
          _cou[28][2][3) = _wmp[1]*_cou[12][2][4) +  fakac * _cou[12][0][4)
          _cou[22][2][3) = _wmp[2]*_cou[12][2][4)  & 
               +  faka3 * (_cou[9][2][3) - fakaac * _cou[9][2][4) ) 
          _cou[26][3][3) = _wmp[0]*_cou[12][3][4) 
          _cou[28][3][3) = _wmp[1]*_cou[12][3][4) 
          _cou[22][3][3) = _wmp[2]*_cou[12][3][4)  & 
               +  faka3 * (_cou[9][3][3) - fakaac * _cou[9][3][4) ) + fakac * _cou[12][0][4)
          _cou[31][1][3) = _wmp[1]*_cou[13][1][4)  & 
               +  faka * (_cou[7][1][3) - fakaac * _cou[7][1][4) ) 
          _cou[32][1][3) = _wmp[2]*_cou[13][1][4) 
          _cou[31][2][3) = _wmp[1]*_cou[13][2][4)  & 
               +  faka * (_cou[7][2][3) - fakaac * _cou[7][2][4) ) + fakac * _cou[13][0][4)
          _cou[32][2][3) = _wmp[2]*_cou[13][2][4) 
          _cou[31][3][3) = _wmp[1]*_cou[13][3][4)  & 
               +  faka * (_cou[7][3][3) - fakaac * _cou[7][3][4) ) 
          _cou[32][3][3) = _wmp[2]*_cou[13][3][4) +  fakac * _cou[13][0][4)
          _cou[33][1][3) = _wmp[2]*_cou[14][1][4) 
          _cou[33][2][3) = _wmp[2]*_cou[14][2][4) 
          _cou[33][3][3) = _wmp[2]*_cou[14][3][4) +  fakac * _cou[14][0][4)
          _cou[29][1][3) = _wmp[2]*_cou[15][1][4)  & 
               +  faka * (_cou[7][1][3) - fakaac * _cou[7][1][4) ) 
          _cou[29][2][3) = _wmp[2]*_cou[15][2][4)  & 
               +  faka * (_cou[7][2][3) - fakaac * _cou[7][2][4) ) 
          _cou[29][3][3) = _wmp[2]*_cou[15][3][4)  & 
               +  faka * (_cou[7][3][3) - fakaac * _cou[7][3][4) ) + fakac * _cou[15][0][4)
          _cou[34][1][3) = _wmp[1]*_cou[16][1][4) 
          _cou[34][2][3) = _wmp[1]*_cou[16][2][4) +  fakac * _cou[16][0][4)
          _cou[34][3][3) = _wmp[1]*_cou[16][3][4) 
          _cou[30][1][3) = _wmp[2]*_cou[17][1][4)  & 
               +  faka * (_cou[8][1][3) - fakaac * _cou[8][1][4) ) 
          _cou[30][2][3) = _wmp[2]*_cou[17][2][4)  & 
               +  faka * (_cou[8][2][3) - fakaac * _cou[8][2][4) ) 
          _cou[30][3][3) = _wmp[2]*_cou[17][3][4)  & 
               +  faka * (_cou[8][3][3) - fakaac * _cou[8][3][4) ) + fakac * _cou[17][0][4)


          _cou[7][20][2) = _wmp[0]*_cou[1][20][3)  & 
               +  faka * (_cou[0][20][2) - fakaac * _cou[0][20][3) ) + fakac4 * _cou[1][10][3)
          _cou[4][20][2) = _wmp[1]*_cou[1][20][3) 
          _cou[5][20][2) = _wmp[2]*_cou[1][20][3) 
          _cou[7][21][2) = _wmp[0]*_cou[1][21][3)  & 
               +  faka * (_cou[0][21][2) - fakaac * _cou[0][21][3) ) 
          _cou[4][21][2) = _wmp[1]*_cou[1][21][3) +  fakac4 * _cou[1][11][3)
          _cou[5][21][2) = _wmp[2]*_cou[1][21][3) 
          _cou[7][22][2) = _wmp[0]*_cou[1][22][3)  & 
               +  faka * (_cou[0][22][2) - fakaac * _cou[0][22][3) ) 
          _cou[4][22][2) = _wmp[1]*_cou[1][22][3) 
          _cou[5][22][2) = _wmp[2]*_cou[1][22][3) +  fakac4 * _cou[1][12][3)
          _cou[7][23][2) = _wmp[0]*_cou[1][23][3)  & 
               +  faka * (_cou[0][23][2) - fakaac * _cou[0][23][3) ) + fakac3 * _cou[1][13][3)
          _cou[4][23][2) = _wmp[1]*_cou[1][23][3) +  fakac * _cou[1][10][3)
          _cou[5][23][2) = _wmp[2]*_cou[1][23][3) 
          _cou[7][24][2) = _wmp[0]*_cou[1][24][3)  & 
               +  faka * (_cou[0][24][2) - fakaac * _cou[0][24][3) ) + fakac * _cou[1][11][3)
          _cou[4][24][2) = _wmp[1]*_cou[1][24][3) +  fakac3 * _cou[1][14][3)
          _cou[5][24][2) = _wmp[2]*_cou[1][24][3) 
          _cou[7][25][2) = _wmp[0]*_cou[1][25][3)  & 
               +  faka * (_cou[0][25][2) - fakaac * _cou[0][25][3) ) + fakac3 * _cou[1][15][3)
          _cou[4][25][2) = _wmp[1]*_cou[1][25][3) 
          _cou[5][25][2) = _wmp[2]*_cou[1][25][3) +  fakac * _cou[1][10][3)
          _cou[7][26][2) = _wmp[0]*_cou[1][26][3)  & 
               +  faka * (_cou[0][26][2) - fakaac * _cou[0][26][3) ) + fakac * _cou[1][12][3)
          _cou[4][26][2) = _wmp[1]*_cou[1][26][3) 
          _cou[5][26][2) = _wmp[2]*_cou[1][26][3) +  fakac3 * _cou[1][16][3)
          _cou[7][27][2) = _wmp[0]*_cou[1][27][3)  & 
               +  faka * (_cou[0][27][2) - fakaac * _cou[0][27][3) ) 
          _cou[4][27][2) = _wmp[1]*_cou[1][27][3) +  fakac3 * _cou[1][17][3)
          _cou[5][27][2) = _wmp[2]*_cou[1][27][3) +  fakac * _cou[1][11][3)
          _cou[7][28][2) = _wmp[0]*_cou[1][28][3)  & 
               +  faka * (_cou[0][28][2) - fakaac * _cou[0][28][3) ) 
          _cou[4][28][2) = _wmp[1]*_cou[1][28][3) +  fakac * _cou[1][12][3)
          _cou[5][28][2) = _wmp[2]*_cou[1][28][3) +  fakac3 * _cou[1][18][3)
          _cou[7][29][2) = _wmp[0]*_cou[1][29][3)  & 
               +  faka * (_cou[0][29][2) - fakaac * _cou[0][29][3) ) + fakac2 * _cou[1][16][3)
          _cou[4][29][2) = _wmp[1]*_cou[1][29][3) 
          _cou[5][29][2) = _wmp[2]*_cou[1][29][3) +  fakac2 * _cou[1][15][3)
          _cou[7][30][2) = _wmp[0]*_cou[1][30][3)  & 
               +  faka * (_cou[0][30][2) - fakaac * _cou[0][30][3) ) 
          _cou[4][30][2) = _wmp[1]*_cou[1][30][3) +  fakac2 * _cou[1][18][3)
          _cou[5][30][2) = _wmp[2]*_cou[1][30][3) +  fakac2 * _cou[1][17][3)
          _cou[7][31][2) = _wmp[0]*_cou[1][31][3)  & 
               +  faka * (_cou[0][31][2) - fakaac * _cou[0][31][3) ) + fakac2 * _cou[1][14][3)
          _cou[4][31][2) = _wmp[1]*_cou[1][31][3) +  fakac2 * _cou[1][13][3)
          _cou[5][31][2) = _wmp[2]*_cou[1][31][3) 
          _cou[7][32][2) = _wmp[0]*_cou[1][32][3)  & 
               +  faka * (_cou[0][32][2) - fakaac * _cou[0][32][3) ) + fakac2 * _cou[1][19][3)
          _cou[4][32][2) = _wmp[1]*_cou[1][32][3) +  fakac * _cou[1][15][3)
          _cou[5][32][2) = _wmp[2]*_cou[1][32][3) +  fakac * _cou[1][13][3)
          _cou[7][33][2) = _wmp[0]*_cou[1][33][3)  & 
               +  faka * (_cou[0][33][2) - fakaac * _cou[0][33][3) ) + fakac * _cou[1][17][3)
          _cou[4][33][2) = _wmp[1]*_cou[1][33][3) +  fakac2 * _cou[1][19][3)
          _cou[5][33][2) = _wmp[2]*_cou[1][33][3) +  fakac * _cou[1][14][3)
          _cou[7][34][2) = _wmp[0]*_cou[1][34][3)  & 
               +  faka * (_cou[0][34][2) - fakaac * _cou[0][34][3) ) + fakac * _cou[1][18][3)
          _cou[4][34][2) = _wmp[1]*_cou[1][34][3) +  fakac * _cou[1][16][3)
          _cou[5][34][2) = _wmp[2]*_cou[1][34][3) +  fakac2 * _cou[1][19][3)
          _cou[8][20][2) = _wmp[1]*_cou[2][20][3)  & 
               +  faka * (_cou[0][20][2) - fakaac * _cou[0][20][3) ) 
          _cou[6][20][2) = _wmp[2]*_cou[2][20][3) 
          _cou[8][21][2) = _wmp[1]*_cou[2][21][3)  & 
               +  faka * (_cou[0][21][2) - fakaac * _cou[0][21][3) ) + fakac4 * _cou[2][11][3)
          _cou[6][21][2) = _wmp[2]*_cou[2][21][3) 
          _cou[8][22][2) = _wmp[1]*_cou[2][22][3)  & 
               +  faka * (_cou[0][22][2) - fakaac * _cou[0][22][3) ) 
          _cou[6][22][2) = _wmp[2]*_cou[2][22][3) +  fakac4 * _cou[2][12][3)
          _cou[8][23][2) = _wmp[1]*_cou[2][23][3)  & 
               +  faka * (_cou[0][23][2) - fakaac * _cou[0][23][3) ) + fakac * _cou[2][10][3)
          _cou[6][23][2) = _wmp[2]*_cou[2][23][3) 
          _cou[8][24][2) = _wmp[1]*_cou[2][24][3)  & 
               +  faka * (_cou[0][24][2) - fakaac * _cou[0][24][3) ) + fakac3 * _cou[2][14][3)
          _cou[6][24][2) = _wmp[2]*_cou[2][24][3) 
          _cou[8][25][2) = _wmp[1]*_cou[2][25][3)  & 
               +  faka * (_cou[0][25][2) - fakaac * _cou[0][25][3) ) 
          _cou[6][25][2) = _wmp[2]*_cou[2][25][3) +  fakac * _cou[2][10][3)
          _cou[8][26][2) = _wmp[1]*_cou[2][26][3)  & 
               +  faka * (_cou[0][26][2) - fakaac * _cou[0][26][3) ) 
          _cou[6][26][2) = _wmp[2]*_cou[2][26][3) +  fakac3 * _cou[2][16][3)
          _cou[8][27][2) = _wmp[1]*_cou[2][27][3)  & 
               +  faka * (_cou[0][27][2) - fakaac * _cou[0][27][3) ) + fakac3 * _cou[2][17][3)
          _cou[6][27][2) = _wmp[2]*_cou[2][27][3) +  fakac * _cou[2][11][3)
          _cou[8][28][2) = _wmp[1]*_cou[2][28][3)  & 
               +  faka * (_cou[0][28][2) - fakaac * _cou[0][28][3) ) + fakac * _cou[2][12][3)
          _cou[6][28][2) = _wmp[2]*_cou[2][28][3) +  fakac3 * _cou[2][18][3)
          _cou[8][29][2) = _wmp[1]*_cou[2][29][3)  & 
               +  faka * (_cou[0][29][2) - fakaac * _cou[0][29][3) ) 
          _cou[6][29][2) = _wmp[2]*_cou[2][29][3) +  fakac2 * _cou[2][15][3)
          _cou[8][30][2) = _wmp[1]*_cou[2][30][3)  & 
               +  faka * (_cou[0][30][2) - fakaac * _cou[0][30][3) ) + fakac2 * _cou[2][18][3)
          _cou[6][30][2) = _wmp[2]*_cou[2][30][3) +  fakac2 * _cou[2][17][3)
          _cou[8][31][2) = _wmp[1]*_cou[2][31][3)  & 
               +  faka * (_cou[0][31][2) - fakaac * _cou[0][31][3) ) + fakac2 * _cou[2][13][3)
          _cou[6][31][2) = _wmp[2]*_cou[2][31][3) 
          _cou[8][32][2) = _wmp[1]*_cou[2][32][3)  & 
               +  faka * (_cou[0][32][2) - fakaac * _cou[0][32][3) ) + fakac * _cou[2][15][3)
          _cou[6][32][2) = _wmp[2]*_cou[2][32][3) +  fakac * _cou[2][13][3)
          _cou[8][33][2) = _wmp[1]*_cou[2][33][3)  & 
               +  faka * (_cou[0][33][2) - fakaac * _cou[0][33][3) ) + fakac2 * _cou[2][19][3)
          _cou[6][33][2) = _wmp[2]*_cou[2][33][3) +  fakac * _cou[2][14][3)
          _cou[8][34][2) = _wmp[1]*_cou[2][34][3)  & 
               +  faka * (_cou[0][34][2) - fakaac * _cou[0][34][3) ) + fakac * _cou[2][16][3)
          _cou[6][34][2) = _wmp[2]*_cou[2][34][3) +  fakac2 * _cou[2][19][3)
          _cou[9][20][2) = _wmp[2]*_cou[3][20][3)  & 
               +  faka * (_cou[0][20][2) - fakaac * _cou[0][20][3) ) 
          _cou[9][21][2) = _wmp[2]*_cou[3][21][3)  & 
               +  faka * (_cou[0][21][2) - fakaac * _cou[0][21][3) ) 
          _cou[9][22][2) = _wmp[2]*_cou[3][22][3)  & 
               +  faka * (_cou[0][22][2) - fakaac * _cou[0][22][3) ) + fakac4 * _cou[3][12][3)
          _cou[9][23][2) = _wmp[2]*_cou[3][23][3)  & 
               +  faka * (_cou[0][23][2) - fakaac * _cou[0][23][3) ) 
          _cou[9][24][2) = _wmp[2]*_cou[3][24][3)  & 
               +  faka * (_cou[0][24][2) - fakaac * _cou[0][24][3) ) 
          _cou[9][25][2) = _wmp[2]*_cou[3][25][3)  & 
               +  faka * (_cou[0][25][2) - fakaac * _cou[0][25][3) ) + fakac * _cou[3][10][3)
          _cou[9][26][2) = _wmp[2]*_cou[3][26][3)  & 
               +  faka * (_cou[0][26][2) - fakaac * _cou[0][26][3) ) + fakac3 * _cou[3][16][3)
          _cou[9][27][2) = _wmp[2]*_cou[3][27][3)  & 
               +  faka * (_cou[0][27][2) - fakaac * _cou[0][27][3) ) + fakac * _cou[3][11][3)
          _cou[9][28][2) = _wmp[2]*_cou[3][28][3)  & 
               +  faka * (_cou[0][28][2) - fakaac * _cou[0][28][3) ) + fakac3 * _cou[3][18][3)
          _cou[9][29][2) = _wmp[2]*_cou[3][29][3)  & 
               +  faka * (_cou[0][29][2) - fakaac * _cou[0][29][3) ) + fakac2 * _cou[3][15][3)
          _cou[9][30][2) = _wmp[2]*_cou[3][30][3)  & 
               +  faka * (_cou[0][30][2) - fakaac * _cou[0][30][3) ) + fakac2 * _cou[3][17][3)
          _cou[9][31][2) = _wmp[2]*_cou[3][31][3)  & 
               +  faka * (_cou[0][31][2) - fakaac * _cou[0][31][3) ) 
          _cou[9][32][2) = _wmp[2]*_cou[3][32][3)  & 
               +  faka * (_cou[0][32][2) - fakaac * _cou[0][32][3) ) + fakac * _cou[3][13][3)
          _cou[9][33][2) = _wmp[2]*_cou[3][33][3)  & 
               +  faka * (_cou[0][33][2) - fakaac * _cou[0][33][3) ) + fakac * _cou[3][14][3)
          _cou[9][34][2) = _wmp[2]*_cou[3][34][3)  & 
               +  faka * (_cou[0][34][2) - fakaac * _cou[0][34][3) ) + fakac2 * _cou[3][19][3)


          _cou[13][10][2) = _wmp[0]*_cou[4][10][3)  & 
               +  faka * (_cou[2][10][2) - fakaac * _cou[2][10][3) ) + fakac3 * _cou[4][7][3)
          _cou[14][10][2) = _wmp[1]*_cou[4][10][3)  & 
               +  faka * (_cou[1][10][2) - fakaac * _cou[1][10][3) ) 
          _cou[19][10][2) = _wmp[2]*_cou[4][10][3) 
          _cou[13][11][2) = _wmp[0]*_cou[4][11][3)  & 
               +  faka * (_cou[2][11][2) - fakaac * _cou[2][11][3) ) 
          _cou[14][11][2) = _wmp[1]*_cou[4][11][3)  & 
               +  faka * (_cou[1][11][2) - fakaac * _cou[1][11][3) ) + fakac3 * _cou[4][8][3)
          _cou[19][11][2) = _wmp[2]*_cou[4][11][3) 
          _cou[13][12][2) = _wmp[0]*_cou[4][12][3)  & 
               +  faka * (_cou[2][12][2) - fakaac * _cou[2][12][3) ) 
          _cou[14][12][2) = _wmp[1]*_cou[4][12][3)  & 
               +  faka * (_cou[1][12][2) - fakaac * _cou[1][12][3) ) 
          _cou[19][12][2) = _wmp[2]*_cou[4][12][3) +  fakac3 * _cou[4][9][3)
          _cou[13][13][2) = _wmp[0]*_cou[4][13][3)  & 
               +  faka * (_cou[2][13][2) - fakaac * _cou[2][13][3) ) + fakac2 * _cou[4][4][3)
          _cou[14][13][2) = _wmp[1]*_cou[4][13][3)  & 
               +  faka * (_cou[1][13][2) - fakaac * _cou[1][13][3) ) + fakac * _cou[4][7][3)
          _cou[19][13][2) = _wmp[2]*_cou[4][13][3) 
          _cou[13][14][2) = _wmp[0]*_cou[4][14][3)  & 
               +  faka * (_cou[2][14][2) - fakaac * _cou[2][14][3) ) + fakac * _cou[4][8][3)
          _cou[14][14][2) = _wmp[1]*_cou[4][14][3)  & 
               +  faka * (_cou[1][14][2) - fakaac * _cou[1][14][3) ) + fakac2 * _cou[4][4][3)
          _cou[19][14][2) = _wmp[2]*_cou[4][14][3) 
          _cou[13][15][2) = _wmp[0]*_cou[4][15][3)  & 
               +  faka * (_cou[2][15][2) - fakaac * _cou[2][15][3) ) + fakac2 * _cou[4][5][3)
          _cou[14][15][2) = _wmp[1]*_cou[4][15][3)  & 
               +  faka * (_cou[1][15][2) - fakaac * _cou[1][15][3) ) 
          _cou[19][15][2) = _wmp[2]*_cou[4][15][3) +  fakac * _cou[4][7][3)
          _cou[13][16][2) = _wmp[0]*_cou[4][16][3)  & 
               +  faka * (_cou[2][16][2) - fakaac * _cou[2][16][3) ) + fakac * _cou[4][9][3)
          _cou[14][16][2) = _wmp[1]*_cou[4][16][3)  & 
               +  faka * (_cou[1][16][2) - fakaac * _cou[1][16][3) ) 
          _cou[19][16][2) = _wmp[2]*_cou[4][16][3) +  fakac2 * _cou[4][5][3)
          _cou[13][17][2) = _wmp[0]*_cou[4][17][3)  & 
               +  faka * (_cou[2][17][2) - fakaac * _cou[2][17][3) ) 
          _cou[14][17][2) = _wmp[1]*_cou[4][17][3)  & 
               +  faka * (_cou[1][17][2) - fakaac * _cou[1][17][3) ) + fakac2 * _cou[4][6][3)
          _cou[19][17][2) = _wmp[2]*_cou[4][17][3) +  fakac * _cou[4][8][3)
          _cou[13][18][2) = _wmp[0]*_cou[4][18][3)  & 
               +  faka * (_cou[2][18][2) - fakaac * _cou[2][18][3) ) 
          _cou[14][18][2) = _wmp[1]*_cou[4][18][3)  & 
               +  faka * (_cou[1][18][2) - fakaac * _cou[1][18][3) ) + fakac * _cou[4][9][3)
          _cou[19][18][2) = _wmp[2]*_cou[4][18][3) +  fakac2 * _cou[4][6][3)
          _cou[13][19][2) = _wmp[0]*_cou[4][19][3)  & 
               +  faka * (_cou[2][19][2) - fakaac * _cou[2][19][3) ) + fakac * _cou[4][6][3)
          _cou[14][19][2) = _wmp[1]*_cou[4][19][3)  & 
               +  faka * (_cou[1][19][2) - fakaac * _cou[1][19][3) ) + fakac * _cou[4][5][3)
          _cou[19][19][2) = _wmp[2]*_cou[4][19][3) +  fakac * _cou[4][4][3)
          _cou[15][10][2) = _wmp[0]*_cou[5][10][3)  & 
               +  faka * (_cou[3][10][2) - fakaac * _cou[3][10][3) ) + fakac3 * _cou[5][7][3)
          _cou[16][10][2) = _wmp[2]*_cou[5][10][3)  & 
               +  faka * (_cou[1][10][2) - fakaac * _cou[1][10][3) ) 
          _cou[15][11][2) = _wmp[0]*_cou[5][11][3)  & 
               +  faka * (_cou[3][11][2) - fakaac * _cou[3][11][3) ) 
          _cou[16][11][2) = _wmp[2]*_cou[5][11][3)  & 
               +  faka * (_cou[1][11][2) - fakaac * _cou[1][11][3) ) 
          _cou[15][12][2) = _wmp[0]*_cou[5][12][3)  & 
               +  faka * (_cou[3][12][2) - fakaac * _cou[3][12][3) ) 
          _cou[16][12][2) = _wmp[2]*_cou[5][12][3)  & 
               +  faka * (_cou[1][12][2) - fakaac * _cou[1][12][3) ) + fakac3 * _cou[5][9][3)
          _cou[15][13][2) = _wmp[0]*_cou[5][13][3)  & 
               +  faka * (_cou[3][13][2) - fakaac * _cou[3][13][3) ) + fakac2 * _cou[5][4][3)
          _cou[16][13][2) = _wmp[2]*_cou[5][13][3)  & 
               +  faka * (_cou[1][13][2) - fakaac * _cou[1][13][3) ) 
          _cou[15][14][2) = _wmp[0]*_cou[5][14][3)  & 
               +  faka * (_cou[3][14][2) - fakaac * _cou[3][14][3) ) + fakac * _cou[5][8][3)
          _cou[16][14][2) = _wmp[2]*_cou[5][14][3)  & 
               +  faka * (_cou[1][14][2) - fakaac * _cou[1][14][3) ) 
          _cou[15][15][2) = _wmp[0]*_cou[5][15][3)  & 
               +  faka * (_cou[3][15][2) - fakaac * _cou[3][15][3) ) + fakac2 * _cou[5][5][3)
          _cou[16][15][2) = _wmp[2]*_cou[5][15][3)  & 
               +  faka * (_cou[1][15][2) - fakaac * _cou[1][15][3) ) + fakac * _cou[5][7][3)
          _cou[15][16][2) = _wmp[0]*_cou[5][16][3)  & 
               +  faka * (_cou[3][16][2) - fakaac * _cou[3][16][3) ) + fakac * _cou[5][9][3)
          _cou[16][16][2) = _wmp[2]*_cou[5][16][3)  & 
               +  faka * (_cou[1][16][2) - fakaac * _cou[1][16][3) ) + fakac2 * _cou[5][5][3)
          _cou[15][17][2) = _wmp[0]*_cou[5][17][3)  & 
               +  faka * (_cou[3][17][2) - fakaac * _cou[3][17][3) ) 
          _cou[16][17][2) = _wmp[2]*_cou[5][17][3)  & 
               +  faka * (_cou[1][17][2) - fakaac * _cou[1][17][3) ) + fakac * _cou[5][8][3)
          _cou[15][18][2) = _wmp[0]*_cou[5][18][3)  & 
               +  faka * (_cou[3][18][2) - fakaac * _cou[3][18][3) ) 
          _cou[16][18][2) = _wmp[2]*_cou[5][18][3)  & 
               +  faka * (_cou[1][18][2) - fakaac * _cou[1][18][3) ) + fakac2 * _cou[5][6][3)
          _cou[15][19][2) = _wmp[0]*_cou[5][19][3)  & 
               +  faka * (_cou[3][19][2) - fakaac * _cou[3][19][3) ) + fakac * _cou[5][6][3)
          _cou[16][19][2) = _wmp[2]*_cou[5][19][3)  & 
               +  faka * (_cou[1][19][2) - fakaac * _cou[1][19][3) ) + fakac * _cou[5][4][3)
          _cou[17][10][2) = _wmp[1]*_cou[6][10][3)  & 
               +  faka * (_cou[3][10][2) - fakaac * _cou[3][10][3) ) 
          _cou[18][10][2) = _wmp[2]*_cou[6][10][3)  & 
               +  faka * (_cou[2][10][2) - fakaac * _cou[2][10][3) ) 
          _cou[17][11][2) = _wmp[1]*_cou[6][11][3)  & 
               +  faka * (_cou[3][11][2) - fakaac * _cou[3][11][3) ) + fakac3 * _cou[6][8][3)
          _cou[18][11][2) = _wmp[2]*_cou[6][11][3)  & 
               +  faka * (_cou[2][11][2) - fakaac * _cou[2][11][3) ) 
          _cou[17][12][2) = _wmp[1]*_cou[6][12][3)  & 
               +  faka * (_cou[3][12][2) - fakaac * _cou[3][12][3) ) 
          _cou[18][12][2) = _wmp[2]*_cou[6][12][3)  & 
               +  faka * (_cou[2][12][2) - fakaac * _cou[2][12][3) ) + fakac3 * _cou[6][9][3)
          _cou[17][13][2) = _wmp[1]*_cou[6][13][3)  & 
               +  faka * (_cou[3][13][2) - fakaac * _cou[3][13][3) ) + fakac * _cou[6][7][3)
          _cou[18][13][2) = _wmp[2]*_cou[6][13][3)  & 
               +  faka * (_cou[2][13][2) - fakaac * _cou[2][13][3) ) 
          _cou[17][14][2) = _wmp[1]*_cou[6][14][3)  & 
               +  faka * (_cou[3][14][2) - fakaac * _cou[3][14][3) ) + fakac2 * _cou[6][4][3)
          _cou[18][14][2) = _wmp[2]*_cou[6][14][3)  & 
               +  faka * (_cou[2][14][2) - fakaac * _cou[2][14][3) ) 
          _cou[17][15][2) = _wmp[1]*_cou[6][15][3)  & 
               +  faka * (_cou[3][15][2) - fakaac * _cou[3][15][3) ) 
          _cou[18][15][2) = _wmp[2]*_cou[6][15][3)  & 
               +  faka * (_cou[2][15][2) - fakaac * _cou[2][15][3) ) + fakac * _cou[6][7][3)
          _cou[17][16][2) = _wmp[1]*_cou[6][16][3)  & 
               +  faka * (_cou[3][16][2) - fakaac * _cou[3][16][3) ) 
          _cou[18][16][2) = _wmp[2]*_cou[6][16][3)  & 
               +  faka * (_cou[2][16][2) - fakaac * _cou[2][16][3) ) + fakac2 * _cou[6][5][3)
          _cou[17][17][2) = _wmp[1]*_cou[6][17][3)  & 
               +  faka * (_cou[3][17][2) - fakaac * _cou[3][17][3) ) + fakac2 * _cou[6][6][3)
          _cou[18][17][2) = _wmp[2]*_cou[6][17][3)  & 
               +  faka * (_cou[2][17][2) - fakaac * _cou[2][17][3) ) + fakac * _cou[6][8][3)
          _cou[17][18][2) = _wmp[1]*_cou[6][18][3)  & 
               +  faka * (_cou[3][18][2) - fakaac * _cou[3][18][3) ) + fakac * _cou[6][9][3)
          _cou[18][18][2) = _wmp[2]*_cou[6][18][3)  & 
               +  faka * (_cou[2][18][2) - fakaac * _cou[2][18][3) ) + fakac2 * _cou[6][6][3)
          _cou[17][19][2) = _wmp[1]*_cou[6][19][3)  & 
               +  faka * (_cou[3][19][2) - fakaac * _cou[3][19][3) ) + fakac * _cou[6][5][3)
          _cou[18][19][2) = _wmp[2]*_cou[6][19][3)  & 
               +  faka * (_cou[2][19][2) - fakaac * _cou[2][19][3) ) + fakac * _cou[6][4][3)
          _cou[10][10][2) = _wmp[0]*_cou[7][10][3)  & 
               + faka2 * (_cou[1][10][2) - fakaac * _cou[1][10][3) ) + fakac3 * _cou[7][7][3)
          _cou[10][11][2) = _wmp[0]*_cou[7][11][3)  & 
               +  faka2 * (_cou[1][11][2) - fakaac * _cou[1][11][3) ) 
          _cou[10][12][2) = _wmp[0]*_cou[7][12][3)  & 
               +  faka2 * (_cou[1][12][2) - fakaac * _cou[1][12][3) ) 
          _cou[10][13][2) = _wmp[0]*_cou[7][13][3)  & 
               + faka2 * (_cou[1][13][2) - fakaac * _cou[1][13][3) ) + fakac2 * _cou[7][4][3)
          _cou[10][14][2) = _wmp[0]*_cou[7][14][3)  & 
               +  faka2 * (_cou[1][14][2) - fakaac * _cou[1][14][3) ) + fakac * _cou[7][8][3)
          _cou[10][15][2) = _wmp[0]*_cou[7][15][3)  & 
               + faka2 * (_cou[1][15][2) - fakaac * _cou[1][15][3) ) + fakac2 * _cou[7][5][3)
          _cou[10][16][2) = _wmp[0]*_cou[7][16][3)  & 
               +  faka2 * (_cou[1][16][2) - fakaac * _cou[1][16][3) ) + fakac * _cou[7][9][3)
          _cou[10][17][2) = _wmp[0]*_cou[7][17][3)  & 
               +  faka2 * (_cou[1][17][2) - fakaac * _cou[1][17][3) ) 
          _cou[10][18][2) = _wmp[0]*_cou[7][18][3)  & 
               +  faka2 * (_cou[1][18][2) - fakaac * _cou[1][18][3) ) 
          _cou[10][19][2) = _wmp[0]*_cou[7][19][3)  & 
               +  faka2 * (_cou[1][19][2) - fakaac * _cou[1][19][3) ) + fakac * _cou[7][6][3)
          _cou[11][10][2) = _wmp[1]*_cou[8][10][3)  & 
               +  faka2 * (_cou[2][10][2) - fakaac * _cou[2][10][3) ) 
          _cou[11][11][2) = _wmp[1]*_cou[8][11][3)  & 
               + faka2 * (_cou[2][11][2) - fakaac * _cou[2][11][3) ) + fakac3 * _cou[8][8][3)
          _cou[11][12][2) = _wmp[1]*_cou[8][12][3)  & 
               +  faka2 * (_cou[2][12][2) - fakaac * _cou[2][12][3) ) 
          _cou[11][13][2) = _wmp[1]*_cou[8][13][3)  & 
               +  faka2 * (_cou[2][13][2) - fakaac * _cou[2][13][3) ) + fakac * _cou[8][7][3)
          _cou[11][14][2) = _wmp[1]*_cou[8][14][3)  & 
               + faka2 * (_cou[2][14][2) - fakaac * _cou[2][14][3) ) + fakac2 * _cou[8][4][3)
          _cou[11][15][2) = _wmp[1]*_cou[8][15][3)  & 
               +  faka2 * (_cou[2][15][2) - fakaac * _cou[2][15][3) ) 
          _cou[11][16][2) = _wmp[1]*_cou[8][16][3)  & 
               +  faka2 * (_cou[2][16][2) - fakaac * _cou[2][16][3) ) 
          _cou[11][17][2) = _wmp[1]*_cou[8][17][3)  & 
               + faka2 * (_cou[2][17][2) - fakaac * _cou[2][17][3) ) + fakac2 * _cou[8][6][3)
          _cou[11][18][2) = _wmp[1]*_cou[8][18][3)  & 
               +  faka2 * (_cou[2][18][2) - fakaac * _cou[2][18][3) ) + fakac * _cou[8][9][3)
          _cou[11][19][2) = _wmp[1]*_cou[8][19][3)  & 
               +  faka2 * (_cou[2][19][2) - fakaac * _cou[2][19][3) ) + fakac * _cou[8][5][3)
          _cou[12][10][2) = _wmp[2]*_cou[9][10][3)  & 
               +  faka2 * (_cou[3][10][2) - fakaac * _cou[3][10][3) ) 
          _cou[12][11][2) = _wmp[2]*_cou[9][11][3)  & 
               +  faka2 * (_cou[3][11][2) - fakaac * _cou[3][11][3) ) 
          _cou[12][12][2) = _wmp[2]*_cou[9][12][3)  & 
               + faka2 * (_cou[3][12][2) - fakaac * _cou[3][12][3) ) + fakac3 * _cou[9][9][3)
          _cou[12][13][2) = _wmp[2]*_cou[9][13][3)  & 
               +  faka2 * (_cou[3][13][2) - fakaac * _cou[3][13][3) ) 
          _cou[12][14][2) = _wmp[2]*_cou[9][14][3)  & 
               +  faka2 * (_cou[3][14][2) - fakaac * _cou[3][14][3) ) 
          _cou[12][15][2) = _wmp[2]*_cou[9][15][3)  & 
               +  faka2 * (_cou[3][15][2) - fakaac * _cou[3][15][3) ) + fakac * _cou[9][7][3)
          _cou[12][16][2) = _wmp[2]*_cou[9][16][3)  & 
               + faka2 * (_cou[3][16][2) - fakaac * _cou[3][16][3) ) + fakac2 * _cou[9][5][3)
          _cou[12][17][2) = _wmp[2]*_cou[9][17][3)  & 
               +  faka2 * (_cou[3][17][2) - fakaac * _cou[3][17][3) ) + fakac * _cou[9][8][3)
          _cou[12][18][2) = _wmp[2]*_cou[9][18][3)  & 
               + faka2 * (_cou[3][18][2) - fakaac * _cou[3][18][3) ) + fakac2 * _cou[9][6][3)
          _cou[12][19][2) = _wmp[2]*_cou[9][19][3)  & 
               +  faka2 * (_cou[3][19][2) - fakaac * _cou[3][19][3) ) + fakac * _cou[9][4][3)

          _cou[20][4][2) = _wmp[0]*_cou[10][4][3)  & 
               +  faka3 * (_cou[7][4][2) - fakaac * _cou[7][4][3) ) + fakac * _cou[10][2][3)
          _cou[23][4][2) = _wmp[1]*_cou[10][4][3) +  fakac * _cou[10][1][3)
          _cou[25][4][2) = _wmp[2]*_cou[10][4][3) 
          _cou[20][5][2) = _wmp[0]*_cou[10][5][3)  & 
               +  faka3 * (_cou[7][5][2) - fakaac * _cou[7][5][3) ) + fakac * _cou[10][3][3)
          _cou[23][5][2) = _wmp[1]*_cou[10][5][3) 
          _cou[25][5][2) = _wmp[2]*_cou[10][5][3) +  fakac * _cou[10][1][3)
          _cou[20][6][2) = _wmp[0]*_cou[10][6][3)  & 
               +  faka3 * (_cou[7][6][2) - fakaac * _cou[7][6][3) ) 
          _cou[23][6][2) = _wmp[1]*_cou[10][6][3) +  fakac * _cou[10][3][3)
          _cou[25][6][2) = _wmp[2]*_cou[10][6][3) +  fakac * _cou[10][2][3)
          _cou[20][7][2) = _wmp[0]*_cou[10][7][3)  & 
               + faka3 * (_cou[7][7][2) - fakaac * _cou[7][7][3) ) + fakac2 * _cou[10][1][3)
          _cou[23][7][2) = _wmp[1]*_cou[10][7][3) 
          _cou[25][7][2) = _wmp[2]*_cou[10][7][3) 
          _cou[20][8][2) = _wmp[0]*_cou[10][8][3)  & 
               +  faka3 * (_cou[7][8][2) - fakaac * _cou[7][8][3) ) 
          _cou[23][8][2) = _wmp[1]*_cou[10][8][3) +  fakac2 * _cou[10][2][3)
          _cou[25][8][2) = _wmp[2]*_cou[10][8][3) 
          _cou[20][9][2) = _wmp[0]*_cou[10][9][3)  & 
               +  faka3 * (_cou[7][9][2) - fakaac * _cou[7][9][3) ) 
          _cou[23][9][2) = _wmp[1]*_cou[10][9][3) 
          _cou[25][9][2) = _wmp[2]*_cou[10][9][3) +  fakac2 * _cou[10][3][3)
          _cou[24][4][2) = _wmp[0]*_cou[11][4][3) +  fakac * _cou[11][2][3)
          _cou[21][4][2) = _wmp[1]*_cou[11][4][3)  & 
               +  faka3 * (_cou[8][4][2) - fakaac * _cou[8][4][3) ) + fakac * _cou[11][1][3)
          _cou[27][4][2) = _wmp[2]*_cou[11][4][3) 
          _cou[24][5][2) = _wmp[0]*_cou[11][5][3) +  fakac * _cou[11][3][3)
          _cou[21][5][2) = _wmp[1]*_cou[11][5][3)  & 
               +  faka3 * (_cou[8][5][2) - fakaac * _cou[8][5][3) ) 
          _cou[27][5][2) = _wmp[2]*_cou[11][5][3) +  fakac * _cou[11][1][3)
          _cou[24][6][2) = _wmp[0]*_cou[11][6][3) 
          _cou[21][6][2) = _wmp[1]*_cou[11][6][3)  & 
               +  faka3 * (_cou[8][6][2) - fakaac * _cou[8][6][3) ) + fakac * _cou[11][3][3)
          _cou[27][6][2) = _wmp[2]*_cou[11][6][3) +  fakac * _cou[11][2][3)
          _cou[24][7][2) = _wmp[0]*_cou[11][7][3) +  fakac2 * _cou[11][1][3)
          _cou[21][7][2) = _wmp[1]*_cou[11][7][3)  & 
               +  faka3 * (_cou[8][7][2) - fakaac * _cou[8][7][3) ) 
          _cou[27][7][2) = _wmp[2]*_cou[11][7][3) 
          _cou[24][8][2) = _wmp[0]*_cou[11][8][3) 
          _cou[21][8][2) = _wmp[1]*_cou[11][8][3)  & 
               + faka3 * (_cou[8][8][2) - fakaac * _cou[8][8][3) ) + fakac2 * _cou[11][2][3)
          _cou[27][8][2) = _wmp[2]*_cou[11][8][3) 
          _cou[24][9][2) = _wmp[0]*_cou[11][9][3) 
          _cou[21][9][2) = _wmp[1]*_cou[11][9][3)  & 
               +  faka3 * (_cou[8][9][2) - fakaac * _cou[8][9][3) ) 
          _cou[27][9][2) = _wmp[2]*_cou[11][9][3) +  fakac2 * _cou[11][3][3)
          _cou[26][4][2) = _wmp[0]*_cou[12][4][3) +  fakac * _cou[12][2][3)
          _cou[28][4][2) = _wmp[1]*_cou[12][4][3) +  fakac * _cou[12][1][3)
          _cou[22][4][2) = _wmp[2]*_cou[12][4][3)  & 
               +  faka3 * (_cou[9][4][2) - fakaac * _cou[9][4][3) ) 
          _cou[26][5][2) = _wmp[0]*_cou[12][5][3) +  fakac * _cou[12][3][3)
          _cou[28][5][2) = _wmp[1]*_cou[12][5][3) 
          _cou[22][5][2) = _wmp[2]*_cou[12][5][3)  & 
               +  faka3 * (_cou[9][5][2) - fakaac * _cou[9][5][3) ) + fakac * _cou[12][1][3)
          _cou[26][6][2) = _wmp[0]*_cou[12][6][3) 
          _cou[28][6][2) = _wmp[1]*_cou[12][6][3) +  fakac * _cou[12][3][3)
          _cou[22][6][2) = _wmp[2]*_cou[12][6][3)  & 
               +  faka3 * (_cou[9][6][2) - fakaac * _cou[9][6][3) ) + fakac * _cou[12][2][3)
          _cou[26][7][2) = _wmp[0]*_cou[12][7][3) +  fakac2 * _cou[12][1][3)
          _cou[28][7][2) = _wmp[1]*_cou[12][7][3) 
          _cou[22][7][2) = _wmp[2]*_cou[12][7][3)  & 
               +  faka3 * (_cou[9][7][2) - fakaac * _cou[9][7][3) ) 
          _cou[26][8][2) = _wmp[0]*_cou[12][8][3) 
          _cou[28][8][2) = _wmp[1]*_cou[12][8][3) +  fakac2 * _cou[12][2][3)
          _cou[22][8][2) = _wmp[2]*_cou[12][8][3)  & 
               +  faka3 * (_cou[9][8][2) - fakaac * _cou[9][8][3) ) 
          _cou[26][9][2) = _wmp[0]*_cou[12][9][3) 
          _cou[28][9][2) = _wmp[1]*_cou[12][9][3) 
          _cou[22][9][2) = _wmp[2]*_cou[12][9][3)  & 
               + faka3 * (_cou[9][9][2) - fakaac * _cou[9][9][3) ) + fakac2 * _cou[12][3][3)
          _cou[31][4][2) = _wmp[1]*_cou[13][4][3)  & 
               +  faka * (_cou[7][4][2) - fakaac * _cou[7][4][3) ) + fakac * _cou[13][1][3)
          _cou[32][4][2) = _wmp[2]*_cou[13][4][3) 
          _cou[31][5][2) = _wmp[1]*_cou[13][5][3)  & 
               +  faka * (_cou[7][5][2) - fakaac * _cou[7][5][3) ) 
          _cou[32][5][2) = _wmp[2]*_cou[13][5][3) +  fakac * _cou[13][1][3)
          _cou[31][6][2) = _wmp[1]*_cou[13][6][3)  & 
               +  faka * (_cou[7][6][2) - fakaac * _cou[7][6][3) ) + fakac * _cou[13][3][3)
          _cou[32][6][2) = _wmp[2]*_cou[13][6][3) +  fakac * _cou[13][2][3)
          _cou[31][7][2) = _wmp[1]*_cou[13][7][3)  & 
               +  faka * (_cou[7][7][2) - fakaac * _cou[7][7][3) ) 
          _cou[32][7][2) = _wmp[2]*_cou[13][7][3) 
          _cou[31][8][2) = _wmp[1]*_cou[13][8][3)  & 
               +  faka * (_cou[7][8][2) - fakaac * _cou[7][8][3) ) + fakac2 * _cou[13][2][3)
          _cou[32][8][2) = _wmp[2]*_cou[13][8][3) 
          _cou[31][9][2) = _wmp[1]*_cou[13][9][3)  & 
               +  faka * (_cou[7][9][2) - fakaac * _cou[7][9][3) ) 
          _cou[32][9][2) = _wmp[2]*_cou[13][9][3) +  fakac2 * _cou[13][3][3)
          _cou[33][4][2) = _wmp[2]*_cou[14][4][3) 
          _cou[33][5][2) = _wmp[2]*_cou[14][5][3) +  fakac * _cou[14][1][3)
          _cou[33][6][2) = _wmp[2]*_cou[14][6][3) +  fakac * _cou[14][2][3)
          _cou[33][7][2) = _wmp[2]*_cou[14][7][3) 
          _cou[33][8][2) = _wmp[2]*_cou[14][8][3) 
          _cou[33][9][2) = _wmp[2]*_cou[14][9][3) +  fakac2 * _cou[14][3][3)
          _cou[29][4][2) = _wmp[2]*_cou[15][4][3)  & 
               +  faka * (_cou[7][4][2) - fakaac * _cou[7][4][3) ) 
          _cou[29][5][2) = _wmp[2]*_cou[15][5][3)  & 
               +  faka * (_cou[7][5][2) - fakaac * _cou[7][5][3) ) + fakac * _cou[15][1][3)
          _cou[29][6][2) = _wmp[2]*_cou[15][6][3)  & 
               +  faka * (_cou[7][6][2) - fakaac * _cou[7][6][3) ) + fakac * _cou[15][2][3)
          _cou[29][7][2) = _wmp[2]*_cou[15][7][3)  & 
               +  faka * (_cou[7][7][2) - fakaac * _cou[7][7][3) ) 
          _cou[29][8][2) = _wmp[2]*_cou[15][8][3)  & 
               +  faka * (_cou[7][8][2) - fakaac * _cou[7][8][3) ) 
          _cou[29][9][2) = _wmp[2]*_cou[15][9][3)  & 
               +  faka * (_cou[7][9][2) - fakaac * _cou[7][9][3) ) + fakac2 * _cou[15][3][3)
          _cou[34][4][2) = _wmp[1]*_cou[16][4][3) +  fakac * _cou[16][1][3)
          _cou[34][5][2) = _wmp[1]*_cou[16][5][3) 
          _cou[34][6][2) = _wmp[1]*_cou[16][6][3) +  fakac * _cou[16][3][3)
          _cou[34][7][2) = _wmp[1]*_cou[16][7][3) 
          _cou[34][8][2) = _wmp[1]*_cou[16][8][3) +  fakac2 * _cou[16][2][3)
          _cou[34][9][2) = _wmp[1]*_cou[16][9][3) 
          _cou[30][4][2) = _wmp[2]*_cou[17][4][3)  & 
               +  faka * (_cou[8][4][2) - fakaac * _cou[8][4][3) ) 
          _cou[30][5][2) = _wmp[2]*_cou[17][5][3)  & 
               +  faka * (_cou[8][5][2) - fakaac * _cou[8][5][3) ) + fakac * _cou[17][1][3)
          _cou[30][6][2) = _wmp[2]*_cou[17][6][3)  & 
               +  faka * (_cou[8][6][2) - fakaac * _cou[8][6][3) ) + fakac * _cou[17][2][3)
          _cou[30][7][2) = _wmp[2]*_cou[17][7][3)  & 
               +  faka * (_cou[8][7][2) - fakaac * _cou[8][7][3) ) 
          _cou[30][8][2) = _wmp[2]*_cou[17][8][3)  & 
               +  faka * (_cou[8][8][2) - fakaac * _cou[8][8][3) ) 
          _cou[30][9][2) = _wmp[2]*_cou[17][9][3)  & 
               +  faka * (_cou[8][9][2) - fakaac * _cou[8][9][3) ) + fakac2 * _cou[17][3][3)



          _cou[13][20][1) = _wmp[0]*_cou[4][20][2)  & 
               +  faka * (_cou[2][20][1) - fakaac * _cou[2][20][2) ) + fakac4 * _cou[4][10][2)
          _cou[14][20][1) = _wmp[1]*_cou[4][20][2)  & 
               +  faka * (_cou[1][20][1) - fakaac * _cou[1][20][2) ) 
          _cou[19][20][1) = _wmp[2]*_cou[4][20][2) 
          _cou[13][21][1) = _wmp[0]*_cou[4][21][2)  & 
               +  faka * (_cou[2][21][1) - fakaac * _cou[2][21][2) ) 
          _cou[14][21][1) = _wmp[1]*_cou[4][21][2)  & 
               +  faka * (_cou[1][21][1) - fakaac * _cou[1][21][2) ) + fakac4 * _cou[4][11][2)
          _cou[19][21][1) = _wmp[2]*_cou[4][21][2) 
          _cou[13][22][1) = _wmp[0]*_cou[4][22][2)  & 
               +  faka * (_cou[2][22][1) - fakaac * _cou[2][22][2) ) 
          _cou[14][22][1) = _wmp[1]*_cou[4][22][2)  & 
               +  faka * (_cou[1][22][1) - fakaac * _cou[1][22][2) ) 
          _cou[19][22][1) = _wmp[2]*_cou[4][22][2) +  fakac4 * _cou[4][12][2)
          _cou[13][23][1) = _wmp[0]*_cou[4][23][2)  & 
               +  faka * (_cou[2][23][1) - fakaac * _cou[2][23][2) ) + fakac3 * _cou[4][13][2)
          _cou[14][23][1) = _wmp[1]*_cou[4][23][2)  & 
               +  faka * (_cou[1][23][1) - fakaac * _cou[1][23][2) ) + fakac * _cou[4][10][2)
          _cou[19][23][1) = _wmp[2]*_cou[4][23][2) 
          _cou[13][24][1) = _wmp[0]*_cou[4][24][2)  & 
               +  faka * (_cou[2][24][1) - fakaac * _cou[2][24][2) ) + fakac * _cou[4][11][2)
          _cou[14][24][1) = _wmp[1]*_cou[4][24][2)  & 
               +  faka * (_cou[1][24][1) - fakaac * _cou[1][24][2) ) + fakac3 * _cou[4][14][2)
          _cou[19][24][1) = _wmp[2]*_cou[4][24][2) 
          _cou[13][25][1) = _wmp[0]*_cou[4][25][2)  & 
               +  faka * (_cou[2][25][1) - fakaac * _cou[2][25][2) ) + fakac3 * _cou[4][15][2)
          _cou[14][25][1) = _wmp[1]*_cou[4][25][2)  & 
               +  faka * (_cou[1][25][1) - fakaac * _cou[1][25][2) ) 
          _cou[19][25][1) = _wmp[2]*_cou[4][25][2) +  fakac * _cou[4][10][2)
          _cou[13][26][1) = _wmp[0]*_cou[4][26][2)  & 
               +  faka * (_cou[2][26][1) - fakaac * _cou[2][26][2) ) + fakac * _cou[4][12][2)
          _cou[14][26][1) = _wmp[1]*_cou[4][26][2)  & 
               +  faka * (_cou[1][26][1) - fakaac * _cou[1][26][2) ) 
          _cou[19][26][1) = _wmp[2]*_cou[4][26][2) +  fakac3 * _cou[4][16][2)
          _cou[13][27][1) = _wmp[0]*_cou[4][27][2)  & 
               +  faka * (_cou[2][27][1) - fakaac * _cou[2][27][2) ) 
          _cou[14][27][1) = _wmp[1]*_cou[4][27][2)  & 
               +  faka * (_cou[1][27][1) - fakaac * _cou[1][27][2) ) + fakac3 * _cou[4][17][2)
          _cou[19][27][1) = _wmp[2]*_cou[4][27][2) +  fakac * _cou[4][11][2)
          _cou[13][28][1) = _wmp[0]*_cou[4][28][2)  & 
               +  faka * (_cou[2][28][1) - fakaac * _cou[2][28][2) ) 
          _cou[14][28][1) = _wmp[1]*_cou[4][28][2)  & 
               +  faka * (_cou[1][28][1) - fakaac * _cou[1][28][2) ) + fakac * _cou[4][12][2)
          _cou[19][28][1) = _wmp[2]*_cou[4][28][2) +  fakac3 * _cou[4][18][2)
          _cou[13][29][1) = _wmp[0]*_cou[4][29][2)  & 
               +  faka * (_cou[2][29][1) - fakaac * _cou[2][29][2) ) + fakac2 * _cou[4][16][2)
          _cou[14][29][1) = _wmp[1]*_cou[4][29][2)  & 
               +  faka * (_cou[1][29][1) - fakaac * _cou[1][29][2) ) 
          _cou[19][29][1) = _wmp[2]*_cou[4][29][2) +  fakac2 * _cou[4][15][2)
          _cou[13][30][1) = _wmp[0]*_cou[4][30][2)  & 
               +  faka * (_cou[2][30][1) - fakaac * _cou[2][30][2) ) 
          _cou[14][30][1) = _wmp[1]*_cou[4][30][2)  & 
               +  faka * (_cou[1][30][1) - fakaac * _cou[1][30][2) ) + fakac2 * _cou[4][18][2)
          _cou[19][30][1) = _wmp[2]*_cou[4][30][2) +  fakac2 * _cou[4][17][2)
          _cou[13][31][1) = _wmp[0]*_cou[4][31][2)  & 
               +  faka * (_cou[2][31][1) - fakaac * _cou[2][31][2) ) + fakac2 * _cou[4][14][2)
          _cou[14][31][1) = _wmp[1]*_cou[4][31][2)  & 
               +  faka * (_cou[1][31][1) - fakaac * _cou[1][31][2) ) + fakac2 * _cou[4][13][2)
          _cou[19][31][1) = _wmp[2]*_cou[4][31][2) 
          _cou[13][32][1) = _wmp[0]*_cou[4][32][2)  & 
               +  faka * (_cou[2][32][1) - fakaac * _cou[2][32][2) ) + fakac2 * _cou[4][19][2)
          _cou[14][32][1) = _wmp[1]*_cou[4][32][2)  & 
               +  faka * (_cou[1][32][1) - fakaac * _cou[1][32][2) ) + fakac * _cou[4][15][2)
          _cou[19][32][1) = _wmp[2]*_cou[4][32][2) +  fakac * _cou[4][13][2)
          _cou[13][33][1) = _wmp[0]*_cou[4][33][2)  & 
               +  faka * (_cou[2][33][1) - fakaac * _cou[2][33][2) ) + fakac * _cou[4][17][2)
          _cou[14][33][1) = _wmp[1]*_cou[4][33][2)  & 
               +  faka * (_cou[1][33][1) - fakaac * _cou[1][33][2) ) + fakac2 * _cou[4][19][2)
          _cou[19][33][1) = _wmp[2]*_cou[4][33][2) +  fakac * _cou[4][14][2)
          _cou[13][34][1) = _wmp[0]*_cou[4][34][2)  & 
               +  faka * (_cou[2][34][1) - fakaac * _cou[2][34][2) ) + fakac * _cou[4][18][2)
          _cou[14][34][1) = _wmp[1]*_cou[4][34][2)  & 
               +  faka * (_cou[1][34][1) - fakaac * _cou[1][34][2) ) + fakac * _cou[4][16][2)
          _cou[19][34][1) = _wmp[2]*_cou[4][34][2) +  fakac2 * _cou[4][19][2)
          _cou[15][20][1) = _wmp[0]*_cou[5][20][2)  & 
               +  faka * (_cou[3][20][1) - fakaac * _cou[3][20][2) ) + fakac4 * _cou[5][10][2)
          _cou[16][20][1) = _wmp[2]*_cou[5][20][2)  & 
               +  faka * (_cou[1][20][1) - fakaac * _cou[1][20][2) ) 
          _cou[15][21][1) = _wmp[0]*_cou[5][21][2)  & 
               +  faka * (_cou[3][21][1) - fakaac * _cou[3][21][2) ) 
          _cou[16][21][1) = _wmp[2]*_cou[5][21][2)  & 
               +  faka * (_cou[1][21][1) - fakaac * _cou[1][21][2) ) 
          _cou[15][22][1) = _wmp[0]*_cou[5][22][2)  & 
               +  faka * (_cou[3][22][1) - fakaac * _cou[3][22][2) ) 
          _cou[16][22][1) = _wmp[2]*_cou[5][22][2)  & 
               +  faka * (_cou[1][22][1) - fakaac * _cou[1][22][2) ) + fakac4 * _cou[5][12][2)
          _cou[15][23][1) = _wmp[0]*_cou[5][23][2)  & 
               +  faka * (_cou[3][23][1) - fakaac * _cou[3][23][2) ) + fakac3 * _cou[5][13][2)
          _cou[16][23][1) = _wmp[2]*_cou[5][23][2)  & 
               +  faka * (_cou[1][23][1) - fakaac * _cou[1][23][2) ) 
          _cou[15][24][1) = _wmp[0]*_cou[5][24][2)  & 
               +  faka * (_cou[3][24][1) - fakaac * _cou[3][24][2) ) + fakac * _cou[5][11][2)
          _cou[16][24][1) = _wmp[2]*_cou[5][24][2)  & 
               +  faka * (_cou[1][24][1) - fakaac * _cou[1][24][2) ) 
          _cou[15][25][1) = _wmp[0]*_cou[5][25][2)  & 
               +  faka * (_cou[3][25][1) - fakaac * _cou[3][25][2) ) + fakac3 * _cou[5][15][2)
          _cou[16][25][1) = _wmp[2]*_cou[5][25][2)  & 
               +  faka * (_cou[1][25][1) - fakaac * _cou[1][25][2) ) + fakac * _cou[5][10][2)
          _cou[15][26][1) = _wmp[0]*_cou[5][26][2)  & 
               +  faka * (_cou[3][26][1) - fakaac * _cou[3][26][2) ) + fakac * _cou[5][12][2)
          _cou[16][26][1) = _wmp[2]*_cou[5][26][2)  & 
               +  faka * (_cou[1][26][1) - fakaac * _cou[1][26][2) ) + fakac3 * _cou[5][16][2)
          _cou[15][27][1) = _wmp[0]*_cou[5][27][2)  & 
               +  faka * (_cou[3][27][1) - fakaac * _cou[3][27][2) ) 
          _cou[16][27][1) = _wmp[2]*_cou[5][27][2)  & 
               +  faka * (_cou[1][27][1) - fakaac * _cou[1][27][2) ) + fakac * _cou[5][11][2)
          _cou[15][28][1) = _wmp[0]*_cou[5][28][2)  & 
               +  faka * (_cou[3][28][1) - fakaac * _cou[3][28][2) ) 
          _cou[16][28][1) = _wmp[2]*_cou[5][28][2)  & 
               +  faka * (_cou[1][28][1) - fakaac * _cou[1][28][2) ) + fakac3 * _cou[5][18][2)
          _cou[15][29][1) = _wmp[0]*_cou[5][29][2)  & 
               +  faka * (_cou[3][29][1) - fakaac * _cou[3][29][2) ) + fakac2 * _cou[5][16][2)
          _cou[16][29][1) = _wmp[2]*_cou[5][29][2)  & 
               +  faka * (_cou[1][29][1) - fakaac * _cou[1][29][2) ) + fakac2 * _cou[5][15][2)
          _cou[15][30][1) = _wmp[0]*_cou[5][30][2)  & 
               +  faka * (_cou[3][30][1) - fakaac * _cou[3][30][2) ) 
          _cou[16][30][1) = _wmp[2]*_cou[5][30][2)  & 
               +  faka * (_cou[1][30][1) - fakaac * _cou[1][30][2) ) + fakac2 * _cou[5][17][2)
          _cou[15][31][1) = _wmp[0]*_cou[5][31][2)  & 
               +  faka * (_cou[3][31][1) - fakaac * _cou[3][31][2) ) + fakac2 * _cou[5][14][2)
          _cou[16][31][1) = _wmp[2]*_cou[5][31][2)  & 
               +  faka * (_cou[1][31][1) - fakaac * _cou[1][31][2) ) 
          _cou[15][32][1) = _wmp[0]*_cou[5][32][2)  & 
               +  faka * (_cou[3][32][1) - fakaac * _cou[3][32][2) ) + fakac2 * _cou[5][19][2)
          _cou[16][32][1) = _wmp[2]*_cou[5][32][2)  & 
               +  faka * (_cou[1][32][1) - fakaac * _cou[1][32][2) ) + fakac * _cou[5][13][2)
          _cou[15][33][1) = _wmp[0]*_cou[5][33][2)  & 
               +  faka * (_cou[3][33][1) - fakaac * _cou[3][33][2) ) + fakac * _cou[5][17][2)
          _cou[16][33][1) = _wmp[2]*_cou[5][33][2)  & 
               +  faka * (_cou[1][33][1) - fakaac * _cou[1][33][2) ) + fakac * _cou[5][14][2)
          _cou[15][34][1) = _wmp[0]*_cou[5][34][2)  & 
               +  faka * (_cou[3][34][1) - fakaac * _cou[3][34][2) ) + fakac * _cou[5][18][2)
          _cou[16][34][1) = _wmp[2]*_cou[5][34][2)  & 
               +  faka * (_cou[1][34][1) - fakaac * _cou[1][34][2) ) + fakac2 * _cou[5][19][2)
          _cou[17][20][1) = _wmp[1]*_cou[6][20][2)  & 
               +  faka * (_cou[3][20][1) - fakaac * _cou[3][20][2) ) 
          _cou[18][20][1) = _wmp[2]*_cou[6][20][2)  & 
               +  faka * (_cou[2][20][1) - fakaac * _cou[2][20][2) ) 
          _cou[17][21][1) = _wmp[1]*_cou[6][21][2)  & 
               +  faka * (_cou[3][21][1) - fakaac * _cou[3][21][2) ) + fakac4 * _cou[6][11][2)
          _cou[18][21][1) = _wmp[2]*_cou[6][21][2)  & 
               +  faka * (_cou[2][21][1) - fakaac * _cou[2][21][2) ) 
          _cou[17][22][1) = _wmp[1]*_cou[6][22][2)  & 
               +  faka * (_cou[3][22][1) - fakaac * _cou[3][22][2) ) 
          _cou[18][22][1) = _wmp[2]*_cou[6][22][2)  & 
               +  faka * (_cou[2][22][1) - fakaac * _cou[2][22][2) ) + fakac4 * _cou[6][12][2)
          _cou[17][23][1) = _wmp[1]*_cou[6][23][2)  & 
               +  faka * (_cou[3][23][1) - fakaac * _cou[3][23][2) ) + fakac * _cou[6][10][2)
          _cou[18][23][1) = _wmp[2]*_cou[6][23][2)  & 
               +  faka * (_cou[2][23][1) - fakaac * _cou[2][23][2) ) 
          _cou[17][24][1) = _wmp[1]*_cou[6][24][2)  & 
               +  faka * (_cou[3][24][1) - fakaac * _cou[3][24][2) ) + fakac3 * _cou[6][14][2)
          _cou[18][24][1) = _wmp[2]*_cou[6][24][2)  & 
               +  faka * (_cou[2][24][1) - fakaac * _cou[2][24][2) ) 
          _cou[17][25][1) = _wmp[1]*_cou[6][25][2)  & 
               +  faka * (_cou[3][25][1) - fakaac * _cou[3][25][2) ) 
          _cou[18][25][1) = _wmp[2]*_cou[6][25][2)  & 
               +  faka * (_cou[2][25][1) - fakaac * _cou[2][25][2) ) + fakac * _cou[6][10][2)
          _cou[17][26][1) = _wmp[1]*_cou[6][26][2)  & 
               +  faka * (_cou[3][26][1) - fakaac * _cou[3][26][2) ) 
          _cou[18][26][1) = _wmp[2]*_cou[6][26][2)  & 
               +  faka * (_cou[2][26][1) - fakaac * _cou[2][26][2) ) + fakac3 * _cou[6][16][2)
          _cou[17][27][1) = _wmp[1]*_cou[6][27][2)  & 
               +  faka * (_cou[3][27][1) - fakaac * _cou[3][27][2) ) + fakac3 * _cou[6][17][2)
          _cou[18][27][1) = _wmp[2]*_cou[6][27][2)  & 
               +  faka * (_cou[2][27][1) - fakaac * _cou[2][27][2) ) + fakac * _cou[6][11][2)
          _cou[17][28][1) = _wmp[1]*_cou[6][28][2)  & 
               +  faka * (_cou[3][28][1) - fakaac * _cou[3][28][2) ) + fakac * _cou[6][12][2)
          _cou[18][28][1) = _wmp[2]*_cou[6][28][2)  & 
               +  faka * (_cou[2][28][1) - fakaac * _cou[2][28][2) ) + fakac3 * _cou[6][18][2)
          _cou[17][29][1) = _wmp[1]*_cou[6][29][2)  & 
               +  faka * (_cou[3][29][1) - fakaac * _cou[3][29][2) ) 
          _cou[18][29][1) = _wmp[2]*_cou[6][29][2)  & 
               +  faka * (_cou[2][29][1) - fakaac * _cou[2][29][2) ) + fakac2 * _cou[6][15][2)
          _cou[17][30][1) = _wmp[1]*_cou[6][30][2)  & 
               +  faka * (_cou[3][30][1) - fakaac * _cou[3][30][2) ) + fakac2 * _cou[6][18][2)
          _cou[18][30][1) = _wmp[2]*_cou[6][30][2)  & 
               +  faka * (_cou[2][30][1) - fakaac * _cou[2][30][2) ) + fakac2 * _cou[6][17][2)
          _cou[17][31][1) = _wmp[1]*_cou[6][31][2)  & 
               +  faka * (_cou[3][31][1) - fakaac * _cou[3][31][2) ) + fakac2 * _cou[6][13][2)
          _cou[18][31][1) = _wmp[2]*_cou[6][31][2)  & 
               +  faka * (_cou[2][31][1) - fakaac * _cou[2][31][2) ) 
          _cou[17][32][1) = _wmp[1]*_cou[6][32][2)  & 
               +  faka * (_cou[3][32][1) - fakaac * _cou[3][32][2) ) + fakac * _cou[6][15][2)
          _cou[18][32][1) = _wmp[2]*_cou[6][32][2)  & 
               +  faka * (_cou[2][32][1) - fakaac * _cou[2][32][2) ) + fakac * _cou[6][13][2)
          _cou[17][33][1) = _wmp[1]*_cou[6][33][2)  & 
               +  faka * (_cou[3][33][1) - fakaac * _cou[3][33][2) ) + fakac2 * _cou[6][19][2)
          _cou[18][33][1) = _wmp[2]*_cou[6][33][2)  & 
               +  faka * (_cou[2][33][1) - fakaac * _cou[2][33][2) ) + fakac * _cou[6][14][2)
          _cou[17][34][1) = _wmp[1]*_cou[6][34][2)  & 
               +  faka * (_cou[3][34][1) - fakaac * _cou[3][34][2) ) + fakac * _cou[6][16][2)
          _cou[18][34][1) = _wmp[2]*_cou[6][34][2)  & 
               +  faka * (_cou[2][34][1) - fakaac * _cou[2][34][2) ) + fakac2 * _cou[6][19][2)
          _cou[10][20][1) = _wmp[0]*_cou[7][20][2)  & 
               + faka2 * (_cou[1][20][1) - fakaac * _cou[1][20][2) ) + fakac4 * _cou[7][10][2)
          _cou[10][21][1) = _wmp[0]*_cou[7][21][2)  & 
               +  faka2 * (_cou[1][21][1) - fakaac * _cou[1][21][2) ) 
          _cou[10][22][1) = _wmp[0]*_cou[7][22][2)  & 
               +  faka2 * (_cou[1][22][1) - fakaac * _cou[1][22][2) ) 
          _cou[10][23][1) = _wmp[0]*_cou[7][23][2)  & 
               + faka2 * (_cou[1][23][1) - fakaac * _cou[1][23][2) ) + fakac3 * _cou[7][13][2)
          _cou[10][24][1) = _wmp[0]*_cou[7][24][2)  & 
               +  faka2 * (_cou[1][24][1) - fakaac * _cou[1][24][2) ) + fakac * _cou[7][11][2)
          _cou[10][25][1) = _wmp[0]*_cou[7][25][2)  & 
               + faka2 * (_cou[1][25][1) - fakaac * _cou[1][25][2) ) + fakac3 * _cou[7][15][2)
          _cou[10][26][1) = _wmp[0]*_cou[7][26][2)  & 
               +  faka2 * (_cou[1][26][1) - fakaac * _cou[1][26][2) ) + fakac * _cou[7][12][2)
          _cou[10][27][1) = _wmp[0]*_cou[7][27][2)  & 
               +  faka2 * (_cou[1][27][1) - fakaac * _cou[1][27][2) ) 
          _cou[10][28][1) = _wmp[0]*_cou[7][28][2)  & 
               +  faka2 * (_cou[1][28][1) - fakaac * _cou[1][28][2) ) 
          _cou[10][29][1) = _wmp[0]*_cou[7][29][2)  & 
               + faka2 * (_cou[1][29][1) - fakaac * _cou[1][29][2) ) + fakac2 * _cou[7][16][2)
          _cou[10][30][1) = _wmp[0]*_cou[7][30][2)  & 
               +  faka2 * (_cou[1][30][1) - fakaac * _cou[1][30][2) ) 
          _cou[10][31][1) = _wmp[0]*_cou[7][31][2)  & 
               + faka2 * (_cou[1][31][1) - fakaac * _cou[1][31][2) ) + fakac2 * _cou[7][14][2)
          _cou[10][32][1) = _wmp[0]*_cou[7][32][2)  & 
               + faka2 * (_cou[1][32][1) - fakaac * _cou[1][32][2) ) + fakac2 * _cou[7][19][2)
          _cou[10][33][1) = _wmp[0]*_cou[7][33][2)  & 
               +  faka2 * (_cou[1][33][1) - fakaac * _cou[1][33][2) ) + fakac * _cou[7][17][2)
          _cou[10][34][1) = _wmp[0]*_cou[7][34][2)  & 
               +  faka2 * (_cou[1][34][1) - fakaac * _cou[1][34][2) ) + fakac * _cou[7][18][2)
          _cou[11][20][1) = _wmp[1]*_cou[8][20][2)  & 
               +  faka2 * (_cou[2][20][1) - fakaac * _cou[2][20][2) ) 
          _cou[11][21][1) = _wmp[1]*_cou[8][21][2)  & 
               + faka2 * (_cou[2][21][1) - fakaac * _cou[2][21][2) ) + fakac4 * _cou[8][11][2)
          _cou[11][22][1) = _wmp[1]*_cou[8][22][2)  & 
               +  faka2 * (_cou[2][22][1) - fakaac * _cou[2][22][2) ) 
          _cou[11][23][1) = _wmp[1]*_cou[8][23][2)  & 
               +  faka2 * (_cou[2][23][1) - fakaac * _cou[2][23][2) ) + fakac * _cou[8][10][2)
          _cou[11][24][1) = _wmp[1]*_cou[8][24][2)  & 
               + faka2 * (_cou[2][24][1) - fakaac * _cou[2][24][2) ) + fakac3 * _cou[8][14][2)
          _cou[11][25][1) = _wmp[1]*_cou[8][25][2)  & 
               +  faka2 * (_cou[2][25][1) - fakaac * _cou[2][25][2) ) 
          _cou[11][26][1) = _wmp[1]*_cou[8][26][2)  & 
               +  faka2 * (_cou[2][26][1) - fakaac * _cou[2][26][2) ) 
          _cou[11][27][1) = _wmp[1]*_cou[8][27][2)  & 
               + faka2 * (_cou[2][27][1) - fakaac * _cou[2][27][2) ) + fakac3 * _cou[8][17][2)
          _cou[11][28][1) = _wmp[1]*_cou[8][28][2)  & 
               +  faka2 * (_cou[2][28][1) - fakaac * _cou[2][28][2) ) + fakac * _cou[8][12][2)
          _cou[11][29][1) = _wmp[1]*_cou[8][29][2)  & 
               +  faka2 * (_cou[2][29][1) - fakaac * _cou[2][29][2) ) 
          _cou[11][30][1) = _wmp[1]*_cou[8][30][2)  & 
               + faka2 * (_cou[2][30][1) - fakaac * _cou[2][30][2) ) + fakac2 * _cou[8][18][2)
          _cou[11][31][1) = _wmp[1]*_cou[8][31][2)  & 
               + faka2 * (_cou[2][31][1) - fakaac * _cou[2][31][2) ) + fakac2 * _cou[8][13][2)
          _cou[11][32][1) = _wmp[1]*_cou[8][32][2)  & 
               +  faka2 * (_cou[2][32][1) - fakaac * _cou[2][32][2) ) + fakac * _cou[8][15][2)
          _cou[11][33][1) = _wmp[1]*_cou[8][33][2)  & 
               + faka2 * (_cou[2][33][1) - fakaac * _cou[2][33][2) ) + fakac2 * _cou[8][19][2)
          _cou[11][34][1) = _wmp[1]*_cou[8][34][2)  & 
               +  faka2 * (_cou[2][34][1) - fakaac * _cou[2][34][2) ) + fakac * _cou[8][16][2)
          _cou[12][20][1) = _wmp[2]*_cou[9][20][2)  & 
               +  faka2 * (_cou[3][20][1) - fakaac * _cou[3][20][2) ) 
          _cou[12][21][1) = _wmp[2]*_cou[9][21][2)  & 
               +  faka2 * (_cou[3][21][1) - fakaac * _cou[3][21][2) ) 
          _cou[12][22][1) = _wmp[2]*_cou[9][22][2)  & 
               + faka2 * (_cou[3][22][1) - fakaac * _cou[3][22][2) ) + fakac4 * _cou[9][12][2)
          _cou[12][23][1) = _wmp[2]*_cou[9][23][2)  & 
               +  faka2 * (_cou[3][23][1) - fakaac * _cou[3][23][2) ) 
          _cou[12][24][1) = _wmp[2]*_cou[9][24][2)  & 
               +  faka2 * (_cou[3][24][1) - fakaac * _cou[3][24][2) ) 
          _cou[12][25][1) = _wmp[2]*_cou[9][25][2)  & 
               +  faka2 * (_cou[3][25][1) - fakaac * _cou[3][25][2) ) + fakac * _cou[9][10][2)
          _cou[12][26][1) = _wmp[2]*_cou[9][26][2)  & 
               + faka2 * (_cou[3][26][1) - fakaac * _cou[3][26][2) ) + fakac3 * _cou[9][16][2)
          _cou[12][27][1) = _wmp[2]*_cou[9][27][2)  & 
               +  faka2 * (_cou[3][27][1) - fakaac * _cou[3][27][2) ) + fakac * _cou[9][11][2)
          _cou[12][28][1) = _wmp[2]*_cou[9][28][2)  & 
               + faka2 * (_cou[3][28][1) - fakaac * _cou[3][28][2) ) + fakac3 * _cou[9][18][2)
          _cou[12][29][1) = _wmp[2]*_cou[9][29][2)  & 
               + faka2 * (_cou[3][29][1) - fakaac * _cou[3][29][2) ) + fakac2 * _cou[9][15][2)
          _cou[12][30][1) = _wmp[2]*_cou[9][30][2)  & 
               + faka2 * (_cou[3][30][1) - fakaac * _cou[3][30][2) ) + fakac2 * _cou[9][17][2)
          _cou[12][31][1) = _wmp[2]*_cou[9][31][2)  & 
               +  faka2 * (_cou[3][31][1) - fakaac * _cou[3][31][2) ) 
          _cou[12][32][1) = _wmp[2]*_cou[9][32][2)  & 
               +  faka2 * (_cou[3][32][1) - fakaac * _cou[3][32][2) ) + fakac * _cou[9][13][2)
          _cou[12][33][1) = _wmp[2]*_cou[9][33][2)  & 
               +  faka2 * (_cou[3][33][1) - fakaac * _cou[3][33][2) ) + fakac * _cou[9][14][2)
          _cou[12][34][1) = _wmp[2]*_cou[9][34][2)  & 
               + faka2 * (_cou[3][34][1) - fakaac * _cou[3][34][2) ) + fakac2 * _cou[9][19][2)

          _cou[20][10][1) = _wmp[0]*_cou[10][10][2)  & 
               + faka3 * (_cou[7][10][1) - fakaac * _cou[7][10][2) ) + fakac3 * _cou[10][7][2)
          _cou[23][10][1) = _wmp[1]*_cou[10][10][2) 
          _cou[25][10][1) = _wmp[2]*_cou[10][10][2) 
          _cou[20][11][1) = _wmp[0]*_cou[10][11][2)  & 
               +  faka3 * (_cou[7][11][1) - fakaac * _cou[7][11][2) )
          _cou[23][11][1) = _wmp[1]*_cou[10][11][2) +  fakac3 * _cou[10][8][2)
          _cou[25][11][1) = _wmp[2]*_cou[10][11][2) 
          _cou[20][12][1) = _wmp[0]*_cou[10][12][2)  & 
               +  faka3 * (_cou[7][12][1) - fakaac * _cou[7][12][2) ) 
          _cou[23][12][1) = _wmp[1]*_cou[10][12][2) 
          _cou[25][12][1) = _wmp[2]*_cou[10][12][2) +  fakac3 * _cou[10][9][2)
          _cou[20][13][1) = _wmp[0]*_cou[10][13][2)  & 
               + faka3 * (_cou[7][13][1) - fakaac * _cou[7][13][2) ) + fakac2 * _cou[10][4][2)
          _cou[23][13][1) = _wmp[1]*_cou[10][13][2) +  fakac * _cou[10][7][2)
          _cou[25][13][1) = _wmp[2]*_cou[10][13][2) 
          _cou[20][14][1) = _wmp[0]*_cou[10][14][2)  & 
               +  faka3 * (_cou[7][14][1) - fakaac * _cou[7][14][2) ) + fakac * _cou[10][8][2)
          _cou[23][14][1) = _wmp[1]*_cou[10][14][2) +  fakac2 * _cou[10][4][2)
          _cou[25][14][1) = _wmp[2]*_cou[10][14][2) 
          _cou[20][15][1) = _wmp[0]*_cou[10][15][2)  & 
               + faka3 * (_cou[7][15][1) - fakaac * _cou[7][15][2) ) + fakac2 * _cou[10][5][2)
          _cou[23][15][1) = _wmp[1]*_cou[10][15][2) 
          _cou[25][15][1) = _wmp[2]*_cou[10][15][2) +  fakac * _cou[10][7][2)
          _cou[20][16][1) = _wmp[0]*_cou[10][16][2)  & 
               +  faka3 * (_cou[7][16][1) - fakaac * _cou[7][16][2) ) + fakac * _cou[10][9][2)
          _cou[23][16][1) = _wmp[1]*_cou[10][16][2) 
          _cou[25][16][1) = _wmp[2]*_cou[10][16][2) +  fakac2 * _cou[10][5][2)
          _cou[20][17][1) = _wmp[0]*_cou[10][17][2)  & 
               +  faka3 * (_cou[7][17][1) - fakaac * _cou[7][17][2) ) 
          _cou[23][17][1) = _wmp[1]*_cou[10][17][2) +  fakac2 * _cou[10][6][2)
          _cou[25][17][1) = _wmp[2]*_cou[10][17][2) +  fakac * _cou[10][8][2)
          _cou[20][18][1) = _wmp[0]*_cou[10][18][2)  & 
               +  faka3 * (_cou[7][18][1) - fakaac * _cou[7][18][2) ) 
          _cou[23][18][1) = _wmp[1]*_cou[10][18][2) +  fakac * _cou[10][9][2)
          _cou[25][18][1) = _wmp[2]*_cou[10][18][2) +  fakac2 * _cou[10][6][2)
          _cou[20][19][1) = _wmp[0]*_cou[10][19][2)  & 
               +  faka3 * (_cou[7][19][1) - fakaac * _cou[7][19][2) ) + fakac * _cou[10][6][2)
          _cou[23][19][1) = _wmp[1]*_cou[10][19][2) +  fakac * _cou[10][5][2)
          _cou[25][19][1) = _wmp[2]*_cou[10][19][2) +  fakac * _cou[10][4][2)
          _cou[24][10][1) = _wmp[0]*_cou[11][10][2) +  fakac3 * _cou[11][7][2)
          _cou[21][10][1) = _wmp[1]*_cou[11][10][2)  & 
               +  faka3 * (_cou[8][10][1) - fakaac * _cou[8][10][2) ) 
          _cou[27][10][1) = _wmp[2]*_cou[11][10][2) 
          _cou[24][11][1) = _wmp[0]*_cou[11][11][2) 
          _cou[21][11][1) = _wmp[1]*_cou[11][11][2)  & 
               + faka3 * (_cou[8][11][1) - fakaac * _cou[8][11][2) ) + fakac3 * _cou[11][8][2)
          _cou[27][11][1) = _wmp[2]*_cou[11][11][2) 
          _cou[24][12][1) = _wmp[0]*_cou[11][12][2) 
          _cou[21][12][1) = _wmp[1]*_cou[11][12][2)  & 
               +  faka3 * (_cou[8][12][1) - fakaac * _cou[8][12][2) ) 
          _cou[27][12][1) = _wmp[2]*_cou[11][12][2) +  fakac3 * _cou[11][9][2)
          _cou[24][13][1) = _wmp[0]*_cou[11][13][2) +  fakac2 * _cou[11][4][2)
          _cou[21][13][1) = _wmp[1]*_cou[11][13][2)  & 
               +  faka3 * (_cou[8][13][1) - fakaac * _cou[8][13][2) ) + fakac * _cou[11][7][2)
          _cou[27][13][1) = _wmp[2]*_cou[11][13][2) 
          _cou[24][14][1) = _wmp[0]*_cou[11][14][2) +  fakac * _cou[11][8][2)
          _cou[21][14][1) = _wmp[1]*_cou[11][14][2)  & 
               + faka3 * (_cou[8][14][1) - fakaac * _cou[8][14][2) ) + fakac2 * _cou[11][4][2)
          _cou[27][14][1) = _wmp[2]*_cou[11][14][2) 
          _cou[24][15][1) = _wmp[0]*_cou[11][15][2) +  fakac2 * _cou[11][5][2)
          _cou[21][15][1) = _wmp[1]*_cou[11][15][2)  & 
               +  faka3 * (_cou[8][15][1) - fakaac * _cou[8][15][2) ) 
          _cou[27][15][1) = _wmp[2]*_cou[11][15][2) +  fakac * _cou[11][7][2)
          _cou[24][16][1) = _wmp[0]*_cou[11][16][2) +  fakac * _cou[11][9][2)
          _cou[21][16][1) = _wmp[1]*_cou[11][16][2)  & 
               +  faka3 * (_cou[8][16][1) - fakaac * _cou[8][16][2) ) 
          _cou[27][16][1) = _wmp[2]*_cou[11][16][2) +  fakac2 * _cou[11][5][2)
          _cou[24][17][1) = _wmp[0]*_cou[11][17][2) 
          _cou[21][17][1) = _wmp[1]*_cou[11][17][2)  & 
               + faka3 * (_cou[8][17][1) - fakaac * _cou[8][17][2) ) + fakac2 * _cou[11][6][2)
          _cou[27][17][1) = _wmp[2]*_cou[11][17][2) +  fakac * _cou[11][8][2)
          _cou[24][18][1) = _wmp[0]*_cou[11][18][2) 
          _cou[21][18][1) = _wmp[1]*_cou[11][18][2)  & 
               +  faka3 * (_cou[8][18][1) - fakaac * _cou[8][18][2) ) + fakac * _cou[11][9][2)
          _cou[27][18][1) = _wmp[2]*_cou[11][18][2) +  fakac2 * _cou[11][6][2)
          _cou[24][19][1) = _wmp[0]*_cou[11][19][2) +  fakac * _cou[11][6][2)
          _cou[21][19][1) = _wmp[1]*_cou[11][19][2)  & 
               +  faka3 * (_cou[8][19][1) - fakaac * _cou[8][19][2) ) + fakac * _cou[11][5][2)
          _cou[27][19][1) = _wmp[2]*_cou[11][19][2) +  fakac * _cou[11][4][2)
          _cou[26][10][1) = _wmp[0]*_cou[12][10][2) +  fakac3 * _cou[12][7][2)
          _cou[28][10][1) = _wmp[1]*_cou[12][10][2) 
          _cou[22][10][1) = _wmp[2]*_cou[12][10][2)  & 
               +  faka3 * (_cou[9][10][1) - fakaac * _cou[9][10][2) ) 
          _cou[26][11][1) = _wmp[0]*_cou[12][11][2) 
          _cou[28][11][1) = _wmp[1]*_cou[12][11][2) +  fakac3 * _cou[12][8][2)
          _cou[22][11][1) = _wmp[2]*_cou[12][11][2)  & 
               +  faka3 * (_cou[9][11][1) - fakaac * _cou[9][11][2) ) 
          _cou[26][12][1) = _wmp[0]*_cou[12][12][2) 
          _cou[28][12][1) = _wmp[1]*_cou[12][12][2) 
          _cou[22][12][1) = _wmp[2]*_cou[12][12][2)  & 
               + faka3 * (_cou[9][12][1) - fakaac * _cou[9][12][2) ) + fakac3 * _cou[12][9][2)
          _cou[26][13][1) = _wmp[0]*_cou[12][13][2) +  fakac2 * _cou[12][4][2)
          _cou[28][13][1) = _wmp[1]*_cou[12][13][2) +  fakac * _cou[12][7][2)
          _cou[22][13][1) = _wmp[2]*_cou[12][13][2)  & 
               +  faka3 * (_cou[9][13][1) - fakaac * _cou[9][13][2) ) 
          _cou[26][14][1) = _wmp[0]*_cou[12][14][2) +  fakac * _cou[12][8][2)
          _cou[28][14][1) = _wmp[1]*_cou[12][14][2) +  fakac2 * _cou[12][4][2)
          _cou[22][14][1) = _wmp[2]*_cou[12][14][2)  & 
               +  faka3 * (_cou[9][14][1) - fakaac * _cou[9][14][2) ) 
          _cou[26][15][1) = _wmp[0]*_cou[12][15][2) +  fakac2 * _cou[12][5][2)
          _cou[28][15][1) = _wmp[1]*_cou[12][15][2) 
          _cou[22][15][1) = _wmp[2]*_cou[12][15][2)  & 
               +  faka3 * (_cou[9][15][1) - fakaac * _cou[9][15][2) ) + fakac * _cou[12][7][2)
          _cou[26][16][1) = _wmp[0]*_cou[12][16][2) +  fakac * _cou[12][9][2)
          _cou[28][16][1) = _wmp[1]*_cou[12][16][2) 
          _cou[22][16][1) = _wmp[2]*_cou[12][16][2)  & 
               + faka3 * (_cou[9][16][1) - fakaac * _cou[9][16][2) ) + fakac2 * _cou[12][5][2)
          _cou[26][17][1) = _wmp[0]*_cou[12][17][2) 
          _cou[28][17][1) = _wmp[1]*_cou[12][17][2) +  fakac2 * _cou[12][6][2)
          _cou[22][17][1) = _wmp[2]*_cou[12][17][2)  & 
               +  faka3 * (_cou[9][17][1) - fakaac * _cou[9][17][2) ) + fakac * _cou[12][8][2)
          _cou[26][18][1) = _wmp[0]*_cou[12][18][2) 
          _cou[28][18][1) = _wmp[1]*_cou[12][18][2) +  fakac * _cou[12][9][2)
          _cou[22][18][1) = _wmp[2]*_cou[12][18][2)  & 
               + faka3 * (_cou[9][18][1) - fakaac * _cou[9][18][2) ) + fakac2 * _cou[12][6][2)
          _cou[26][19][1) = _wmp[0]*_cou[12][19][2) +  fakac * _cou[12][6][2)
          _cou[28][19][1) = _wmp[1]*_cou[12][19][2) +  fakac * _cou[12][5][2)
          _cou[22][19][1) = _wmp[2]*_cou[12][19][2)  & 
               +  faka3 * (_cou[9][19][1) - fakaac * _cou[9][19][2) ) + fakac * _cou[12][4][2)
          _cou[31][10][1) = _wmp[1]*_cou[13][10][2)  & 
               +  faka * (_cou[7][10][1) - fakaac * _cou[7][10][2) ) 
          _cou[32][10][1) = _wmp[2]*_cou[13][10][2) 
          _cou[31][11][1) = _wmp[1]*_cou[13][11][2)  & 
               +  faka * (_cou[7][11][1) - fakaac * _cou[7][11][2) ) + fakac3 * _cou[13][8][2)
          _cou[32][11][1) = _wmp[2]*_cou[13][11][2) 
          _cou[31][12][1) = _wmp[1]*_cou[13][12][2)  & 
               +  faka * (_cou[7][12][1) - fakaac * _cou[7][12][2) ) 
          _cou[32][12][1) = _wmp[2]*_cou[13][12][2) +  fakac3 * _cou[13][9][2)
          _cou[31][13][1) = _wmp[1]*_cou[13][13][2)  & 
               +  faka * (_cou[7][13][1) - fakaac * _cou[7][13][2) ) + fakac * _cou[13][7][2)
          _cou[32][13][1) = _wmp[2]*_cou[13][13][2) 
          _cou[31][14][1) = _wmp[1]*_cou[13][14][2)  & 
               +  faka * (_cou[7][14][1) - fakaac * _cou[7][14][2) ) + fakac2 * _cou[13][4][2)
          _cou[32][14][1) = _wmp[2]*_cou[13][14][2) 
          _cou[31][15][1) = _wmp[1]*_cou[13][15][2)  & 
               +  faka * (_cou[7][15][1) - fakaac * _cou[7][15][2) ) 
          _cou[32][15][1) = _wmp[2]*_cou[13][15][2) +  fakac * _cou[13][7][2)
          _cou[31][16][1) = _wmp[1]*_cou[13][16][2)  & 
               +  faka * (_cou[7][16][1) - fakaac * _cou[7][16][2) ) 
          _cou[32][16][1) = _wmp[2]*_cou[13][16][2) +  fakac2 * _cou[13][5][2)
          _cou[31][17][1) = _wmp[1]*_cou[13][17][2)  & 
               +  faka * (_cou[7][17][1) - fakaac * _cou[7][17][2) ) + fakac2 * _cou[13][6][2)
          _cou[32][17][1) = _wmp[2]*_cou[13][17][2) +  fakac * _cou[13][8][2)
          _cou[31][18][1) = _wmp[1]*_cou[13][18][2)  & 
               +  faka * (_cou[7][18][1) - fakaac * _cou[7][18][2) ) + fakac * _cou[13][9][2)
          _cou[32][18][1) = _wmp[2]*_cou[13][18][2) +  fakac2 * _cou[13][6][2)
          _cou[31][19][1) = _wmp[1]*_cou[13][19][2)  & 
               +  faka * (_cou[7][19][1) - fakaac * _cou[7][19][2) ) + fakac * _cou[13][5][2)
          _cou[32][19][1) = _wmp[2]*_cou[13][19][2) +  fakac * _cou[13][4][2)
          _cou[33][10][1) = _wmp[2]*_cou[14][10][2) 
          _cou[33][11][1) = _wmp[2]*_cou[14][11][2) 
          _cou[33][12][1) = _wmp[2]*_cou[14][12][2) +  fakac3 * _cou[14][9][2)
          _cou[33][13][1) = _wmp[2]*_cou[14][13][2) 
          _cou[33][14][1) = _wmp[2]*_cou[14][14][2) 
          _cou[33][15][1) = _wmp[2]*_cou[14][15][2) +  fakac * _cou[14][7][2)
          _cou[33][16][1) = _wmp[2]*_cou[14][16][2) +  fakac2 * _cou[14][5][2)
          _cou[33][17][1) = _wmp[2]*_cou[14][17][2) +  fakac * _cou[14][8][2)
          _cou[33][18][1) = _wmp[2]*_cou[14][18][2) +  fakac2 * _cou[14][6][2)
          _cou[33][19][1) = _wmp[2]*_cou[14][19][2) +  fakac * _cou[14][4][2)
          _cou[29][10][1) = _wmp[2]*_cou[15][10][2)  & 
               +  faka * (_cou[7][10][1) - fakaac * _cou[7][10][2) ) 
          _cou[29][11][1) = _wmp[2]*_cou[15][11][2)  & 
               +  faka * (_cou[7][11][1) - fakaac * _cou[7][11][2) ) 
          _cou[29][12][1) = _wmp[2]*_cou[15][12][2)  & 
               +  faka * (_cou[7][12][1) - fakaac * _cou[7][12][2) ) + fakac3 * _cou[15][9][2)
          _cou[29][13][1) = _wmp[2]*_cou[15][13][2)  & 
               +  faka * (_cou[7][13][1) - fakaac * _cou[7][13][2) ) 
          _cou[29][14][1) = _wmp[2]*_cou[15][14][2)  & 
               +  faka * (_cou[7][14][1) - fakaac * _cou[7][14][2) ) 
          _cou[29][15][1) = _wmp[2]*_cou[15][15][2)  & 
               +  faka * (_cou[7][15][1) - fakaac * _cou[7][15][2) ) + fakac * _cou[15][7][2)
          _cou[29][16][1) = _wmp[2]*_cou[15][16][2)  & 
               +  faka * (_cou[7][16][1) - fakaac * _cou[7][16][2) ) + fakac2 * _cou[15][5][2)
          _cou[29][17][1) = _wmp[2]*_cou[15][17][2)  & 
               +  faka * (_cou[7][17][1) - fakaac * _cou[7][17][2) ) + fakac * _cou[15][8][2)
          _cou[29][18][1) = _wmp[2]*_cou[15][18][2)  & 
               +  faka * (_cou[7][18][1) - fakaac * _cou[7][18][2) ) + fakac2 * _cou[15][6][2)
          _cou[29][19][1) = _wmp[2]*_cou[15][19][2)  & 
               +  faka * (_cou[7][19][1) - fakaac * _cou[7][19][2) ) + fakac * _cou[15][4][2)
          _cou[34][10][1) = _wmp[1]*_cou[16][10][2) 
          _cou[34][11][1) = _wmp[1]*_cou[16][11][2) +  fakac3 * _cou[16][8][2)
          _cou[34][12][1) = _wmp[1]*_cou[16][12][2) 
          _cou[34][13][1) = _wmp[1]*_cou[16][13][2) +  fakac * _cou[16][7][2)
          _cou[34][14][1) = _wmp[1]*_cou[16][14][2) +  fakac2 * _cou[16][4][2)
          _cou[34][15][1) = _wmp[1]*_cou[16][15][2) 
          _cou[34][16][1) = _wmp[1]*_cou[16][16][2) 
          _cou[34][17][1) = _wmp[1]*_cou[16][17][2) +  fakac2 * _cou[16][6][2)
          _cou[34][18][1) = _wmp[1]*_cou[16][18][2) +  fakac * _cou[16][9][2)
          _cou[34][19][1) = _wmp[1]*_cou[16][19][2) +  fakac * _cou[16][5][2)
          _cou[30][10][1) = _wmp[2]*_cou[17][10][2)  & 
               +  faka * (_cou[8][10][1) - fakaac * _cou[8][10][2) ) 
          _cou[30][11][1) = _wmp[2]*_cou[17][11][2)  & 
               +  faka * (_cou[8][11][1) - fakaac * _cou[8][11][2) ) 
          _cou[30][12][1) = _wmp[2]*_cou[17][12][2)  & 
               +  faka * (_cou[8][12][1) - fakaac * _cou[8][12][2) ) + fakac3 * _cou[17][9][2)
          _cou[30][13][1) = _wmp[2]*_cou[17][13][2)  & 
               +  faka * (_cou[8][13][1) - fakaac * _cou[8][13][2) ) 
          _cou[30][14][1) = _wmp[2]*_cou[17][14][2)  & 
               +  faka * (_cou[8][14][1) - fakaac * _cou[8][14][2) ) 
          _cou[30][15][1) = _wmp[2]*_cou[17][15][2)  & 
               +  faka * (_cou[8][15][1) - fakaac * _cou[8][15][2) ) + fakac * _cou[17][7][2)
          _cou[30][16][1) = _wmp[2]*_cou[17][16][2)  & 
               +  faka * (_cou[8][16][1) - fakaac * _cou[8][16][2) ) + fakac2 * _cou[17][5][2)
          _cou[30][17][1) = _wmp[2]*_cou[17][17][2)  & 
               +  faka * (_cou[8][17][1) - fakaac * _cou[8][17][2) ) + fakac * _cou[17][8][2)
          _cou[30][18][1) = _wmp[2]*_cou[17][18][2)  & 
               +  faka * (_cou[8][18][1) - fakaac * _cou[8][18][2) ) + fakac2 * _cou[17][6][2)
          _cou[30][19][1) = _wmp[2]*_cou[17][19][2)  & 
               +  faka * (_cou[8][19][1) - fakaac * _cou[8][19][2) ) + fakac * _cou[17][4][2)


          _cou[20][20][0) = _wmp[0]*_cou[10][20][1)  & 
               + faka3 * (_cou[7][20][0) - fakaac * _cou[7][20][1) ) + fakac4 * _cou[10][10][1)
          _cou[23][20][0) = _wmp[1]*_cou[10][20][1) 
          _cou[25][20][0) = _wmp[2]*_cou[10][20][1) 
          _cou[20][21][0) = _wmp[0]*_cou[10][21][1)  & 
               +  faka3 * (_cou[7][21][0) - fakaac * _cou[7][21][1) )
          _cou[23][21][0) = _wmp[1]*_cou[10][21][1) +  fakac4 * _cou[10][11][1)
          _cou[25][21][0) = _wmp[2]*_cou[10][21][1) 
          _cou[20][22][0) = _wmp[0]*_cou[10][22][1)  & 
               +  faka3 * (_cou[7][22][0) - fakaac * _cou[7][22][1) )
          _cou[23][22][0) = _wmp[1]*_cou[10][22][1) 
          _cou[25][22][0) = _wmp[2]*_cou[10][22][1) +  fakac4 * _cou[10][12][1)
          _cou[20][23][0) = _wmp[0]*_cou[10][23][1)  & 
               + faka3 * (_cou[7][23][0) - fakaac * _cou[7][23][1) ) + fakac3 * _cou[10][13][1)
          _cou[23][23][0) = _wmp[1]*_cou[10][23][1) +  fakac * _cou[10][10][1)
          _cou[25][23][0) = _wmp[2]*_cou[10][23][1) 
          _cou[20][24][0) = _wmp[0]*_cou[10][24][1)  & 
               +  faka3 * (_cou[7][24][0) - fakaac * _cou[7][24][1) ) + fakac * _cou[10][11][1)
          _cou[23][24][0) = _wmp[1]*_cou[10][24][1) +  fakac3 * _cou[10][14][1)
          _cou[25][24][0) = _wmp[2]*_cou[10][24][1) 
          _cou[20][25][0) = _wmp[0]*_cou[10][25][1)  & 
               + faka3 * (_cou[7][25][0) - fakaac * _cou[7][25][1) ) + fakac3 * _cou[10][15][1)
          _cou[23][25][0) = _wmp[1]*_cou[10][25][1) 
          _cou[25][25][0) = _wmp[2]*_cou[10][25][1) +  fakac * _cou[10][10][1)
          _cou[20][26][0) = _wmp[0]*_cou[10][26][1)  & 
               +  faka3 * (_cou[7][26][0) - fakaac * _cou[7][26][1) ) + fakac * _cou[10][12][1)
          _cou[23][26][0) = _wmp[1]*_cou[10][26][1) 
          _cou[25][26][0) = _wmp[2]*_cou[10][26][1) +  fakac3 * _cou[10][16][1)
          _cou[20][27][0) = _wmp[0]*_cou[10][27][1)  & 
               +  faka3 * (_cou[7][27][0) - fakaac * _cou[7][27][1) )
          _cou[23][27][0) = _wmp[1]*_cou[10][27][1) +  fakac3 * _cou[10][17][1)
          _cou[25][27][0) = _wmp[2]*_cou[10][27][1) +  fakac * _cou[10][11][1)
          _cou[20][28][0) = _wmp[0]*_cou[10][28][1)  & 
               +  faka3 * (_cou[7][28][0) - fakaac * _cou[7][28][1) )
          _cou[23][28][0) = _wmp[1]*_cou[10][28][1) +  fakac * _cou[10][12][1)
          _cou[25][28][0) = _wmp[2]*_cou[10][28][1) +  fakac3 * _cou[10][18][1)
          _cou[20][29][0) = _wmp[0]*_cou[10][29][1)  & 
               + faka3 * (_cou[7][29][0) - fakaac * _cou[7][29][1) ) + fakac2 * _cou[10][16][1)
          _cou[23][29][0) = _wmp[1]*_cou[10][29][1) 
          _cou[25][29][0) = _wmp[2]*_cou[10][29][1) +  fakac2 * _cou[10][15][1)
          _cou[20][30][0) = _wmp[0]*_cou[10][30][1)  & 
               +  faka3 * (_cou[7][30][0) - fakaac * _cou[7][30][1) )
          _cou[23][30][0) = _wmp[1]*_cou[10][30][1) +  fakac2 * _cou[10][18][1)
          _cou[25][30][0) = _wmp[2]*_cou[10][30][1) +  fakac2 * _cou[10][17][1)
          _cou[20][31][0) = _wmp[0]*_cou[10][31][1)  & 
               + faka3 * (_cou[7][31][0) - fakaac * _cou[7][31][1) ) + fakac2 * _cou[10][14][1)
          _cou[23][31][0) = _wmp[1]*_cou[10][31][1) +  fakac2 * _cou[10][13][1)
          _cou[25][31][0) = _wmp[2]*_cou[10][31][1) 
          _cou[20][32][0) = _wmp[0]*_cou[10][32][1)  & 
               + faka3 * (_cou[7][32][0) - fakaac * _cou[7][32][1) ) + fakac2 * _cou[10][19][1)
          _cou[23][32][0) = _wmp[1]*_cou[10][32][1) +  fakac * _cou[10][15][1)
          _cou[25][32][0) = _wmp[2]*_cou[10][32][1) +  fakac * _cou[10][13][1)
          _cou[20][33][0) = _wmp[0]*_cou[10][33][1)  & 
               +  faka3 * (_cou[7][33][0) - fakaac * _cou[7][33][1) ) + fakac * _cou[10][17][1)
          _cou[23][33][0) = _wmp[1]*_cou[10][33][1) +  fakac2 * _cou[10][19][1)
          _cou[25][33][0) = _wmp[2]*_cou[10][33][1) +  fakac * _cou[10][14][1)
          _cou[20][34][0) = _wmp[0]*_cou[10][34][1)  & 
               +  faka3 * (_cou[7][34][0) - fakaac * _cou[7][34][1) ) + fakac * _cou[10][18][1)
          _cou[23][34][0) = _wmp[1]*_cou[10][34][1) +  fakac * _cou[10][16][1)
          _cou[25][34][0) = _wmp[2]*_cou[10][34][1) +  fakac2 * _cou[10][19][1)
          _cou[24][20][0) = _wmp[0]*_cou[11][20][1) +  fakac4 * _cou[11][10][1)
          _cou[21][20][0) = _wmp[1]*_cou[11][20][1)  & 
               +  faka3 * (_cou[8][20][0) - fakaac * _cou[8][20][1) )
          _cou[27][20][0) = _wmp[2]*_cou[11][20][1) 
          _cou[24][21][0) = _wmp[0]*_cou[11][21][1) 
          _cou[21][21][0) = _wmp[1]*_cou[11][21][1)  & 
               + faka3 * (_cou[8][21][0) - fakaac * _cou[8][21][1) ) + fakac4 * _cou[11][11][1)
          _cou[27][21][0) = _wmp[2]*_cou[11][21][1) 
          _cou[24][22][0) = _wmp[0]*_cou[11][22][1) 
          _cou[21][22][0) = _wmp[1]*_cou[11][22][1)  & 
               +  faka3 * (_cou[8][22][0) - fakaac * _cou[8][22][1) )
          _cou[27][22][0) = _wmp[2]*_cou[11][22][1) +  fakac4 * _cou[11][12][1)
          _cou[24][23][0) = _wmp[0]*_cou[11][23][1) +  fakac3 * _cou[11][13][1)
          _cou[21][23][0) = _wmp[1]*_cou[11][23][1)  & 
               +  faka3 * (_cou[8][23][0) - fakaac * _cou[8][23][1) ) + fakac * _cou[11][10][1)
          _cou[27][23][0) = _wmp[2]*_cou[11][23][1) 
          _cou[24][24][0) = _wmp[0]*_cou[11][24][1) +  fakac * _cou[11][11][1)
          _cou[21][24][0) = _wmp[1]*_cou[11][24][1)  & 
               + faka3 * (_cou[8][24][0) - fakaac * _cou[8][24][1) ) + fakac3 * _cou[11][14][1)
          _cou[27][24][0) = _wmp[2]*_cou[11][24][1) 
          _cou[24][25][0) = _wmp[0]*_cou[11][25][1) +  fakac3 * _cou[11][15][1)
          _cou[21][25][0) = _wmp[1]*_cou[11][25][1)  & 
               +  faka3 * (_cou[8][25][0) - fakaac * _cou[8][25][1) ) 
          _cou[27][25][0) = _wmp[2]*_cou[11][25][1) +  fakac * _cou[11][10][1)
          _cou[24][26][0) = _wmp[0]*_cou[11][26][1) +  fakac * _cou[11][12][1)
          _cou[21][26][0) = _wmp[1]*_cou[11][26][1)  & 
               +  faka3 * (_cou[8][26][0) - fakaac * _cou[8][26][1) ) 
          _cou[27][26][0) = _wmp[2]*_cou[11][26][1) +  fakac3 * _cou[11][16][1)
          _cou[24][27][0) = _wmp[0]*_cou[11][27][1) 
          _cou[21][27][0) = _wmp[1]*_cou[11][27][1)  & 
               + faka3 * (_cou[8][27][0) - fakaac * _cou[8][27][1) ) + fakac3 * _cou[11][17][1)
          _cou[27][27][0) = _wmp[2]*_cou[11][27][1) +  fakac * _cou[11][11][1)
          _cou[24][28][0) = _wmp[0]*_cou[11][28][1) 
          _cou[21][28][0) = _wmp[1]*_cou[11][28][1)  & 
               +  faka3 * (_cou[8][28][0) - fakaac * _cou[8][28][1) ) + fakac * _cou[11][12][1)
          _cou[27][28][0) = _wmp[2]*_cou[11][28][1) +  fakac3 * _cou[11][18][1)
          _cou[24][29][0) = _wmp[0]*_cou[11][29][1) +  fakac2 * _cou[11][16][1)
          _cou[21][29][0) = _wmp[1]*_cou[11][29][1)  & 
               +  faka3 * (_cou[8][29][0) - fakaac * _cou[8][29][1) ) 
          _cou[27][29][0) = _wmp[2]*_cou[11][29][1) +  fakac2 * _cou[11][15][1)
          _cou[24][30][0) = _wmp[0]*_cou[11][30][1) 
          _cou[21][30][0) = _wmp[1]*_cou[11][30][1)  & 
               + faka3 * (_cou[8][30][0) - fakaac * _cou[8][30][1) ) + fakac2 * _cou[11][18][1)
          _cou[27][30][0) = _wmp[2]*_cou[11][30][1) +  fakac2 * _cou[11][17][1)
          _cou[24][31][0) = _wmp[0]*_cou[11][31][1) +  fakac2 * _cou[11][14][1)
          _cou[21][31][0) = _wmp[1]*_cou[11][31][1)  & 
               + faka3 * (_cou[8][31][0) - fakaac * _cou[8][31][1) ) + fakac2 * _cou[11][13][1)
          _cou[27][31][0) = _wmp[2]*_cou[11][31][1) 
          _cou[24][32][0) = _wmp[0]*_cou[11][32][1) +  fakac2 * _cou[11][19][1)
          _cou[21][32][0) = _wmp[1]*_cou[11][32][1)  & 
               +  faka3 * (_cou[8][32][0) - fakaac * _cou[8][32][1) ) + fakac * _cou[11][15][1)
          _cou[27][32][0) = _wmp[2]*_cou[11][32][1) +  fakac * _cou[11][13][1)
          _cou[24][33][0) = _wmp[0]*_cou[11][33][1) +  fakac * _cou[11][17][1)
          _cou[21][33][0) = _wmp[1]*_cou[11][33][1)  & 
               + faka3 * (_cou[8][33][0) - fakaac * _cou[8][33][1) ) + fakac2 * _cou[11][19][1)
          _cou[27][33][0) = _wmp[2]*_cou[11][33][1) +  fakac * _cou[11][14][1)
          _cou[24][34][0) = _wmp[0]*_cou[11][34][1) +  fakac * _cou[11][18][1)
          _cou[21][34][0) = _wmp[1]*_cou[11][34][1)  & 
               +  faka3 * (_cou[8][34][0) - fakaac * _cou[8][34][1) ) + fakac * _cou[11][16][1)
          _cou[27][34][0) = _wmp[2]*_cou[11][34][1) +  fakac2* _cou[11][19][1)
          _cou[26][20][0) = _wmp[0]*_cou[12][20][1) +  fakac4* _cou[12][10][1)
          _cou[28][20][0) = _wmp[1]*_cou[12][20][1) 
          _cou[22][20][0) = _wmp[2]*_cou[12][20][1)  & 
               +  faka3 * (_cou[9][20][0) - fakaac * _cou[9][20][1) ) 
          _cou[26][21][0) = _wmp[0]*_cou[12][21][1) 
          _cou[28][21][0) = _wmp[1]*_cou[12][21][1) +  fakac4 * _cou[12][11][1)
          _cou[22][21][0) = _wmp[2]*_cou[12][21][1)  & 
               +  faka3 * (_cou[9][21][0) - fakaac * _cou[9][21][1) ) 
          _cou[26][22][0) = _wmp[0]*_cou[12][22][1) 
          _cou[28][22][0) = _wmp[1]*_cou[12][22][1) 
          _cou[22][22][0) = _wmp[2]*_cou[12][22][1)  & 
               + faka3 * (_cou[9][22][0) - fakaac * _cou[9][22][1) ) + fakac4 * _cou[12][12][1)
          _cou[26][23][0) = _wmp[0]*_cou[12][23][1) +  fakac3 * _cou[12][13][1)
          _cou[28][23][0) = _wmp[1]*_cou[12][23][1) +  fakac * _cou[12][10][1)
          _cou[22][23][0) = _wmp[2]*_cou[12][23][1)  & 
               +  faka3 * (_cou[9][23][0) - fakaac * _cou[9][23][1) ) 
          _cou[26][24][0) = _wmp[0]*_cou[12][24][1) +  fakac * _cou[12][11][1)
          _cou[28][24][0) = _wmp[1]*_cou[12][24][1) +  fakac3 * _cou[12][14][1)
          _cou[22][24][0) = _wmp[2]*_cou[12][24][1)  & 
               +  faka3 * (_cou[9][24][0) - fakaac * _cou[9][24][1) ) 
          _cou[26][25][0) = _wmp[0]*_cou[12][25][1) +  fakac3 * _cou[12][15][1)
          _cou[28][25][0) = _wmp[1]*_cou[12][25][1) 
          _cou[22][25][0) = _wmp[2]*_cou[12][25][1)  & 
               +  faka3 * (_cou[9][25][0) - fakaac * _cou[9][25][1) ) + fakac * _cou[12][10][1)
          _cou[26][26][0) = _wmp[0]*_cou[12][26][1) +  fakac * _cou[12][12][1)
          _cou[28][26][0) = _wmp[1]*_cou[12][26][1) 
          _cou[22][26][0) = _wmp[2]*_cou[12][26][1)  & 
               + faka3 * (_cou[9][26][0) - fakaac * _cou[9][26][1) ) + fakac3 * _cou[12][16][1)
          _cou[26][27][0) = _wmp[0]*_cou[12][27][1) 
          _cou[28][27][0) = _wmp[1]*_cou[12][27][1) +  fakac3 * _cou[12][17][1)
          _cou[22][27][0) = _wmp[2]*_cou[12][27][1)  & 
               +  faka3 * (_cou[9][27][0) - fakaac * _cou[9][27][1) ) + fakac * _cou[12][11][1)
          _cou[26][28][0) = _wmp[0]*_cou[12][28][1) 
          _cou[28][28][0) = _wmp[1]*_cou[12][28][1) +  fakac * _cou[12][12][1)
          _cou[22][28][0) = _wmp[2]*_cou[12][28][1)  & 
               + faka3 * (_cou[9][28][0) - fakaac * _cou[9][28][1) ) + fakac3 * _cou[12][18][1)
          _cou[26][29][0) = _wmp[0]*_cou[12][29][1) +  fakac2 * _cou[12][16][1)
          _cou[28][29][0) = _wmp[1]*_cou[12][29][1) 
          _cou[22][29][0) = _wmp[2]*_cou[12][29][1)  & 
               + faka3 * (_cou[9][29][0) - fakaac * _cou[9][29][1) ) + fakac2 * _cou[12][15][1)
          _cou[26][30][0) = _wmp[0]*_cou[12][30][1) 
          _cou[28][30][0) = _wmp[1]*_cou[12][30][1) +  fakac2 * _cou[12][18][1)
          _cou[22][30][0) = _wmp[2]*_cou[12][30][1)  & 
               + faka3 * (_cou[9][30][0) - fakaac * _cou[9][30][1) ) + fakac2 * _cou[12][17][1)
          _cou[26][31][0) = _wmp[0]*_cou[12][31][1) +  fakac2 * _cou[12][14][1)
          _cou[28][31][0) = _wmp[1]*_cou[12][31][1) +  fakac2 * _cou[12][13][1)
          _cou[22][31][0) = _wmp[2]*_cou[12][31][1)  & 
               +  faka3 * (_cou[9][31][0) - fakaac * _cou[9][31][1) ) 
          _cou[26][32][0) = _wmp[0]*_cou[12][32][1) +  fakac2 * _cou[12][19][1)
          _cou[28][32][0) = _wmp[1]*_cou[12][32][1) +  fakac * _cou[12][15][1)
          _cou[22][32][0) = _wmp[2]*_cou[12][32][1)  & 
               +  faka3 * (_cou[9][32][0) - fakaac * _cou[9][32][1) ) + fakac * _cou[12][13][1)
          _cou[26][33][0) = _wmp[0]*_cou[12][33][1) +  fakac * _cou[12][17][1)
          _cou[28][33][0) = _wmp[1]*_cou[12][33][1) +  fakac2 * _cou[12][19][1)
          _cou[22][33][0) = _wmp[2]*_cou[12][33][1)  & 
               +  faka3 * (_cou[9][33][0) - fakaac * _cou[9][33][1) ) + fakac * _cou[12][14][1)
          _cou[26][34][0) = _wmp[0]*_cou[12][34][1) +  fakac * _cou[12][18][1)
          _cou[28][34][0) = _wmp[1]*_cou[12][34][1) +  fakac * _cou[12][16][1)
          _cou[22][34][0) = _wmp[2]*_cou[12][34][1)  & 
               + faka3 * (_cou[9][34][0) - fakaac * _cou[9][34][1) ) + fakac2 * _cou[12][19][1)
          _cou[31][20][0) = _wmp[1]*_cou[13][20][1)  & 
               +  faka * (_cou[7][20][0) - fakaac * _cou[7][20][1) ) 
          _cou[32][20][0) = _wmp[2]*_cou[13][20][1) 
          _cou[31][21][0) = _wmp[1]*_cou[13][21][1)  & 
               +  faka * (_cou[7][21][0) - fakaac * _cou[7][21][1) ) + fakac4 * _cou[13][11][1)
          _cou[32][21][0) = _wmp[2]*_cou[13][21][1) 
          _cou[31][22][0) = _wmp[1]*_cou[13][22][1)  & 
               +  faka * (_cou[7][22][0) - fakaac * _cou[7][22][1) ) 
          _cou[32][22][0) = _wmp[2]*_cou[13][22][1) +  fakac4 * _cou[13][12][1)
          _cou[31][23][0) = _wmp[1]*_cou[13][23][1)  & 
               +  faka * (_cou[7][23][0) - fakaac * _cou[7][23][1) ) + fakac * _cou[13][10][1)
          _cou[32][23][0) = _wmp[2]*_cou[13][23][1) 
          _cou[31][24][0) = _wmp[1]*_cou[13][24][1)  & 
               +  faka * (_cou[7][24][0) - fakaac * _cou[7][24][1) ) + fakac3 * _cou[13][14][1)
          _cou[32][24][0) = _wmp[2]*_cou[13][24][1) 
          _cou[31][25][0) = _wmp[1]*_cou[13][25][1)  & 
               +  faka * (_cou[7][25][0) - fakaac * _cou[7][25][1) ) 
          _cou[32][25][0) = _wmp[2]*_cou[13][25][1) +  fakac * _cou[13][10][1)
          _cou[31][26][0) = _wmp[1]*_cou[13][26][1)  & 
               +  faka * (_cou[7][26][0) - fakaac * _cou[7][26][1) ) 
          _cou[32][26][0) = _wmp[2]*_cou[13][26][1) +  fakac3* _cou[13][16][1)
          _cou[31][27][0) = _wmp[1]*_cou[13][27][1)  & 
               +  faka * (_cou[7][27][0) - fakaac * _cou[7][27][1) ) + fakac3* _cou[13][17][1)
          _cou[32][27][0) = _wmp[2]*_cou[13][27][1) +  fakac * _cou[13][11][1)
          _cou[31][28][0) = _wmp[1]*_cou[13][28][1)  & 
               +  faka * (_cou[7][28][0) - fakaac * _cou[7][28][1) ) + fakac * _cou[13][12][1)
          _cou[32][28][0) = _wmp[2]*_cou[13][28][1) +  fakac3 * _cou[13][18][1)
          _cou[31][29][0) = _wmp[1]*_cou[13][29][1)  & 
               +  faka * (_cou[7][29][0) - fakaac * _cou[7][29][1) ) 
          _cou[32][29][0) = _wmp[2]*_cou[13][29][1) +  fakac2 * _cou[13][15][1)
          _cou[31][30][0) = _wmp[1]*_cou[13][30][1)  & 
               +  faka * (_cou[7][30][0) - fakaac * _cou[7][30][1) ) + fakac2 * _cou[13][18][1)
          _cou[32][30][0) = _wmp[2]*_cou[13][30][1) +  fakac2 * _cou[13][17][1)
          _cou[31][31][0) = _wmp[1]*_cou[13][31][1)  & 
               +  faka * (_cou[7][31][0) - fakaac * _cou[7][31][1) ) + fakac2 * _cou[13][13][1)
          _cou[32][31][0) = _wmp[2]*_cou[13][31][1) 
          _cou[31][32][0) = _wmp[1]*_cou[13][32][1)  & 
               +  faka * (_cou[7][32][0) - fakaac * _cou[7][32][1) ) + fakac * _cou[13][15][1)
          _cou[32][32][0) = _wmp[2]*_cou[13][32][1) +  fakac * _cou[13][13][1)
          _cou[31][33][0) = _wmp[1]*_cou[13][33][1)  & 
               +  faka * (_cou[7][33][0) - fakaac * _cou[7][33][1) ) + fakac2 * _cou[13][19][1)
          _cou[32][33][0) = _wmp[2]*_cou[13][33][1) +  fakac * _cou[13][14][1)
          _cou[31][34][0) = _wmp[1]*_cou[13][34][1)  & 
               +  faka * (_cou[7][34][0) - fakaac * _cou[7][34][1) ) + fakac * _cou[13][16][1)
          _cou[32][34][0) = _wmp[2]*_cou[13][34][1) +  fakac2 * _cou[13][19][1)
          _cou[33][20][0) = _wmp[2]*_cou[14][20][1) 
          _cou[33][21][0) = _wmp[2]*_cou[14][21][1) 
          _cou[33][22][0) = _wmp[2]*_cou[14][22][1) +  fakac4 * _cou[14][12][1)
          _cou[33][23][0) = _wmp[2]*_cou[14][23][1) 
          _cou[33][24][0) = _wmp[2]*_cou[14][24][1) 
          _cou[33][25][0) = _wmp[2]*_cou[14][25][1) +  fakac * _cou[14][10][1)
          _cou[33][26][0) = _wmp[2]*_cou[14][26][1) +  fakac3 * _cou[14][16][1)
          _cou[33][27][0) = _wmp[2]*_cou[14][27][1) +  fakac * _cou[14][11][1)
          _cou[33][28][0) = _wmp[2]*_cou[14][28][1) +  fakac3 * _cou[14][18][1)
          _cou[33][29][0) = _wmp[2]*_cou[14][29][1) +  fakac2 * _cou[14][15][1)
          _cou[33][30][0) = _wmp[2]*_cou[14][30][1) +  fakac2 * _cou[14][17][1)
          _cou[33][31][0) = _wmp[2]*_cou[14][31][1) 
          _cou[33][32][0) = _wmp[2]*_cou[14][32][1) +  fakac * _cou[14][13][1)
          _cou[33][33][0) = _wmp[2]*_cou[14][33][1) +  fakac * _cou[14][14][1)
          _cou[33][34][0) = _wmp[2]*_cou[14][34][1) +  fakac2 * _cou[14][19][1)
          _cou[29][20][0) = _wmp[2]*_cou[15][20][1)  & 
               +  faka * (_cou[7][20][0) - fakaac * _cou[7][20][1) ) 
          _cou[29][21][0) = _wmp[2]*_cou[15][21][1)  & 
               +  faka * (_cou[7][21][0) - fakaac * _cou[7][21][1) ) 
          _cou[29][22][0) = _wmp[2]*_cou[15][22][1)  & 
               +  faka * (_cou[7][22][0) - fakaac * _cou[7][22][1) ) + fakac4 * _cou[15][12][1)
          _cou[29][23][0) = _wmp[2]*_cou[15][23][1)  & 
               +  faka * (_cou[7][23][0) - fakaac * _cou[7][23][1) ) 
          _cou[29][24][0) = _wmp[2]*_cou[15][24][1)  & 
               +  faka * (_cou[7][24][0) - fakaac * _cou[7][24][1) ) 
          _cou[29][25][0) = _wmp[2]*_cou[15][25][1)  & 
               +  faka * (_cou[7][25][0) - fakaac * _cou[7][25][1) ) + fakac * _cou[15][10][1)
          _cou[29][26][0) = _wmp[2]*_cou[15][26][1)  & 
               +  faka * (_cou[7][26][0) - fakaac * _cou[7][26][1) ) + fakac3 * _cou[15][16][1)
          _cou[29][27][0) = _wmp[2]*_cou[15][27][1)  & 
               +  faka * (_cou[7][27][0) - fakaac * _cou[7][27][1) ) + fakac * _cou[15][11][1)
          _cou[29][28][0) = _wmp[2]*_cou[15][28][1)  & 
               +  faka * (_cou[7][28][0) - fakaac * _cou[7][28][1) ) + fakac3 * _cou[15][18][1)
          _cou[29][29][0) = _wmp[2]*_cou[15][29][1)  & 
               +  faka * (_cou[7][29][0) - fakaac * _cou[7][29][1) ) + fakac2 * _cou[15][15][1)
          _cou[29][30][0) = _wmp[2]*_cou[15][30][1)  & 
               +  faka * (_cou[7][30][0) - fakaac * _cou[7][30][1) ) + fakac2 * _cou[15][17][1)
          _cou[29][31][0) = _wmp[2]*_cou[15][31][1)  & 
               +  faka * (_cou[7][31][0) - fakaac * _cou[7][31][1) ) 
          _cou[29][32][0) = _wmp[2]*_cou[15][32][1)  & 
               +  faka * (_cou[7][32][0) - fakaac * _cou[7][32][1) ) + fakac * _cou[15][13][1)
          _cou[29][33][0) = _wmp[2]*_cou[15][33][1)  & 
               +  faka * (_cou[7][33][0) - fakaac * _cou[7][33][1) ) + fakac * _cou[15][14][1)
          _cou[29][34][0) = _wmp[2]*_cou[15][34][1)  & 
               +  faka * (_cou[7][34][0) - fakaac * _cou[7][34][1) ) + fakac2 * _cou[15][19][1)
          _cou[34][20][0) = _wmp[1]*_cou[16][20][1) 
          _cou[34][21][0) = _wmp[1]*_cou[16][21][1) +  fakac4 * _cou[16][11][1)
          _cou[34][22][0) = _wmp[1]*_cou[16][22][1) 
          _cou[34][23][0) = _wmp[1]*_cou[16][23][1) +  fakac * _cou[16][10][1)
          _cou[34][24][0) = _wmp[1]*_cou[16][24][1) +  fakac3 * _cou[16][14][1)
          _cou[34][25][0) = _wmp[1]*_cou[16][25][1) 
          _cou[34][26][0) = _wmp[1]*_cou[16][26][1) 
          _cou[34][27][0) = _wmp[1]*_cou[16][27][1) +  fakac3 * _cou[16][17][1)
          _cou[34][28][0) = _wmp[1]*_cou[16][28][1) +  fakac * _cou[16][12][1)
          _cou[34][29][0) = _wmp[1]*_cou[16][29][1) 
          _cou[34][30][0) = _wmp[1]*_cou[16][30][1) +  fakac2 * _cou[16][18][1)
          _cou[34][31][0) = _wmp[1]*_cou[16][31][1) +  fakac2 * _cou[16][13][1)
          _cou[34][32][0) = _wmp[1]*_cou[16][32][1) +  fakac * _cou[16][15][1)
          _cou[34][33][0) = _wmp[1]*_cou[16][33][1) +  fakac2* _cou[16][19][1)
          _cou[34][34][0) = _wmp[1]*_cou[16][34][1) +  fakac * _cou[16][16][1)
          _cou[30][20][0) = _wmp[2]*_cou[17][20][1)  & 
               +  faka * (_cou[8][20][0) - fakaac * _cou[8][20][1) ) 
          _cou[30][21][0) = _wmp[2]*_cou[17][21][1)  & 
               +  faka * (_cou[8][21][0) - fakaac * _cou[8][21][1) ) 
          _cou[30][22][0) = _wmp[2]*_cou[17][22][1)  & 
               +  faka * (_cou[8][22][0) - fakaac * _cou[8][22][1) ) + fakac4 * _cou[17][12][1)
          _cou[30][23][0) = _wmp[2]*_cou[17][23][1)  & 
               +  faka * (_cou[8][23][0) - fakaac * _cou[8][23][1) ) 
          _cou[30][24][0) = _wmp[2]*_cou[17][24][1)  & 
               +  faka * (_cou[8][24][0) - fakaac * _cou[8][24][1) ) 
          _cou[30][25][0) = _wmp[2]*_cou[17][25][1)  & 
               +  faka * (_cou[8][25][0) - fakaac * _cou[8][25][1) ) + fakac * _cou[17][10][1)
          _cou[30][26][0) = _wmp[2]*_cou[17][26][1)  & 
               +  faka * (_cou[8][26][0) - fakaac * _cou[8][26][1) ) + fakac3 * _cou[17][16][1)
          _cou[30][27][0) = _wmp[2]*_cou[17][27][1)  & 
               +  faka * (_cou[8][27][0) - fakaac * _cou[8][27][1) ) + fakac * _cou[17][11][1)
          _cou[30][28][0) = _wmp[2]*_cou[17][28][1)  & 
               +  faka * (_cou[8][28][0) - fakaac * _cou[8][28][1) ) + fakac3 * _cou[17][18][1)
          _cou[30][29][0) = _wmp[2]*_cou[17][29][1)  & 
               +  faka * (_cou[8][29][0) - fakaac * _cou[8][29][1) ) + fakac2 * _cou[17][15][1)
          _cou[30][30][0) = _wmp[2]*_cou[17][30][1)  & 
               +  faka * (_cou[8][30][0) - fakaac * _cou[8][30][1) ) + fakac2 * _cou[17][17][1)
          _cou[30][31][0) = _wmp[2]*_cou[17][31][1)  & 
               +  faka * (_cou[8][31][0) - fakaac * _cou[8][31][1) ) 
          _cou[30][32][0) = _wmp[2]*_cou[17][32][1)  & 
               +  faka * (_cou[8][32][0) - fakaac * _cou[8][32][1) ) + fakac * _cou[17][13][1)
          _cou[30][33][0) = _wmp[2]*_cou[17][33][1)  & 
               +  faka * (_cou[8][33][0) - fakaac * _cou[8][33][1) ) + fakac * _cou[17][14][1)
          _cou[30][34][0) = _wmp[2]*_cou[17][34][1)  & 
               +  faka * (_cou[8][34][0) - fakaac * _cou[8][34][1) ) + fakac2 * _cou[17][19][1)

       endif

 
 
 
             */
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

