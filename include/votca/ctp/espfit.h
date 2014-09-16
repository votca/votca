/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __CTP_ESPFIT__H
#define	__CTP_ESPFIT__H


#include <votca/ctp/elements.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include <votca/ctp/grid.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/qmatom.h>


/**
* \brief Takes a list of atoms, and the corresponding density matrix and puts out a table of partial charges
*
* 
* 
*/
using namespace std;
using namespace votca::tools;


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
class Espfit{
public:
    
    Espfit(){};
    //Espfit(vector< QMAtom* >& Atomlist, ub::matrix<double> &DMAT, AOBasis &dftbasis)
     //   :_atomlist(Atomlist), _dmat(DMAT), _dftbasis(dftbasis)
       // {};
   ~Espfit(){};
    
   void setLog(Logger *log) { _log = log; }
    
    
    
    

    
    void FittoDensity(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_dftbasis) {
        
    
   
        
    Elements _elements;    
    // setting up grid    
    Grid _grid;
    _grid.setupCHELPgrid(_atomlist);
    std::vector< ub::vector<double> > &_gridpoints=_grid.getGrid();
    LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << _gridpoints.size() << " points " << endl;
    _grid.printGridtofile("grid.xyz");
    
    // Calculating nuclear charge at gridpoints
    
    ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_gridpoints.size());
    ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_gridpoints.size());
    
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points"  << flush;
    
    for ( int i = 0 ; i < _gridpoints.size(); i++){
      double x_k = _gridpoints[i](0);
      double y_k = _gridpoints[i](1);
      double z_k = _gridpoints[i](2);


      for ( int j = 0; j < _atomlist.size(); j++){

            double x_j = _atomlist[j]->x;
            double y_j = _atomlist[j]->y;
            double z_j = _atomlist[j]->z;
	    double Znuc = _elements.getNucCrgECP(_atomlist[j]->type);  

            double dist_j = sqrt( (x_j - x_k)*(x_j - x_k) +  (y_j - y_k)*(y_j - y_k) + (z_j - z_k)*(z_j - z_k)     )*1.8897259886;

	    _NucPatGrid(i) += Znuc/dist_j;
        }
    }

    double Ztot = 0.0;
    for ( int j = 0; j < _atomlist.size(); j++){

           Ztot += _elements.getNucCrgECP(_atomlist[j]->type);  
    }
    
    ub::vector<double> DMATGSasarray=_dmat.data();
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;  
    #pragma omp parallel for
    for ( int i = 0 ; i < _gridpoints.size(); i++){
        // AOESP matrix
         AOESP _aoesp;
         _aoesp.Initialize(_dftbasis._AOBasisSize);
         _aoesp.Fill(&_dftbasis, _gridpoints[i]*1.8897259886);
        ub::vector<double> AOESPasarray=_aoesp._aomatrix.data();
      
        for ( int _i =0; _i < DMATGSasarray.size(); _i++ ){
            _ESPatGrid(i) -= DMATGSasarray(_i)*AOESPasarray(_i);
        }   
    }
    
    
    _ESPatGrid += _NucPatGrid;

    // Fitting atomic partial charges
    ub::matrix<double> _Amat = ub::zero_matrix<double>(_atomlist.size()+1,_atomlist.size()+1);
    ub::matrix<double> _Bvec = ub::zero_matrix<double>(_atomlist.size()+1);    
    
    // setting up _Amat
    LOG(logDEBUG, *_log) << TimeStamp() << "  Setting up Amat of size "  << _Amat.size1() << " by " << _Amat.size2() << flush; 
    for ( int _i =0 ; _i < _Amat.size1()-1; _i++){
        double x_i = _atomlist[_i]->x;
        double y_i = _atomlist[_i]->y;
        double z_i = _atomlist[_i]->z;
        
        for ( int _j=_i; _j<_Amat.size2()-1; _j++){
            double x_j = _atomlist[_j]->x;
            double y_j = _atomlist[_j]->y;
            double z_j = _atomlist[_j]->z;
            for ( int _k=0; _k < _gridpoints.size(); _k++){
            
                double x_k = _gridpoints[_k](0);
                double y_k = _gridpoints[_k](1);
                double z_k = _gridpoints[_k](2);
                
                double dist_i = sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     )*1.8897259886;
                double dist_j = sqrt( (x_j - x_k)*(x_j - x_k) +  (y_j - y_k)*(y_j - y_k) + (z_j - z_k)*(z_j - z_k)     )*1.8897259886;
                
                 _Amat(_i,_j) += 1.0/dist_i/dist_j; 
                
            }
            _Amat(_j,_i) = _Amat(_i,_j);
        }
        
    }
    
    for ( int _i =0 ; _i < _Amat.size1(); _i++){
      _Amat(_i,_Amat.size1()-1) = 1.0;
      _Amat(_Amat.size1()-1,_i) = 1.0;
    }
    _Amat(_Amat.size1()-1,_Amat.size1()-1) = 0.0;

    

    // setting up Bvec
      LOG(logDEBUG, *_log) << TimeStamp() << "  Setting up Bvec"  << flush; 
    for ( int _i =0 ; _i < _Bvec.size1()-1; _i++){
        double x_i = _atomlist[_i]->x;
        double y_i = _atomlist[_i]->y;
        double z_i = _atomlist[_i]->z;
        for ( int _k=0; _k < _gridpoints.size(); _k++){
            
                double x_k = _gridpoints[_k](0);
                double y_k = _gridpoints[_k](1);
                double z_k = _gridpoints[_k](2);
                
                double dist_i = sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     )*1.8897259886;
                _Bvec(_i,0) += _ESPatGrid(_k)/dist_i;
                
        }
        
        
        
       }
    
    _Bvec(_Bvec.size1()-1,0) = 0.0; //netcharge!!!!
    
    // invert _Amat
     LOG(logDEBUG, *_log) << TimeStamp() << "  Inverting _Amat"  << flush; 
    ub::matrix<double> _Amat_inverse = ub::zero_matrix<double>(_atomlist.size()+1,_atomlist.size()+1);
    linalg_invert( _Amat , _Amat_inverse);
    //_Amat.resize(0,0);
    LOG(logDEBUG, *_log) << TimeStamp() << "  Calculating CHELPG charges"  << flush; 
    ub::matrix<double> _charges = ub::prod(_Amat_inverse,_Bvec);
    
    double _sumcrg = 0.0;
    for ( int _i =0 ; _i < _Bvec.size1()-1; _i++){
        
        LOG(logDEBUG, *_log) << " Atom " << _i << " Type " << _atomlist[_i]->type << " Charge: " << _charges(_i,0) << flush;
        _sumcrg += _charges(_i,0);
        
    }
    
    LOG(logDEBUG, *_log) << " Sum of partial charges: " << _sumcrg << flush;
    
    // get RMSE
    double _rmse = 0.0;
    double _totalESPsq = 0.0;
    for ( int _k=0 ; _k < _gridpoints.size(); _k++ ){
        double x_k = _gridpoints[_k](0);
        double y_k = _gridpoints[_k](1);
        double z_k = _gridpoints[_k](2);
        double temp = 0.0;
        for ( int _i=0; _i < _atomlist.size(); _i++ ){
            double x_i = _atomlist[_i]->x;
            double y_i = _atomlist[_i]->y;
            double z_i = _atomlist[_i]->z;
            
            double dist =  sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     )*1.8897259886;
            temp += _charges(_i,0)/dist;
        }
        _rmse += (_ESPatGrid(_k) - temp)*(_ESPatGrid(_k) - temp);
        _totalESPsq += _ESPatGrid(_k)*_ESPatGrid(_k);
        
    }
    LOG(logDEBUG, *_log) << " RMSE of fit:  " << sqrt(_rmse/_gridpoints.size()) << flush;
    LOG(logDEBUG, *_log) << " RRMSE of fit: " << sqrt(_rmse/_totalESPsq) << flush;
    
    
    
    
    
    } 
    
private:
    
     Logger *_log;

 
    
    
};    
    
    
    
    
    
    
    
}}

#endif	/* ESPFIT_H */
