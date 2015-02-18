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
   ~Espfit(){};
    
   void setLog(Logger *log) { _log = log; }
    
   
    
   
   
  
    

   
   
   
   

    
    void Fit2Density(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_dftbasis, double _netcharge=0.0) {
        
    
   
        
 
    // setting up grid    
    Grid _grid;
    _grid.setupCHELPgrid(_atomlist);
    
    LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl;
        
    // Calculating nuclear potential at gridpoints
    
    ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());
    // ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_gridpoints.size());
    
    ub::vector<double> _NucPatGrid = EvalNuclearPotential(  _atomlist,  _grid );


    double Ztot = 0.0;
    for ( int j = 0; j < _atomlist.size(); j++){

           Ztot += _elements.getNucCrgECP(_atomlist[j]->type);  
    }
    
    ub::vector<double> DMATGSasarray=_dmat.data();
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;  
    #pragma omp parallel for
    for ( int i = 0 ; i < _grid.getsize(); i++){
        // AOESP matrix
         AOESP _aoesp;
         _aoesp.Initialize(_dftbasis._AOBasisSize);
         _aoesp.Fill(&_dftbasis, _grid.getGrid()[i]*NmtoBohr);
        ub::vector<double> AOESPasarray=_aoesp._aomatrix.data();
      
        for ( int _i =0; _i < DMATGSasarray.size(); _i++ ){
            _ESPatGrid(i) -= DMATGSasarray(_i)*AOESPasarray(_i);
        }   
    }
    
    
    _ESPatGrid += _NucPatGrid;

    std::vector< ub::vector<double> > _fitcenters;
    
          for ( int j = 0; j < _atomlist.size(); j++){

             ub::vector<double> _pos(3);
            _pos(0) = _atomlist[j]->x;
            _pos(1) = _atomlist[j]->y;
            _pos(2) = _atomlist[j]->z;
            _fitcenters.push_back(_pos);
            
          }
    
    
    
    
    std::vector<double> _charges = FitPartialCharges(_fitcenters,_grid, _ESPatGrid, _netcharge);
    
    //Write charges to qmatoms
        for ( int _i =0 ; _i < _atomlist.size(); _i++){
            _atomlist[_i]->charge=_charges[_i];
        }
    
    
    } 
    
private:
    
     Logger *_log;
     Elements _elements; 
     const double AtoBohr=1.8897259886;
     const double NmtoBohr=10*AtoBohr;
 
     
     
        
   ub::vector<double> EvalNuclearPotential( vector< QMAtom* >& _atoms, Grid _grid ){
    std::vector< ub::vector<double> >& _gridpoints =_grid.getGrid();   
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points"  << flush;
    ub::vector<double> _NucPatGrid(_gridpoints.size());
    
    for ( int i = 0 ; i < _gridpoints.size(); i++){
      double x_k = _gridpoints[i](0);
      double y_k = _gridpoints[i](1);
      double z_k = _gridpoints[i](2);


      for ( int j = 0; j < _atoms.size(); j++){

            double x_j = _atoms[j]->x;
            double y_j = _atoms[j]->y;
            double z_j = _atoms[j]->z;
	    double Znuc = _elements.getNucCrgECP(_atoms[j]->type);  

            double dist_j = sqrt( (x_j - x_k)*(x_j - x_k) +  (y_j - y_k)*(y_j - y_k) + (z_j - z_k)*(z_j - z_k)     )*NmtoBohr;

	    _NucPatGrid(i) += Znuc/dist_j;
        }
    }
       
       
       
    
    return _NucPatGrid;
       
       
   }
     
     
     
     
     
   void FitPartialCharges(Grid& _targetpotential,Grid& _chargepositions, double& netcharge){
    double AtoNm=0.1;
    double NmtoA=10.0;
    double HartreetoJoule;
    double JouletoHartree;

    if(_chargepositions.getsize() >_targetpotential.getsize()){
        throw std::runtime_error("Fit underdetermined, change grid options");
    }

    std::vector< APolarSite > _charges= _chargepositions.getSites();
    std::vector< APolarSite > _target= _targetpotential.getSites();

    std::vector< ub::vector<double> > _chargepos;

    std::vector< APolarSite >::iterator sit;
    for (sit=_charges.begin(); sit!=_charges.end(); ++sit) {
        ub::vector<double> temp= (sit->getPos()).converttoub();
        _chargepos.push_back(NmtoA*temp);    
    }
    
    
    ub::vector<double> _potential=ub::zero_vector<double>(_targetpotential.getsize());
    for( int i=0; i<_targetpotential.getsize();i++){
    _potential(i)=JouletoHartree*_target[i].getPhi();    
    }
    LOG(logDEBUG, *_log) << " Fitting APE to Chargeshell " << flush;
    std::vector<double>_chargesfromfit=FitPartialCharges(_chargepos,_targetpotential,_potential,netcharge);
    int state=0;
    for (int i=0; i<_charges.size();i++) {
        _charges[i].setQ00(_chargesfromfit[i],state);
    }   

       LOG(logDEBUG, *_log) << " Fitting completed " << flush;
   }
     
   
   \

     
     // Fits partial charges to Potential on a grid, constrains net charge
     std::vector<double> FitPartialCharges( std::vector< ub::vector<double> >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge ){
         
    std::vector< ub::vector<double> >& _gridpoints=_grid.getGrid();   
    // Fitting atomic partial charges
    ub::matrix<double> _Amat = ub::zero_matrix<double>(_fitcenters.size()+1,_fitcenters.size()+1);
    ub::matrix<double> _Bvec = ub::zero_matrix<double>(_fitcenters.size()+1);    
    
    // setting up _Amat
    for ( int _i =0 ; _i < _Amat.size1()-1; _i++){
        double x_i = _fitcenters[_i](0);
        double y_i = _fitcenters[_i](1);
        double z_i = _fitcenters[_i](2);
        
        for ( int _j=_i; _j<_Amat.size2()-1; _j++){
            double x_j = _fitcenters[_j](0);
            double y_j = _fitcenters[_j](1);
            double z_j = _fitcenters[_j](2);
            for ( int _k=0; _k < _gridpoints.size(); _k++){
            
                double x_k = _gridpoints[_k](0);
                double y_k = _gridpoints[_k](1);
                double z_k = _gridpoints[_k](2);
                
                double dist_i = sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     )*NmtoBohr;
                double dist_j = sqrt( (x_j - x_k)*(x_j - x_k) +  (y_j - y_k)*(y_j - y_k) + (z_j - z_k)*(z_j - z_k)     )*NmtoBohr;
                
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
    for ( int _i =0 ; _i < _Bvec.size1()-1; _i++){
        double x_i = _fitcenters[_i](0);
        double y_i = _fitcenters[_i](1);
        double z_i = _fitcenters[_i](2);
        for ( int _k=0; _k < _gridpoints.size(); _k++){
            
                double x_k = _gridpoints[_k](0);
                double y_k = _gridpoints[_k](1);
                double z_k = _gridpoints[_k](2);
                
                double dist_i = sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     )*NmtoBohr;
                _Bvec(_i,0) += _potential(_k)/dist_i;
                
        }
        
        
        
       }
    
    _Bvec(_Bvec.size1()-1,0) = _netcharge; //netcharge!!!!
    
    // invert _Amat
    ub::matrix<double> _Amat_inverse = ub::zero_matrix<double>(_fitcenters.size()+1,_fitcenters.size()+1);
    linalg_invert( _Amat , _Amat_inverse);
    //_Amat.resize(0,0);
    ub::matrix<double> _charges = ub::prod(_Amat_inverse,_Bvec);
    
   
    std::vector<double> _result;
    for ( int _i = 0; _i < _charges.size1(); _i++ ){
        
        _result.push_back(_charges(_i,0));
        
    }
       
    double _sumcrg = 0.0;
    for ( int _i =0 ; _i < _fitcenters.size(); _i++){
        
        LOG(logDEBUG, *_log) << " Center " << _i << " FitCharge: " << _result[_i] << flush;
        _sumcrg += _result[_i];
        
    }
    
    LOG(logDEBUG, *_log) << " Sum of fitted charges: " << _sumcrg << flush;
    
    // get RMSE
    double _rmse = 0.0;
    double _totalPotSq = 0.0;
    for ( int _k=0 ; _k < _gridpoints.size(); _k++ ){
        double x_k = _gridpoints[_k](0);
        double y_k = _gridpoints[_k](1);
        double z_k = _gridpoints[_k](2);
        double temp = 0.0;
        for ( int _i=0; _i < _fitcenters.size(); _i++ ){
            double x_i = _fitcenters[_i](0);
            double y_i = _fitcenters[_i](1);
            double z_i = _fitcenters[_i](2);
            
            double dist =  sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     )*NmtoBohr;
            temp += _result[_i]/dist;
        }
        _rmse += (_potential(_k) - temp)*(_potential(_k) - temp);
        _totalPotSq += _potential(_k)*_potential(_k);
        
    }
    LOG(logDEBUG, *_log) << " RMSE of fit:  " << sqrt(_rmse/_gridpoints.size()) << flush;
    LOG(logDEBUG, *_log) << " RRMSE of fit: " << sqrt(_rmse/_totalPotSq) << flush;
    
    return _result;
       
       
   } 
    
};    
    
    
    
    
    
    
    
}}

#endif	/* ESPFIT_H */
