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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/numerical_integrations.h>
#include <votca/ctp/logger.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/aomatrix.h>
#include <votca/tools/linalg.h>
//#include <boost/progress.hpp>

#include <math.h>
#include <votca/tools/constants.h>

using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;





void Espfit::Fit2Density(std::vector< ctp::QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,BasisSet &bs,string gridsize) {


    // setting up grid
    Grid _grid;
    _grid.setAtomlist(&_atomlist);
    _grid.setupCHELPgrid();
    //_grid.printGridtoxyzfile("grid.xyz");
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl;

    // Calculating nuclear potential at gridpoints

    ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());

    AOOverlap overlap;
    overlap.Initialize(_basis.AOBasisSize());
    overlap.Fill(_basis);
    ub::vector<double> DMATasarray=_dmat.data();
    ub::vector<double> AOOasarray=overlap.Matrix().data();
    double N_comp=0.0;
    #pragma omp parallel for reduction(+:N_comp)
    for ( unsigned _i =0; _i < DMATasarray.size(); _i++ ){
            N_comp =N_comp+ DMATasarray(_i)*AOOasarray(_i);
        }

    NumericalIntegration numway;

    numway.GridSetup(gridsize,&bs,_atomlist,&_basis);
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculate Densities at Numerical Grid with gridsize "<<gridsize  << flush;
    double N=numway.IntegrateDensity(_dmat);
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculated Densities at Numerical Grid, Number of electrons is "<< N << flush;

    if(std::abs(N-N_comp)>0.001){
        CTP_LOG(ctp::logDEBUG, *_log) <<"=======================" << flush;
        CTP_LOG(ctp::logDEBUG, *_log) <<"WARNING: Calculated Densities at Numerical Grid, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid."<< flush;
        N=N_comp;
        CTP_LOG(ctp::logDEBUG, *_log) <<"WARNING: Electronnumber set to "<< N << flush;
        CTP_LOG(ctp::logDEBUG, *_log) <<"=======================" << flush;
    }

    double netcharge=getNetcharge( _atomlist,N );

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;
    //boost::progress_display show_progress( _grid.getsize() );
    #pragma omp parallel for
    for ( int i = 0 ; i < _grid.getsize(); i++){
        _ESPatGrid(i)=numway.IntegratePotential(_grid.getGrid()[i]*tools::conv::nm2bohr);
        //++show_progress;
    }

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Electron contribution calculated"  << flush;
    if (!_do_Transition){
    ub::vector<double> _NucPatGrid = EvalNuclearPotential(  _atomlist,  _grid );
    _ESPatGrid += _NucPatGrid;
    }

    std::vector< tools::vec > _fitcenters;

    for ( unsigned j = 0; j < _atomlist.size(); j++){
       tools::vec _pos=_atomlist[j]->getPos()*tools::conv::ang2nm;
      _fitcenters.push_back(_pos);
    }

    std::vector<double> _charges = FitPartialCharges(_fitcenters,_grid, _ESPatGrid, netcharge);

    //Write charges to qmatoms
    for ( unsigned _i =0 ; _i < _atomlist.size(); _i++){
        _atomlist[_i]->charge=_charges[_i];
    }
    return;
    }


ub::vector<double> Espfit::EvalNuclearPotential(std::vector< ctp::QMAtom* >& _atoms, Grid _grid) {
    ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_grid.getsize());

    double Znuc=0.0;
    const std::vector< vec >& _gridpoints = _grid.getGrid();
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points" << flush;

    for (unsigned i = 0; i < _gridpoints.size(); i++) {
        for (unsigned j = 0; j < _atoms.size(); j++) {
            vec posatom=_atoms[j]->getPos()*tools::conv::ang2nm;
            if (_ECP) {
                Znuc = _elements.getNucCrgECP(_atoms[j]->type);
            } else {
                Znuc = _elements.getNucCrg(_atoms[j]->type);
            }
            double dist_j = tools::abs(_gridpoints[i]-posatom) * tools::conv::nm2bohr;
            _NucPatGrid(i) += Znuc / dist_j;
        }

    }
    return _NucPatGrid;
}

double Espfit::getNetcharge( std::vector< ctp::QMAtom* >& _atoms, double N ){
    double netcharge=0.0;
    if( std::abs(N)<0.05){
        //CTP_LOG(ctp::logDEBUG, *_log) << "Number of Electrons is "<<N<< " transitiondensity is used for fit"  << flush;
        _do_Transition=true;
    }
    else{
    double Znuc_ECP = 0.0;
    double Znuc=0.0;
    for ( unsigned j = 0; j < _atoms.size(); j++){
           Znuc_ECP += _elements.getNucCrgECP(_atoms[j]->type);
           Znuc+= _elements.getNucCrg(_atoms[j]->type);
    }

    if (_ECP){
        if (std::abs(Znuc_ECP-N)<4){
            CTP_LOG(ctp::logDEBUG, *_log) <<"Number of Electrons minus ECP_Nucleus charge is "<<Znuc_ECP-N<< " you use ECPs, sounds good"  << flush;
        }
        else if (std::abs(Znuc-N)<4){
            CTP_LOG(ctp::logDEBUG, *_log) <<"Number of Electrons minus real Nucleus charge is "<<Znuc-N<< " you are sure you want ECPs?"  << flush;
        }
        else{
            CTP_LOG(ctp::logDEBUG, *_log) <<"Warning: Your molecule is highly ionized and you want ECPs, sounds interesting" << flush;
        }
        netcharge=Znuc_ECP-N;
    }
    else{
        if (std::abs(Znuc-N)<4){
            CTP_LOG(ctp::logDEBUG, *_log) <<"Number of Electrons minus Nucleus charge is "<<Znuc-N<< " you probably do not use ECPs, if you do use ECPs please use the option. Otherwise you are fine"  << flush;
        }
        else if(std::abs(Znuc_ECP-N)<4){
            CTP_LOG(ctp::logDEBUG, *_log) <<"Number of Electrons minus ECP_Nucleus charge is "<<Znuc_ECP-N<< " you probably use ECPs, if you do use ECPs please use the option to switch on"  << flush;
        }
        else{
            CTP_LOG(ctp::logDEBUG, *_log) <<"Warning: Your molecule is highly ionized and you use real core potentials, sounds interesting" << flush;
        }
    }
    _do_Transition=false;
    }
    netcharge=round(netcharge);
    CTP_LOG(ctp::logDEBUG, *_log) <<"Netcharge constrained to " << netcharge<< flush;

    return netcharge;
}



void Espfit::Fit2Density_analytic(std::vector< ctp::QMAtom* >& _atomlist, ub::matrix<double> &_dmat,AOBasis &_basis) {
     double Nm2Bohr=tools::conv::nm2bohr;
     double A2nm=tools::conv::ang2nm;
    // setting up grid
    Grid _grid;
    _grid.setAtomlist(&_atomlist);
    _grid.setupCHELPgrid();

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl;
    // Calculating nuclear potential at gridpoints

    ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());

    AOOverlap overlap;
    overlap.Initialize(_basis.AOBasisSize());
    overlap.Fill(_basis);
    const ub::vector<double> DMATasarray=_dmat.data();
    const ub::vector<double> AOOasarray=overlap.Matrix().data();
    double N=0.0;
    #pragma omp parallel for reduction(+:N)
    for ( unsigned _i =0; _i < DMATasarray.size(); _i++ ){
            N =N+ DMATasarray(_i)*AOOasarray(_i);
        }

    double netcharge=getNetcharge( _atomlist,N );
    if(!_do_Transition){
        ub::vector<double> _NucPatGrid = EvalNuclearPotential(  _atomlist,  _grid);
        _ESPatGrid += _NucPatGrid;
    }

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;
    #pragma omp parallel for
    for ( int i = 0 ; i < _grid.getsize(); i++){
        // AOESP matrix
         AOESP _aoesp;
         _aoesp.Initialize(_basis.AOBasisSize());
         _aoesp.Fill(_basis, _grid.getGrid()[i]*Nm2Bohr);
        const ub::vector<double> AOESPasarray=_aoesp.Matrix().data();

        for ( unsigned _i =0; _i < DMATasarray.size(); _i++ ){
            _ESPatGrid(i) -= DMATasarray(_i)*AOESPasarray(_i);
        }
    }

    std::vector< tools::vec > _fitcenters;

          for ( unsigned j = 0; j < _atomlist.size(); j++){
             tools::vec _pos=A2nm*_atomlist[j]->getPos();
            _fitcenters.push_back(_pos);
          }
    std::vector<double> _charges = FitPartialCharges(_fitcenters,_grid, _ESPatGrid, netcharge);

    //Write charges to qmatoms
        for ( unsigned _i =0 ; _i < _atomlist.size(); _i++){
            _atomlist[_i]->charge=_charges[_i];
        }
    return;
    }

std::vector<double> Espfit::FitPartialCharges( std::vector< tools::vec >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge ){
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Setting up Matrices for fitting of size "<< _fitcenters.size()+1 <<" x " << _fitcenters.size()+1<< flush;

    const std::vector< tools::vec >& _gridpoints=_grid.getGrid();

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Using "<< _fitcenters.size() <<" Fittingcenters and " << _gridpoints.size()<< " Gridpoints."<< flush;

    ub::matrix<double> _Amat = ub::zero_matrix<double>(_fitcenters.size()+1,_fitcenters.size()+1);
    ub::matrix<double> _Bvec = ub::zero_matrix<double>(_fitcenters.size()+1,1);
    //boost::progress_display show_progress( _fitcenters.size() );
    // setting up _Amat
    #pragma omp parallel for
    for ( unsigned _i =0 ; _i < _Amat.size1()-1; _i++){
        for ( unsigned _j=_i; _j<_Amat.size2()-1; _j++){
            for ( unsigned _k=0; _k < _gridpoints.size(); _k++){
                double dist_i = tools::abs(_fitcenters[_i]-_gridpoints[_k])*tools::conv::nm2bohr;
                double dist_j = tools::abs(_fitcenters[_j]-_gridpoints[_k])*tools::conv::nm2bohr;

                 _Amat(_i,_j) += 1.0/dist_i/dist_j;
            }
            _Amat(_j,_i) = _Amat(_i,_j);
        }
    }

    for ( unsigned _i =0 ; _i < _Amat.size1(); _i++){
      _Amat(_i,_Amat.size1()-1) = 1.0;
      _Amat(_Amat.size1()-1,_i) = 1.0;
    }
    _Amat(_Amat.size1()-1,_Amat.size1()-1) = 0.0;

    // setting up Bvec
    #pragma omp parallel for
    for ( unsigned _i =0 ; _i < _Bvec.size1()-1; _i++){
        for ( unsigned _k=0; _k < _gridpoints.size(); _k++){
                double dist_i = tools::abs(_fitcenters[_i]-_gridpoints[_k])*tools::conv::nm2bohr;
                _Bvec(_i,0) += _potential(_k)/dist_i;
        }
       }

    _Bvec(_Bvec.size1()-1,0) = _netcharge; //netcharge!!!!
    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << "  Inverting Matrices "<< flush;
    // invert _Amat
    ub::matrix<double> _Amat_inverse = ub::zero_matrix<double>(_fitcenters.size()+1,_fitcenters.size()+1);



    if(_do_svd){
        int notfittedatoms=linalg_invert_svd( _Amat , _Amat_inverse,_conditionnumber);
        CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << "SVD Done. "<<notfittedatoms<<" could not be fitted and are set to zero."<< flush;
    }
    else{
        linalg_invert( _Amat , _Amat_inverse);
    }

    CTP_LOG(ctp::logDEBUG, *_log) << ctp::TimeStamp() << " Inverting Matrices done."<< flush;
    //_Amat.resize(0,0);



    ub::matrix<double> _charges = ub::prod(_Amat_inverse,_Bvec);

    std::vector<double> _result;
    for ( unsigned _i = 0; _i < _charges.size1(); _i++ ){
        _result.push_back(_charges(_i,0));
    }

    double _sumcrg = 0.0;
    for ( unsigned _i =0 ; _i < _fitcenters.size(); _i++){

        //CTP_LOG(ctp::logDEBUG, *_log) << " Center " << _i << " FitCharge: " << _result[_i] << "pos " << _fitcenters[_i] << flush;
        _sumcrg += _result[_i];
    }

    CTP_LOG(ctp::logDEBUG, *_log) << " Sum of fitted charges: " << _sumcrg << flush;

    // get RMSE
    double _rmse = 0.0;
    double _totalPotSq = 0.0;
    for ( unsigned _k=0 ; _k < _gridpoints.size(); _k++ ){
        double temp = 0.0;
        for ( unsigned _i=0; _i < _fitcenters.size(); _i++ ){
            double dist =  tools::abs(_gridpoints[_k]-_fitcenters[_i])*tools::conv::nm2bohr;
            temp += _result[_i]/dist;
        }
        _rmse += (_potential(_k) - temp)*(_potential(_k) - temp);
        _totalPotSq += _potential(_k)*_potential(_k);
    }
    CTP_LOG(ctp::logDEBUG, *_log) << " RMSE of fit:  " << sqrt(_rmse/_gridpoints.size()) << flush;
    CTP_LOG(ctp::logDEBUG, *_log) << " RRMSE of fit: " << sqrt(_rmse/_totalPotSq) << flush;

    return _result;
   }


}}
