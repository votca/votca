
#include <votca/ctp/espfit.h>
#include <votca/ctp/aomatrix.h>
#include <votca/tools/linalg.h>
#include <boost/progress.hpp>
#include <votca/ctp/numerical_integrations.h>
#include <math.h> 

using namespace std;
using namespace votca::tools;


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
void Espfit::EvaluateAPECharges(Grid& _targetgrid, Grid& _chargepositions){
    vector<APolarSite*> charges=_chargepositions.Sites();
    vector<APolarSite*> positions=_targetgrid.Sites();
    vector<APolarSite*>::iterator qit;
    vector<APolarSite*>::iterator pit;
    for (pit=positions.begin();pit!=positions.end();++pit){
        double potential=0.0;
        vec pos=(*pit)->getPos();
        for (qit=charges.begin();qit!=charges.end();++qit){
            double dist=abs((*qit)->getPos()-pos);
            potential+=((*qit)->getQ00())/dist;                              
            }
        (*pit)->setPhi(potential,0.0);
       }
    }

void Espfit::FitAPECharges(Grid& _targetgrid_fg, Grid& _targetgrid_bg, Grid& _chargepositions, double& netcharge){
    //double A2Bohr=1.8897259886;
     double Nm2Bohr=18.8972598860;
     //double Nm2A=10.0;
     //double A2nm=0.1;
    double Int2Hartree=Nm2Bohr;
    
    if(_chargepositions.getsize() >_targetgrid_fg.getsize()){
        throw std::runtime_error("Fit underdetermined, change grid options");
    }

    std::vector< APolarSite* > _charges= _chargepositions.Sites();
    std::vector< APolarSite* > _target_fg= _targetgrid_fg.Sites();
    std::vector< APolarSite* > _target_bg= _targetgrid_bg.Sites();
    LOG(logDEBUG, *_log) << " Grid of size fg:" << _targetgrid_fg.getsize() << flush;
    LOG(logDEBUG, *_log) << " Grid of size bg:" << _targetgrid_bg.getsize() << flush;
    LOG(logDEBUG, *_log) << " Chargepositions:" << _chargepositions.getsize() << flush;

    std::vector< ub::vector<double> > _chargepos;
    std::vector< APolarSite* >::iterator sit;
    for (sit=_charges.begin(); sit!=_charges.end(); ++sit) {
        ub::vector<double> temp= ((*sit)->getPos()).converttoub();
        _chargepos.push_back(temp);    
    }
    ub::vector<double> _potential=ub::zero_vector<double>(_targetgrid_fg.getsize());
    for( int i=0; i<_targetgrid_fg.getsize();i++){
    _potential(i)=Int2Hartree*(_target_fg[i]->getPhi()+_target_bg[i]->getPhi());    
    }

    LOG(logDEBUG, *_log) << " Fitting APE to Chargeshell with netcharge " << netcharge <<"e."<< flush;
    std::vector<double>_chargesfromfit=FitPartialCharges(_chargepos,_targetgrid_fg,_potential,netcharge);
    int state=0;
    for (unsigned i=0; i<_charges.size();i++) {
        _charges[i]->setQ00(_chargesfromfit[i],state);
    }   
       LOG(logDEBUG, *_log) << " Fitting completed " << flush;
   }
       


void Espfit::Fit2Density(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat, AOBasis &_basis,BasisSet &bs,string gridsize) { 
   
    double Nm2Bohr=18.8972598860;
    double A2nm=0.1;
    // setting up grid    
    Grid _grid;
    _grid.setAtomlist(&_atomlist);
    _grid.setupCHELPgrid();
    //_grid.printGridtoxyzfile("grid.xyz");
    LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl;
        
    // Calculating nuclear potential at gridpoints
    
    ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());
    // ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_gridpoints.size());
    
    AOOverlap overlap;
    overlap.Initialize(_basis._AOBasisSize);
    overlap.Fill(&_basis);
    ub::vector<double> DMATasarray=_dmat.data();
    ub::vector<double> AOOasarray=overlap._aomatrix.data();
    double N_comp=0.0;
    #pragma omp parallel for reduction(+:N_comp) 
    for ( unsigned _i =0; _i < DMATasarray.size(); _i++ ){
            N_comp =N_comp+ DMATasarray(_i)*AOOasarray(_i);
        } 
    

    NumericalIntegration numway;

    numway.GridSetup(gridsize,&bs,_atomlist);
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculate Densities at Numerical Grid with gridsize "<<gridsize  << flush; 
    double N=numway.IntegrateDensity_Atomblock(_dmat,&_basis);
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculated Densities at Numerical Grid, Number of electrons is "<< N << flush; 
    
    if(std::abs(N-N_comp)>0.001){
        LOG(logDEBUG, *_log) <<"=======================" << flush; 
        LOG(logDEBUG, *_log) <<"WARNING: Calculated Densities at Numerical Grid, Number of electrons "<< N <<" is far away from the the real value "<< N_comp<<", you should increase the accuracy of the integration grid."<< flush; 
        N=N_comp;
        LOG(logDEBUG, *_log) <<"WARNING: Electronnumber set to "<< N << flush; 
        LOG(logDEBUG, *_log) <<"=======================" << flush; 
    }
           
    double netcharge=getNetcharge( _atomlist,N );   
        
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush;     
    boost::progress_display show_progress( _grid.getsize() );
    #pragma omp parallel for
    for ( int i = 0 ; i < _grid.getsize(); i++){
        _ESPatGrid(i)=numway.IntegratePotential(_grid.getGrid()[i]*Nm2Bohr);
        ++show_progress;
    }
    
    LOG(logDEBUG, *_log) << TimeStamp() << " Electron contribution calculated"  << flush; 
    if (!_do_Transition){
    ub::vector<double> _NucPatGrid = EvalNuclearPotential(  _atomlist,  _grid );
    _ESPatGrid += _NucPatGrid;
    }
    
    std::vector< ub::vector<double> > _fitcenters;
    
    for ( unsigned j = 0; j < _atomlist.size(); j++){
       ub::vector<double> _pos(3);
      _pos(0) = A2nm*_atomlist[j]->x;
      _pos(1) = A2nm*_atomlist[j]->y;
      _pos(2) = A2nm*_atomlist[j]->z;
      _fitcenters.push_back(_pos);            
    }
  
    std::vector<double> _charges = FitPartialCharges(_fitcenters,_grid, _ESPatGrid, netcharge);
    
    //Write charges to qmatoms
    for ( unsigned _i =0 ; _i < _atomlist.size(); _i++){
        _atomlist[_i]->charge=_charges[_i];
    } 
    return;
    }


ub::vector<double> Espfit::EvalNuclearPotential(vector< QMAtom* >& _atoms, Grid _grid) {
    ub::vector<double> _NucPatGrid = ub::zero_vector<double>(_grid.getsize());
    double Nm2Bohr = 18.8972598860;
    double A2nm = 0.1;
    double Znuc=0.0;
    std::vector< ub::vector<double> >& _gridpoints = _grid.getGrid();
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points" << flush;
    
    for (unsigned i = 0; i < _gridpoints.size(); i++) {
        double x_k = _gridpoints[i](0);
        double y_k = _gridpoints[i](1);
        double z_k = _gridpoints[i](2);
        for (unsigned j = 0; j < _atoms.size(); j++) {

            double x_j = A2nm * _atoms[j]->x;
            double y_j = A2nm * _atoms[j]->y;
            double z_j = A2nm * _atoms[j]->z;
            if (_ECP) {
                Znuc = _elements.getNucCrgECP(_atoms[j]->type);
            } else {
                Znuc = _elements.getNucCrg(_atoms[j]->type);
            }
            double dist_j = sqrt((x_j - x_k)*(x_j - x_k) + (y_j - y_k)*(y_j - y_k) + (z_j - z_k)*(z_j - z_k)) * Nm2Bohr;
            _NucPatGrid(i) += Znuc / dist_j;
        }
    
    }
    return _NucPatGrid;
}

double Espfit::getNetcharge( vector< QMAtom* >& _atoms, double N ){
    double netcharge=0.0;
    if( std::abs(N)<0.05){
        LOG(logDEBUG, *_log) << "Number of Electrons is "<<N<< " transitiondensity is used for fit"  << flush;
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
            LOG(logDEBUG, *_log) <<"Number of Electrons minus ECP_Nucleus charge is "<<Znuc_ECP-N<< " you use ECPs, sounds good"  << flush;    
        }
        else if (std::abs(Znuc-N)<4){
            LOG(logDEBUG, *_log) <<"Number of Electrons minus real Nucleus charge is "<<Znuc-N<< " you are sure you want ECPs?"  << flush;    
        }
        else{
            LOG(logDEBUG, *_log) <<"Warning: Your molecule is highly ionized and you want ECPs, sounds interesting" << flush;    
        }
        netcharge=Znuc_ECP-N;          
    }
    else{
        if (std::abs(Znuc-N)<4){
            LOG(logDEBUG, *_log) <<"Number of Electrons minus Nucleus charge is "<<Znuc-N<< " you probably do not use ECPs, if you do use ECPs please use the option. Otherwise you are fine"  << flush;
        }
        else if(std::abs(Znuc_ECP-N)<4){
            LOG(logDEBUG, *_log) <<"Number of Electrons minus ECP_Nucleus charge is "<<Znuc_ECP-N<< " you probably use ECPs, if you do use ECPs please use the option to switch on"  << flush;
        }
        else{
            LOG(logDEBUG, *_log) <<"Warning: Your molecule is highly ionized and you use real core potentials, sounds interesting" << flush;    
        }
    }
    _do_Transition=false;
    }
    netcharge=round(netcharge);
    LOG(logDEBUG, *_log) <<"Netcharge constrained to " << netcharge<< flush; 
    
    return netcharge;
}    
    


void Espfit::Fit2Density_analytic(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat,AOBasis &_basis) { 
     double Nm2Bohr=18.8972598860;
     double A2nm=0.1;
    // setting up grid    
    Grid _grid;
    _grid.setAtomlist(&_atomlist);
    _grid.setupCHELPgrid();
    //_grid.printGridtoxyzfile("grid.xyz");
    LOG(logDEBUG, *_log) << TimeStamp() <<  " Done setting up CHELPG grid with " << _grid.getsize() << " points " << endl; 
    // Calculating nuclear potential at gridpoints
    
    ub::vector<double> _ESPatGrid = ub::zero_vector<double>(_grid.getsize());
    
    AOOverlap overlap;
    overlap.Initialize(_basis._AOBasisSize);
    overlap.Fill(&_basis);
    ub::vector<double> DMATasarray=_dmat.data();
    ub::vector<double> AOOasarray=overlap._aomatrix.data();
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
    
    LOG(logDEBUG, *_log) << TimeStamp() << " Calculating ESP at CHELPG grid points"  << flush; 
    #pragma omp parallel for
    for ( int i = 0 ; i < _grid.getsize(); i++){
        // AOESP matrix
         AOESP _aoesp;
         _aoesp.Initialize(_basis._AOBasisSize);
         _aoesp.Fill(&_basis, _grid.getGrid()[i]*Nm2Bohr);
        ub::vector<double> AOESPasarray=_aoesp._aomatrix.data();
      
        for ( unsigned _i =0; _i < DMATasarray.size(); _i++ ){
            _ESPatGrid(i) -= DMATasarray(_i)*AOESPasarray(_i);
        }   
    }
   
    

    

    std::vector< ub::vector<double> > _fitcenters;
    
          for ( unsigned j = 0; j < _atomlist.size(); j++){
             ub::vector<double> _pos(3);
            _pos(0) = A2nm*_atomlist[j]->x;
            _pos(1) = A2nm*_atomlist[j]->y;
            _pos(2) = A2nm*_atomlist[j]->z;
            _fitcenters.push_back(_pos);            
          } 
    std::vector<double> _charges = FitPartialCharges(_fitcenters,_grid, _ESPatGrid, netcharge);
    
    //Write charges to qmatoms
        for ( unsigned _i =0 ; _i < _atomlist.size(); _i++){
            _atomlist[_i]->charge=_charges[_i];
        }  
    return;
    } 

std::vector<double> Espfit::FitPartialCharges( std::vector< ub::vector<double> >& _fitcenters, Grid& _grid, ub::vector<double>& _potential, double& _netcharge ){
    LOG(logDEBUG, *_log) << TimeStamp() << " Setting up Matrices for fitting of size "<< _fitcenters.size()+1 <<" x " << _fitcenters.size()+1<< flush;    
    

    std::vector< ub::vector<double> >& _gridpoints=_grid.getGrid();   
    //cout << "x " << _gridpoints[0](0)<< " y " << _gridpoints[0](1)<< " z " << _gridpoints[0](1);
    //cout << "x " << _fitcenters[0](0)<< " y " << _fitcenters[0](1)<< " z " << _fitcenters[0](1);
    // Fitting atomic partial charges
      
    LOG(logDEBUG, *_log) << TimeStamp() << " Using "<< _fitcenters.size() <<" Fittingcenters and " << _gridpoints.size()<< " Gridpoints."<< flush;  
    
    ub::matrix<double> _Amat = ub::zero_matrix<double>(_fitcenters.size()+1,_fitcenters.size()+1);
    ub::matrix<double> _Bvec = ub::zero_matrix<double>(_fitcenters.size()+1,1);
    boost::progress_display show_progress( _fitcenters.size() );
    // setting up _Amat
    #pragma omp parallel for
    for ( unsigned _i =0 ; _i < _Amat.size1()-1; _i++){
        double x_i = _fitcenters[_i](0);
        double y_i = _fitcenters[_i](1);
        double z_i = _fitcenters[_i](2);
        ++show_progress;
        for ( unsigned _j=_i; _j<_Amat.size2()-1; _j++){
            double x_j = _fitcenters[_j](0);
            double y_j = _fitcenters[_j](1);
            double z_j = _fitcenters[_j](2);
            for ( unsigned _k=0; _k < _gridpoints.size(); _k++){
            
                double x_k = _gridpoints[_k](0);
                double y_k = _gridpoints[_k](1);
                double z_k = _gridpoints[_k](2);
                
                double dist_i = sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     );
                double dist_j = sqrt( (x_j - x_k)*(x_j - x_k) +  (y_j - y_k)*(y_j - y_k) + (z_j - z_k)*(z_j - z_k)     );
                
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
        double x_i = _fitcenters[_i](0);
        double y_i = _fitcenters[_i](1);
        double z_i = _fitcenters[_i](2);
        for ( unsigned _k=0; _k < _gridpoints.size(); _k++){
            
                double x_k = _gridpoints[_k](0);
                double y_k = _gridpoints[_k](1);
                double z_k = _gridpoints[_k](2);
                
                double dist_i = sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     );                
                _Bvec(_i,0) += _potential(_k)/dist_i;                
        }
       }
    
    _Bvec(_Bvec.size1()-1,0) = _netcharge; //netcharge!!!!
    LOG(logDEBUG, *_log) << TimeStamp() << "  Inverting Matrices "<< flush;      
    // invert _Amat
    ub::matrix<double> _Amat_inverse = ub::zero_matrix<double>(_fitcenters.size()+1,_fitcenters.size()+1);
    
    
    
    if(_do_svd){
        int notfittedatoms=linalg_invert_svd( _Amat , _Amat_inverse,_conditionnumber);
        LOG(logDEBUG, *_log) << TimeStamp() << "SVD Done. "<<notfittedatoms<<" could not be fitted and are set to zero."<< flush;
    }
    else{
        linalg_invert( _Amat , _Amat_inverse);
    }

    LOG(logDEBUG, *_log) << TimeStamp() << " Inverting Matrices done."<< flush;    
    //_Amat.resize(0,0);

    
    
    ub::matrix<double> _charges = ub::prod(_Amat_inverse,_Bvec);
   
    std::vector<double> _result;
    for ( unsigned _i = 0; _i < _charges.size1(); _i++ ){       
        _result.push_back(_charges(_i,0));        
    }
       
    double _sumcrg = 0.0;
    for ( unsigned _i =0 ; _i < _fitcenters.size(); _i++){
        
        LOG(logDEBUG, *_log) << " Center " << _i << " FitCharge: " << _result[_i] << flush;
        _sumcrg += _result[_i];       
    }
    
    LOG(logDEBUG, *_log) << " Sum of fitted charges: " << _sumcrg << flush;
    
    // get RMSE
    double _rmse = 0.0;
    double _totalPotSq = 0.0;
    for ( unsigned _k=0 ; _k < _gridpoints.size(); _k++ ){
        double x_k = _gridpoints[_k](0);
        double y_k = _gridpoints[_k](1);
        double z_k = _gridpoints[_k](2);
        double temp = 0.0;
        for ( unsigned _i=0; _i < _fitcenters.size(); _i++ ){
            double x_i = _fitcenters[_i](0);
            double y_i = _fitcenters[_i](1);
            double z_i = _fitcenters[_i](2);
            
            double dist =  sqrt( (x_i - x_k)*(x_i - x_k) +  (y_i - y_k)*(y_i - y_k) + (z_i - z_k)*(z_i - z_k)     );
            temp += _result[_i]/dist;
        }
        _rmse += (_potential(_k) - temp)*(_potential(_k) - temp);
        _totalPotSq += _potential(_k)*_potential(_k);        
    }
    LOG(logDEBUG, *_log) << " RMSE of fit:  " << sqrt(_rmse/_gridpoints.size()) << flush;
    LOG(logDEBUG, *_log) << " RRMSE of fit: " << sqrt(_rmse/_totalPotSq) << flush;
    
    return _result;     
   }     


}}