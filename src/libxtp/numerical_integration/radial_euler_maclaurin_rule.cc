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
#include <votca/xtp/votca_xtp_config.h>

#include <boost/math/constants/constants.hpp>
#include "votca/xtp/radial_euler_maclaurin_rule.h"
#include "votca/xtp/aobasis.h"
#include <votca/xtp/aomatrix.h>

using namespace std;


namespace votca { namespace xtp { 
    
    
    std::vector<double> EulerMaclaurinGrid::getPruningIntervals(string element){
        
        std::vector<double> _r;
        
        // get Bragg-Slater Radius for this element
        double BSradius = _BraggSlaterRadii.at(element);
        
        // row type of element
        int RowType = _period_row.at(element);
        
        if ( RowType == 1 ){
            
            
            _r.push_back( 0.25 * BSradius );
            _r.push_back( 0.5 * BSradius );
            _r.push_back( 1.0 * BSradius );
            _r.push_back( 4.5 * BSradius );
            
        } else if ( RowType == 2 ){
            
            _r.push_back( 0.1667 * BSradius );
            _r.push_back( 0.5 * BSradius );
            _r.push_back( 0.9 * BSradius );
            _r.push_back( 3.5 * BSradius );            
            
            
        } else if ( RowType == 3 ) {
            
            _r.push_back( 0.1 * BSradius );
            _r.push_back( 0.4 * BSradius );
            _r.push_back( 0.8 * BSradius );
            _r.push_back( 2.5 * BSradius );
            
        } else {
            
            cerr << "Pruning unsupported for RowType " << RowType << endl;
            exit(1);
        }

        return _r;
        
        
        /*
         
         
               Data asg1/0.25  , 0.5, 1.0, 4.5,
     &          0.1667, 0.5, 0.9, 3.5,
     &          0.1   , 0.4, 0.8, 2.5/
         
         */
        
        
        
        
        
    }
    
    
    
    
    
    
    
void EulerMaclaurinGrid::getRadialCutoffs(vector<QMAtom* > _atoms, BasisSet* bs, string gridtype) {

            map<string, min_exp>::iterator it;
            vector< QMAtom* > ::iterator ait;
            vector< QMAtom* > ::iterator bit;
            
            double eps = Accuracy[gridtype];
            double _decaymin;
            int _lvalue;
          //  cout << endl << " Setting cutoffs for grid type " << gridtype << " eps = " << eps << endl;
            // 1) is only element based
            // loop over atoms
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get element type of the atom
                string name = (*ait)->type;
                // is this element already in map?
                it = _element_ranges.find(name);
                // only proceed, if element data does not exist yet
                if (it == _element_ranges.end()) {
                    // get first range estimate and add to map
                    min_exp this_atom;
                    double range_max = 0.0;
                    // get the basis set entry for this element
                    Element* _element = bs->getElement(name);
                    // and loop over all shells to figure out minimum decay constant and angular momentum of this function
                    for (Element::ShellIterator its = _element->firstShell(); its != _element->lastShell(); its++) {
                        Shell* shell = (*its);
                        _decaymin = 1e7;
                        _lvalue = 0;
                        int _lmax = shell->getLmax();
                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                            GaussianPrimitive* gaussian = *itg;
                            double _decay = gaussian->decay;
                            if (_decay < _decaymin) {
                                _decaymin = _decay;
                                _lvalue = _lmax;
                            }
                        }
                        double range = DetermineCutoff( 2*_decaymin, 2*_lvalue+2, eps);
                        if ( range > range_max ){
                            this_atom.alpha = _decaymin;
                            this_atom.l     = _lvalue;
                            this_atom.range = range;
                            range_max = range;
                        }
                    } // shells

                 //   cout << "Element " << name << " alpha " << this_atom.alpha << " l " << this_atom.l << " Rcut " << this_atom.range << endl;
                    _element_ranges[name] = this_atom;
                } // new element
            } // atoms

         //   exit(0);

            // calculate overlap matrix
            AOBasis aobasis;
            aobasis.AOBasisFill(bs, _atoms);
            AOOverlap _overlap;
            // initialize overlap matrix
            _overlap.Initialize(aobasis._AOBasisSize);
            // Fill overlap
            _overlap.Fill(&aobasis);
            
            
            
            // refining by going through all atom combinations
            // get collapsed index list
            int atidx = 0;
            std::vector<int> idxstart;
            std::vector<int> idxstop;
            int start = 0;
            int end = 0;
            for (vector< AOShell* >::iterator _row = aobasis.firstShell(); _row != aobasis.lastShell() ; _row++ ) {
                
                 AOShell* _shell_row = aobasis.getShell( _row );
                 
                 if ( _shell_row->getIndex() == atidx ){

                     end   += _shell_row->getNumFunc();
                 } else {
                  
                     idxstart.push_back(start);
                     idxstop.push_back(end);
                     atidx++;
                     start = end;
                     end   += _shell_row->getNumFunc();
                 }
            }

            idxstart.push_back(start);
            idxstop.push_back(end);
                     
            int aidx=0;
            
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                int _a_start = idxstart[aidx];
                int _a_stop  = idxstop[aidx];
                
                double range_max = 0.0;
                //double shiftm_2g = 0.0;
                
                
                // get preset values for this atom type
                double exp_iat = _element_ranges.at((*ait)->type).alpha;
                int    l_iat   = _element_ranges.at((*ait)->type).l;
                
                double x_a = (*ait)->x * 1.8897259886 ;
                double y_a = (*ait)->y * 1.8897259886 ;
                double z_a = (*ait)->z * 1.8897259886 ;
                
                //int iat_diff;
                string type_diff;
                int bidx = 0;
                for (bit = _atoms.begin(); bit < _atoms.end(); ++bit) {
                
                      int _b_start = idxstart[bidx];
                      int _b_stop  = idxstop[bidx];
                      double x_b = (*bit)->x * 1.8897259886 ;
                      double y_b = (*bit)->y * 1.8897259886 ;
                      double z_b = (*bit)->z * 1.8897259886 ;
                      // now do some overlap gymnastics
                      double s_max = 10.0;
                      if ( aidx != bidx ){
                          
                          // find overlap block of these two atoms
                          ub::matrix<double> _overlapblock = ub::project( _overlap._aomatrix ,  ub::range( _a_start , _a_stop ), ub::range(_b_start, _b_stop ) );
                          // determine abs max of this block
                          s_max = 0.0;
                          for ( unsigned i = 0 ; i < _overlapblock.size1(); i++ ){
                          for ( unsigned j = 0 ; j < _overlapblock.size2(); j++ ){
                              s_max = std::max(s_max,std::abs(_overlapblock(i,j)));
                          }
                          }
                      }
                          
                          if ( s_max > 1e-5 ) {
                              // cout << " Atom " << aidx << " neighbor " << bidx << " s-max " << s_max << " exponents " << exp_iat << " and " << _element_ranges.at((*bit)->type).alpha << endl;
                              double range = DetermineCutoff( exp_iat + _element_ranges.at((*bit)->type).alpha, l_iat + _element_ranges.at((*bit)->type).l +2 , eps);
                              // now do some update trickery from Gaussian product formula
                              double dist = sqrt( (x_a - x_b)*(x_a-x_b) + (y_a - y_b)*(y_a-y_b) + (z_a - z_b)*(z_a-z_b)   );
                              double shift_2g = dist*exp_iat/(exp_iat + _element_ranges.at((*bit)->type).alpha );
                              range += shift_2g;
                              
                      
                              if ( aidx != bidx ) range += dist;
                              
                              if ( range > range_max ){
                                  //shiftm_2g = shift_2g;
                                  range_max = range;
                                  //iat_diff = bidx;
                                  type_diff = (*bit)->type;
                              }
                              
             
                              
                          }
                      

                      bidx++;
                }
                
                
                                      
                          
                      // cout << "after atoms :" << endl;
                      //cout << " zmins:     " << exp_iat  << " ---- " << _element_ranges.at(type_diff).alpha << endl;
                      //cout << " lprod:     " << l_iat + _element_ranges.at(type_diff).l << endl;
                      //cout << " range:     " << range_max << endl;
                      //cout << " shiftm_2g: " << shiftm_2g << endl << endl;
                      //exit(0);
                      
                      
                    
                
                if ( round(range_max) > _element_ranges.at((*ait)->type).range ){
                    _element_ranges.at((*ait)->type).range = round(range_max) ;
                    
                }
                
                
                
                aidx++;
            }
            
            
            
            for ( it = _element_ranges.begin() ; it != _element_ranges.end() ; ++it){
                
            //    cout << "Element " << it->first << " alpha " << it->second.alpha << " l " << it->second.l << " Rcut " << it->second.range <<  " Rcut (Ang) " <<  it->second.range  * 0.529177249 << endl;
                
            }
            
        } // getRadialCutoffs
    
    
    
    
    
    
    void EulerMaclaurinGrid::getRadialGrid(BasisSet* bs , vector<QMAtom* > _atoms, string type, GridContainers& _grid) {

        
            map<string, min_exp>::iterator it;
            getRadialCutoffs(_atoms,bs,type);
            
            // go through all elements
            for ( it = _element_ranges.begin() ; it != _element_ranges.end() ; ++it){
                
                 // cout << "Element " << it->first << " alpha " << it->second.alpha << " l " << it->second.l << " Rcut " << it->second.range <<  " Rcut (Ang) " <<  it->second.range  * 0.529177249 << endl;
                
                std::vector<double> points;
                std::vector<double> weights;
                int numberofpoints = getGrid(it->first, type);
                //cout << " Setting grid for element " << it->first <<  " with " << numberofpoints << " points " <<  endl ;
                setGrid( numberofpoints, it->second.range, points, weights );
                
                //grid_element this_element;
                //this_element.gridpoint = points;
                //this_element.weight = weights;
                
                //_element_grids[it->first] = this_element;
                
                _grid._radial_grids[it->first].radius = points; 
                _grid._radial_grids[it->first].weight = weights;
                
            }
              
    }
    
    
    void EulerMaclaurinGrid::setGrid(int np, double cutoff, std::vector<double>& point, std::vector<double>& weight ){
        
        double alpha = -cutoff/(log(1.0 - pow(  (1.0 + double(np))/(2.0 + double(np)),3)) );
        
        
        double factor = 3.0/(1.0+double(np));
        
        for ( int i = 0; i < np; i++){
            double q = double(i+1)/(double(np)+1.0);
            double r = -alpha*log(1.0-pow(q,3));
            double w = factor * alpha * r*r/( 1.0 - pow(q,3) ) * pow(q,2);
            
            point.push_back(r);
            weight.push_back(w);
            
        }
        
        
        
        
    }
    
    
    
    
    double EulerMaclaurinGrid::DetermineCutoff(double alpha, int l, double eps){
        
      // determine norm of function              
                                                                                                                                                                                                     
     /* For a function f(r) = r^k*exp(-alpha*r^2) determine                                                                                                                                                          
        the radial distance r such that the fraction of the                                                                                                                                                          
        function norm that is neglected if the 3D volume                                                                                                                                                             
        integration is terminated at a distance r is less                                                                                                                                                            
        than or equal to eps. */                                                                                                                                                                                       
   
        
        double _cutoff    = 1.0; // initial value
        double _increment = 0.5; // increment

            while (_increment > 0.01) {
                double _neglected = getNeglected(alpha,l ,  _cutoff);
                // cout << "neglected is " << _neglected << endl;
                if (_neglected > eps) {
                    _cutoff += _increment;
                } else {
                    _cutoff -= _increment;
                    if (_cutoff < 0.0) _cutoff = 0.0;
                    _increment = 0.5 * _increment;
                    _cutoff += _increment;
                }
            }

        return _cutoff;
       
    }
    
    
    
    
    double EulerMaclaurinGrid::getNeglected(double alpha, int l, double cutoff){
        
        return RadialIntegral(alpha,l+2,cutoff) / RadialIntegral(alpha,l+2, 0.0); 
        
        
    }


    double EulerMaclaurinGrid::RadialIntegral(double alpha, int l, double cutoff){
        
        const double pi = boost::math::constants::pi<double>();
        int ilo = l % 2;
        double value = 0.0;
        double valexp;
        if ( ilo == 0 ){
            double expo = sqrt(alpha)*cutoff;
            if ( expo > 40.0 ) {
                value = 0.0;
            } else {
                value = 0.5 * sqrt(  pi /alpha  ) * erfc(expo);
            }
        }
        
        double exponent = alpha*cutoff*cutoff;
        if ( exponent > 500.0 ) {
            valexp = 0.0;
            value = 0.0;
        } else {
            valexp = exp(-exponent);
            value = valexp/2.0/alpha;
        }
            
        
        for (  int i = ilo+2; i <= l; i+=2){
            value = ((i-1)*value + pow(cutoff,i-1)*valexp)/2.0/alpha;
        }
        
        return value;
         
        
        
    }
    

    
    int EulerMaclaurinGrid::getGrid(string element, string type){
        
        if ( type == "medium"){
            
            return MediumGrid.at(element);            
            
        }
        else if ( type == "coarse"){
            
            return CoarseGrid.at(element);            
            
        }
        else if ( type == "xcoarse"){
            
            return XcoarseGrid.at(element);            
            
        }
        else if ( type == "fine"){
            
            return FineGrid.at(element);            
            
        }
        else if ( type == "xfine"){
            
            return XfineGrid.at(element);            
            
        }

        throw std::runtime_error("Grid type "+type+" is not implemented");
        return -1;
        
        
    }
    

}
}
