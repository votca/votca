/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include <boost/math/constants/constants.hpp>
#include "votca/xtp/radial_euler_maclaurin_rule.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/qmatom.h"
#include <votca/xtp/aomatrix.h>




namespace votca { namespace xtp { 
    
    
    std::vector<double> EulerMaclaurinGrid::getPruningIntervals(const std::string & element){
        
        std::vector<double> _r;
        
        // get Bragg-Slater Radius for this element
        double BSradius = _BraggSlaterRadii.at(element);
        
        // row type of element
        int RowType = _pruning_set.at(element);
        
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
            
            std::cerr << "Pruning unsupported for RowType " << RowType << std::endl;
            exit(1);
        }

        return _r;
     }
    
     
void EulerMaclaurinGrid::getRadialCutoffs(const AOBasis* aobasis,std::vector<QMAtom* > _atoms, const std::string& gridtype) {

            std::map<std::string, min_exp>::iterator it;
            std::vector< QMAtom* > ::iterator ait;
            std::vector< QMAtom* > ::iterator bit;

            double eps = Accuracy[gridtype];
            double _decaymin;
            int _lvalue;
            //  cout << endl << " Setting cutoffs for grid type " << gridtype << " eps = " << eps << endl;
            // 1) is only element based
            // loop over atoms

            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get element type of the atom
                std::string name = (*ait)->getType();
                // is this element already in map?
                it = _element_ranges.find(name);
                // only proceed, if element data does not exist yet
                if (it == _element_ranges.end()) {
                    // get first range estimate and add to map
                    min_exp this_atom;
                    double range_max = 0.0;
                    const std::vector<const AOShell*>  shells=aobasis->getShellsperAtom((*ait)->getAtomID());
                    std::vector<const AOShell*>::const_iterator its;
                    // and loop over all shells to figure out minimum decay constant and angular momentum of this function
                    for (its =shells.begin(); its !=shells.end() ; its++) {
                        int _lmax = (*its)->getLmax();
                       (*its)->getMinDecay();
                        _decaymin = 1e7;
                        _lvalue = 0;
                       
                       
                            if ((*its)->getMinDecay() < _decaymin) {
                                _decaymin = (*its)->getMinDecay();
                                _lvalue = _lmax;
                            }
                        
                        double range = DetermineCutoff(2 * _decaymin, 2 * _lvalue + 2, eps);
                        if (range > range_max) {
                            this_atom.alpha = _decaymin;
                            this_atom.l = _lvalue;
                            this_atom.range = range;
                            range_max = range;
                        }
                    } // shells

                    //   cout << "Element " << name << " alpha " << this_atom.alpha << " l " << this_atom.l << " Rcut " << this_atom.range << endl;
                    _element_ranges[name] = this_atom;
                } // new element
            } // atoms

            // calculate overlap matrix
           
            AOOverlap _overlap;
            // Fill overlap
            _overlap.Fill(*aobasis);



            // refining by going through all atom combinations
            // get collapsed index list
            int atidx = 0;
            std::vector<int> idxstart;
            std::vector<int> idxstop;
            int start = 0;
            int end = 0;
            for (AOBasis::AOShellIterator _row = aobasis->firstShell(); _row != aobasis->lastShell(); _row++) {

                const AOShell* _shell_row = *_row;

                if (_shell_row->getIndex() == atidx) {

                    end += _shell_row->getNumFunc();
                } else {

                    idxstart.push_back(start);
                    idxstop.push_back(end);
                    atidx++;
                    start = end;
                    end += _shell_row->getNumFunc();
                }
            }

            idxstart.push_back(start);
            idxstop.push_back(end);

            int aidx = 0;

            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {

                int _a_start = idxstart[aidx];
                int _a_stop = idxstop[aidx];

                double range_max = 0.0;

                // get preset values for this atom type
                double exp_iat = _element_ranges.at((*ait)->getType()).alpha;
                int l_iat = _element_ranges.at((*ait)->getType()).l;

                const tools::vec& pos_a = (*ait)->getPos();


                
                int bidx = 0;
                for (bit = _atoms.begin(); bit < _atoms.end(); ++bit) {

                    int _b_start = idxstart[bidx];
                    int _b_stop = idxstop[bidx];
                    const tools::vec& pos_b = (*bit)->getPos();
                    // now do some overlap gymnastics
                    double s_max =  std::numeric_limits<double>::max();
                    if (aidx != bidx) {

                        // find overlap block of these two atoms
                        Eigen::MatrixXd _overlapblock = _overlap.Matrix().block(_a_start,_b_start, _a_stop-_a_start,_b_stop-_b_start);
                        // determine abs max of this block
                        s_max = 0.0;
                        for (unsigned i = 0; i < _overlapblock.rows(); i++) {
                            for (unsigned j = 0; j < _overlapblock.cols(); j++) {
                                s_max = std::max(s_max, std::abs(_overlapblock(i, j)));
                            }
                        }
                    }

                    if (s_max > 1e-5) {
                        // cout << " Atom " << aidx << " neighbor " << bidx << " s-max " << s_max << " exponents " << exp_iat << " and " << _element_ranges.at((*bit)->type).alpha << endl;
                        double range = DetermineCutoff(exp_iat + _element_ranges.at((*bit)->getType()).alpha, l_iat + _element_ranges.at((*bit)->getType()).l + 2, eps);
                        // now do some update trickery from Gaussian product formula
                        double dist = abs(pos_b - pos_a);
                        double shift_2g = dist * exp_iat / (exp_iat + _element_ranges.at((*bit)->getType()).alpha);
                        range += shift_2g;


                        if (aidx != bidx) range += dist;

                        if (range > range_max) {
                            range_max = range;                          
                        }

                    }

                    bidx++;
                }


                if (round(range_max) > _element_ranges.at((*ait)->getType()).range) {
                    _element_ranges.at((*ait)->getType()).range = round(range_max);
                }
                aidx++;
            }
           
          return;  
        } // getRadialCutoffs
    
    
    
    
    
    
    void EulerMaclaurinGrid::getRadialGrid(const AOBasis* aobasis , std::vector<QMAtom* > _atoms,const std::string& type, GridContainers& _grid) {

        
            std::map<std::string, min_exp>::iterator it;
            getRadialCutoffs(aobasis,_atoms,type);
            
            // go through all elements
            for ( it = _element_ranges.begin() ; it != _element_ranges.end() ; ++it){
              
                std::vector<double> points;
                std::vector<double> weights;
                int numberofpoints = getGrid(it->first, type);
                setGrid( numberofpoints, it->second.range, points, weights );
                
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
        
        return;
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
                value = 0.5 * sqrt(  pi /alpha  ) * std::erfc(expo);
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
    

    
    int EulerMaclaurinGrid::getGrid(const std::string& element, const std::string& type){
        
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
