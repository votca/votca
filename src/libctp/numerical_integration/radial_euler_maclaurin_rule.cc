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

#include <boost/math/constants/constants.hpp>
# include "votca/ctp/radial_euler_maclaurin_rule.h"

using namespace std;


namespace votca { namespace ctp { 
void EulerMaclaurinGrid::getRadialCutoffs(vector<QMAtom* > _atoms, BasisSet* bs, string gridtype) {

            map<string, min_exp>::iterator it;
            vector< QMAtom* > ::iterator ait;
            
            double eps = Accuracy[gridtype];
            
            // 1) is only element based
            // loop over atoms
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get element type of the atom
                string name = (*ait)->type;
                // is this element already in map?
                it = _element_ranges.find(name);
                // only proceed, if element data does not exist yet
                if (it == _element_ranges.end()) {
                    // get the basis set entry for this element
                    Element* _element = bs->getElement(name);
                    // and loop over all shells to figure out minimum decay constant and angular momentum of this function
                    double _decaymin = 1e7;
                    int _lvalue = 0;
                    for (Element::ShellIterator its = _element->firstShell(); its != _element->lastShell(); its++) {
                        Shell* shell = (*its);
                        int _lmax = shell->getLmax();
                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                            GaussianPrimitive* gaussian = *itg;
                            double _decay = gaussian->decay;
                            if (_decay < _decaymin) {
                                _decaymin = _decay;
                                _lvalue = _lmax;
                            }
                        }

                    } // shells

                    // get first range estimate and add to map
                    min_exp this_atom;
                    this_atom.alpha = _decaymin;
                    this_atom.l     = _lvalue;
                    this_atom.range = DetermineCutoff( _decaymin, _lvalue, eps);;
                    
                    _element_ranges[name] = this_atom;
                } // new element
            } // atoms


        } // getRadialCutoffs
    
    
    
    
    
    
    void EulerMaclaurinGrid::getRadialGrid(BasisSet* bs , string element, string type, std::vector<double>& _r, std::vector<double>& _weight) {

        
            /* apparently, the procedure is like this
             *
             * 1) get all atom coordinates
             * 2) determine overlap matrix
             * 3) set accuracy
             * 4) for each atom type, get min_exp, l(min_exp) and start Rcut estimate
             * 5) for each atom, refine Rcut taking the overlap matrix with other atoms into account, update Rcut for atom type!!
             * 
             */
        
        
        
            int order = getGrid(element, type);
 
            
            // determine radial cutoff for this (_decaymin,_lvalue) combination
            cout << " For " << element << " decay min is " << _decaymin << " of type " << _lvalue << endl;

            double eps = 1.e-6;
            double cutoff = DetermineCutoff( _decaymin, _lvalue, eps);
            cout << " And a cutoff would be " << cutoff << "for eps= " << eps << endl;
            
                                         
              
              
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
                double _neglected = getNeglected(2 * l + 2, 2.0 * alpha, _cutoff);
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
        
        return RadialIntegral(l+2,alpha,cutoff) / RadialIntegral(l+2,alpha, 0.0); 
        
        
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
        
        
    }
    

}
}