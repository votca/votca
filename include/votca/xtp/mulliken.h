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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __XTP_MULLIKEN__H
#define	__XTP_MULLIKEN__H


#include <votca/xtp/elements.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/qmatom.h>


/**
* \brief Takes a list of atoms, and the corresponding density matrix and puts out a table of mulliken partial charges
*
* 
* 
*/
using namespace std;
using namespace votca::tools;


namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
class Mulliken{
public:
    
    Mulliken(){};
   ~Mulliken(){};
    
    
    void EvaluateMulliken(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat,AOBasis &basis,BasisSet &bs,  bool _do_transition);
  
   
private:
    
     Elements _elements; 
     
     
 
    
};

void Mulliken::EvaluateMulliken(vector< QMAtom* >& _atomlist, ub::matrix<double> &_dmat,AOBasis &basis,BasisSet &bs,  bool _do_transition){
    AOOverlap _overlap;
    // initialize overlap matrix
    _overlap.Initialize(basis._AOBasisSize);
    // Fill overlap
    _overlap.Fill(&basis);
    
    ub::matrix<double> _prodmat = ub::prod( _dmat, _overlap._aomatrix );
    
    vector < QMAtom* > :: iterator atom;

    int id =0;
    for (atom = _atomlist.begin(); atom < _atomlist.end(); ++atom){
                
    
         // get element type and determine its nuclear charge
         if (!_do_transition){
             (*atom)->charge=_elements.getNucCrgECP((*atom)->type);
             //cout << (*atom)->type << " " << (*atom)->charge << endl;
         }
         else {
             (*atom)->charge=0.0;
         }
         // a little messy, have to use basis set to find out which entries in dmat belong to each atom.
         Element* element = bs.getElement((*atom)->type);
         int nooffunc=0;
         for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
             Shell* shell = (*its);
             nooffunc+=shell->getnumofFunc();
         }
         //cout << id << " "<< id+nooffunc << endl;
         for ( int _i = id ; _i < id+nooffunc; _i++){
                (*atom)->charge -= _prodmat(_i,_i);
        }
         id+=nooffunc;
    }
    //cout << id << " " << _dmat.size1() << endl;


}


}}

#endif	/* ESPFIT_H */
