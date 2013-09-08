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

#ifndef __CTP_BASISSET__H
#define	__CTP_BASISSET__H

#include <assert.h> 
#include <string>
#include <map>
#include <vector>
#include <utility>

#include <votca/tools/property.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {

/*
 * Container for AO types and coefficients in a
 * sum ( _contraction_coefficients*exp( -_decay_constants*r ) )
 */
class Shell {

public:

    enum Type{S,P,D};

    Shell( Type type, double fadge_factor = 1 ) : _type(type), _fadge_factor(fadge_factor) { ; }

    void addConstant( double decay_constant, double contraction_coefficient ) 
               { 
                 _decay_constants.push_back( decay_constant );
                 _contraction_coefficients.push_back( contraction_coefficient );
               }
    Type getType(){ return _type; }

private:

    Type _type;
    double _fadge_factor;
    vector<double> _decay_constants;
    vector<double> _contraction_coefficients;
};

/*
 * Container for shells for a specific atom type 
 */
class Element {

public:

    enum Type{H,B,C,N,O,F,Al,Si,P,S,Cl};

    void addShell( Shell* shell ) { _shells.push_back(shell); }

    vector<Shell*> getShells () { return _shells; }

private:

    Type _type;     
    vector<Shell*> _shells;

};


/*
 * Container for all atom types 
 */
class BasisSet {

public:
    
    void Load ( string filename );
    
    void AddShell(Element::Type element_type, Shell* shell) { ; }
 
     vector<Shell*> getShells( Element::Type element_type ) {
         map<Element::Type,Element*>::iterator itm = _elements.find( element_type );
         Element* element = (*itm).second;
         vector<Shell*> shells = element->getShells();
         return shells; 
     }
    
    ~BasisSet();
    
private:    
     
    map<Element::Type, Element*> _elements;
};


inline void BasisSet::Load ( string basis_set_name ) {
    
    Property basis_property;
 
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/basis_sets/") + basis_set_name + string(".xml");
    
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    list<Property*> elementProps = basis_property.Select("basis.element");
        
    for (list<Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) {
        list<Property*> shellProps = (*ite)->Select("shell");
        
        for (list<Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) {
            list<Property*> constProps = (*its)->Select("constant");
            
            for (list<Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) {
                double decay = (*itc)->getAttribute<double>("decay");
                double contraction = (*itc)->getAttribute<double>("contraction");
                cout << decay << " " << contraction;
            }
            
        }
       
    }

}

// adding a shell to a basis set
/*
inline void BasisSet::AddShell(Element::ElementType element_type, Shell* shell) {

    //check if the number of constants and coefficients is the same
    assert( decay_constants.size() == contraction_coefficients.size() );

    // create a new shell and fill it up
    Shell* shell = new Shell( shell_type );
     
    for ( vector<double>::iterator it = decay_constants.begin(); it !=  decay_constants.end(); it++ ) {
       shell->_decay_constants.push_back(*it);
    }

    for ( vector<double>::iterator it = contraction_coefficients.begin(); it !=  contraction_coefficients.end(); it++ ) {
       shell->_contraction_coefficients.push_back(*it);
    }
    
    Element* element ;    
    map<string,Element*>::iterator it = _elements.find( element_type ) ;
    
    if ( it == _elements.end() ) {
        element = new Element();
    } else {
        element = (*it).second;
    }

    element->_type = element_type;
    element->_shells.push_back(shell);
     
    _elements[element_type] =  element;
};
*/

// cleanup the basis set
inline BasisSet::~BasisSet() {
    for ( map<Element::Type,Element*>::iterator it = _elements.begin(); it !=  _elements.end(); it++ ) {

         vector<Shell*> shells = (*it).second->getShells();

         for ( vector<Shell*>::iterator its = shells.begin(); its !=  shells.end(); its++ ) {
             delete *its;
         }

         delete (*it).second;
     }
};


}}

#endif	/* BASISSET_H */

