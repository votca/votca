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
#include <votca/tools/property.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {

/*
 * S, P, or D functions in a Gaussian-basis expansion
 */
class Shell 
{
public:

    // Gaussian described by two constants
    typedef pair<double, double> Gaussian;

    // iterator over pairs (decay constant; contraction coefficient)
    typedef vector< Gaussian* >::iterator GaussianIterator;

    GaussianIterator begin() { return _gaussians.begin(); }
    
    GaussianIterator end(){ return _gaussians.end(); }
    
    // default constructor    
    Shell( string type ) : _type(type) { ; }
    
    // clean up all Gaussians
   ~Shell() { 
       for (GaussianIterator it = _gaussians.begin(); it != _gaussians.end() ; it++ ) {
           delete (*it);
       } 
   }
   
    // returns a constant of a Gaussian
    double getConstant( GaussianIterator it ) { return (*it)->first; }
    
    // returns a contraction coefficient of a Gaussian
    double getContraction( GaussianIterator it ) { return (*it)->second; }
    
    // adds a Gaussian 
    void addGaussian( double decay, double contraction ) {
        Gaussian* _gaussian = new Gaussian(decay, contraction);
        _gaussians.push_back( _gaussian );
    }

private:   
    
    // shell type (S, P, D))
    string _type;

    // vector of pairs of decay constants and contraction coefficients
    vector< Gaussian* > _gaussians;

};

/*
 * A collection of shells associated with a specific element  
 */
class Element 
{   
public:
    Element( string type ) : _type(type) { ; }
    
    typedef vector< Shell* >::iterator ShellIterator;
    ShellIterator begin() { return _shells.begin(); }
    ShellIterator end(){ return _shells.end(); }

    string getType() { return _type; }
    
    Shell* getShell( ShellIterator it ) { return (*it); }
    void addShell( Shell* shell ) { _shells.push_back(shell); }
    
    vector<Shell*> getShells() { return _shells; }
    
private:    
    string _type;     
    vector<Shell*> _shells;
};

/*
 * A collection of elements and shells forms the basis set 
 */
class BasisSet 
{
public:
    
    void Load ( string name );
    
    void AddShell(string element_type, Shell* shell );
 
     vector<Shell*> getShells( string element_type ) {
         map<string,Element*>::iterator itm = _elements.find( element_type );
         Element* element = (*itm).second;
         vector<Shell*> shells = element->getShells();
         return shells; 
     }
    
    ~BasisSet();
    
private:    
    
    map<string,Element*> _elements;
};


inline void BasisSet::Load ( string name ) {
    
    Property basis_property;
 
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/basis_sets/") + name + string(".xml");
    
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    list<Property*> elementProps = basis_property.Select("basis.element");
        
    for (list<Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) {
        
        string elementName = (*ite)->getAttribute<string>("name");
        list<Property*> shellProps = (*ite)->Select("shell");
        
        for (list<Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) {
            
            string shellType = (*ite)->getAttribute<string>("type");
            Shell *shell = new Shell( shellType );
            
            list<Property*> constProps = (*its)->Select("constant");
            
            for (list<Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) {
                double decay = (*itc)->getAttribute<double>("decay");
                double contraction = (*itc)->getAttribute<double>("contraction");
                
                shell->addGaussian(decay, contraction);
                
                cout << decay << " " << contraction << endl;
            }
            AddShell(elementName, shell);
        }
       
    }

}

// adding a shell to a basis set
inline void BasisSet::AddShell(string element_type, Shell* shell ) {

   Element* element ;    
    map<string,Element*>::iterator it = _elements.find( element_type ) ;
    
    if ( it == _elements.end() ) {
        element = new Element( element_type );
    } else {
        element = (*it).second;
    }

    element->addShell( shell );
     
    _elements[element_type] =  element;
};

// cleanup the basis set
inline BasisSet::~BasisSet() {
    
    for ( map<string,Element*>::iterator it = _elements.begin(); it !=  _elements.end(); it++ ) {

         vector<Shell*> shells = (*it).second->getShells();

         for ( vector<Shell*>::iterator its = shells.begin(); its !=  shells.end(); its++ ) {
             delete *its;
         }

         delete (*it).second;
     }
};


}}

#endif	/* BASISSET_H */

