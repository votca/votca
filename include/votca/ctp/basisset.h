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

#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {

class Shell;  
class Element;
class BasisSet;

// Gaussian function: contraction*exp(-decay*r^2)
class GaussianPrimitive 
{
    friend Shell;
public:
    double decay;
    double contraction;
    Shell* shell;
private:
    // private constructor, only a shell can create a primitive
    GaussianPrimitive( double _decay, double _contraction, Shell *_shell = NULL ) 
    : decay(_decay), contraction(_contraction), shell(_shell) { ; }
};      
    
/*
 * S, P, or D functions in a Gaussian-basis expansion
 */
class Shell 
{
    friend Element;   
public:

    string getType() { return _type; }
    
    int getSize() { return _gaussians.size(); }
    
    // iterator over pairs (decay constant; contraction coefficient)
    typedef vector< GaussianPrimitive* >::iterator GaussianIterator;
    GaussianIterator firstGaussian() { return _gaussians.begin(); }
    GaussianIterator lastGaussian(){ return _gaussians.end(); }
   
    // adds a Gaussian 
    GaussianPrimitive*  addGaussian( double decay, double contraction ) 
    {
        GaussianPrimitive* gaussian = new GaussianPrimitive(decay, contraction, this);
        _gaussians.push_back( gaussian );
        return gaussian;
    }

private:   

    // only class Element can construct shells    
    Shell( string type, Element* element = NULL ) : _type(type) { ; }
    
    // only class Element can destruct shells
   ~Shell() 
   { 
       for (vector< GaussianPrimitive* >::iterator it = _gaussians.begin(); it != _gaussians.end() ; it++ ) delete (*it); 
       _gaussians.clear();
   }
    
    // shell type (S, P, D))
    string _type;
    
    Element* _element;

    // vector of pairs of decay constants and contraction coefficients
    vector< GaussianPrimitive* > _gaussians;

};

/*
 * A collection of shells associated with a specific element  
 */
class Element 
{   
    friend BasisSet;
public:
    
    typedef vector< Shell* >::iterator ShellIterator;
    ShellIterator firstShell() { return _shells.begin(); }
    ShellIterator lastShell(){ return _shells.end(); }

    string getType() { return _type; }
    
    Shell* getShell( ShellIterator it ) { return (*it); }
    
    Shell* addShell( string shellType ) 
    { 
        Shell* shell = new Shell( shellType, this );
        _shells.push_back(shell); 
        return shell;
    }
    
    vector<Shell*> getShells() { return _shells; }
    
private:  
    
    // only class BasisSet can create Elements
    Element( string type ) : _type(type) { ; }

    // only class BasisSet can destruct Elements
   ~Element() 
   { 
       for (vector< Shell* >::iterator it = _shells.begin(); it != _shells.end() ; it++ ) delete (*it); 
       _shells.clear();
   }    
   
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
    
    Element* addElement(string elementType );
 
    Element* getElement( string element_type ) {
        
         map<string,Element*>::iterator itm = _elements.find( element_type );
         
         if ( itm == _elements.end() ) throw std::runtime_error( "Basis set does not have element of type " + element_type );
         
         Element* element = (*itm).second;
         return element; 
     }
    
    ~BasisSet();
    
private:    
    
    map<string,Element*> _elements;
};


inline void BasisSet::Load ( string name ) 
{    
    Property basis_property;
 
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/basis_sets/") + name + string(".xml");
    
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    list<Property*> elementProps = basis_property.Select("basis.element");
        
    for (list<Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) 
    {       
        string elementName = (*ite)->getAttribute<string>("name");
        Element *element = addElement( elementName );
        //cout << "\nElement " << elementName;
        
        list<Property*> shellProps = (*ite)->Select("shell");
        for (list<Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) 
        {            
            string shellType = (*its)->getAttribute<string>("type");
            Shell* shell = element->addShell( shellType );
            //cout << "\n\tShell " << shellType;
            
            list<Property*> constProps = (*its)->Select("constant");
            for (list<Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) 
            {
                double decay = (*itc)->getAttribute<double>("decay");
                double contraction = (*itc)->getAttribute<double>("contraction");
                shell->addGaussian(decay, contraction);
                //cout << "\n\t\t" << decay << " " << contraction << endl;
            }
            
        }
       
    }

}

// adding an Element to a Basis Set
inline Element* BasisSet::addElement( string elementType ) {
    Element *element = new Element( elementType );
    _elements[elementType] = element;
    return element;
};

// cleanup the basis set
inline BasisSet::~BasisSet() {
    
    for ( map< string,Element* >::iterator it = _elements.begin(); it !=  _elements.end(); it++ ) {
         delete (*it).second;
     }
    
    _elements.clear();
};


}}

#endif	/* BASISSET_H */

