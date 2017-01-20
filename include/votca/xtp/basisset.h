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

#ifndef __XTP_BASISSET__H
#define	__XTP_BASISSET__H

#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>


using namespace votca::tools;

namespace votca { namespace xtp {
    //namespace CTP = votca::ctp;

class Shell;  
class Element;
class BasisSet;

// Gaussian function: contraction*exp(-decay*r^2)
class GaussianPrimitive 
{
    friend class Shell;
public:
    int power; // used in pseudopotenials only
    double decay;
    std::vector<double> contraction;
    Shell* shell;
private:
    // private constructor, only a shell can create a primitive
    GaussianPrimitive( double _decay, std::vector<double> _contraction, Shell *_shell = NULL ) 
    : decay(_decay),
    contraction(_contraction),
    shell(_shell) { ; }

    GaussianPrimitive( int _power, double _decay, std::vector<double> _contraction, Shell *_shell = NULL ) 
    : power(_power),
    decay(_decay),
    contraction(_contraction),
    shell(_shell) { ; }
};      
    
/*
 * S, P, or D functions in a Gaussian-basis expansion
 */
class Shell 
{
    friend class Element;   
public:

    std::string getType() { return _type; }

  /*  int getLmax(  ) {
        int _lmax;
        if ( _type == "S" ) _lmax = 0;
        if ( _type == "SP" ) _lmax = 1;
        if ( _type == "SPD" ) _lmax = 2;
        if ( _type == "P" ) _lmax = 1;
        if ( _type == "PD" ) _lmax = 2;
        if ( _type == "D" ) _lmax = 2;
        if ( _type == "F" ) _lmax = 3;
        if ( _type == "SPDF" ) _lmax = 3;
        return _lmax;
    }; */
    
    
    bool combined(){
        if (_type.length()>1){
            return true;
        }
        return false;
    }
    
   
    
    int getLmax(  ) {
        return FindLmax(_type);
    }
    
    int FindLmax( string _type ){

        int _lmax;
        // single type shells
        if ( _type.length() == 1 ){
            if ( _type == "S" ){ _lmax = 0;}
            if ( _type == "P" ){ _lmax = 1;}
            if ( _type == "D" ){ _lmax = 2;}
            if ( _type == "F" ){ _lmax = 3;}
            if ( _type == "G" ){ _lmax = 4;}
        } else {
            _lmax = -1;
            for(unsigned i = 0; i < _type.length(); ++i) {
                string local_shell = string( _type, i, 1 );
                int _test = this->FindLmax( local_shell  );
                if ( _test > _lmax ) { _lmax = _test;} 
            }   
        }

        return _lmax;
        
    }
    
    
    
    
   
    int OffsetFuncShell( string shell_type ) {
    int _nbf;
    // single type shells
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 0;}
       if ( shell_type == "P" ){ _nbf = 1;}
       if ( shell_type == "D" ){ _nbf = 4;}
       if ( shell_type == "F" ){ _nbf = 9;}
       if ( shell_type == "G" ){ _nbf = 16;}
    } else {
        // for combined shells, go over all contributions and find minimal offset
        _nbf = 1000;
        for(unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell = string( shell_type, i, 1 );
            int _test = this->OffsetFuncShell( local_shell  );
            if ( _test < _nbf ) { _nbf = _test;} 
        }   
    }
    return _nbf;
        } 
    
    
    
    
    
    
    
    int getnumofFunc() {
        return FindnumofFunc(_type);
    }; 
    
    
    int FindnumofFunc( string shell_type) {
        int _nbf = 0;
        // single type shells
        if ( shell_type.length() == 1 ){
            if ( shell_type == "S" ){ _nbf = 1;}
            if ( shell_type == "P" ){ _nbf = 3;}
            if ( shell_type == "D" ){ _nbf = 5;}
            if ( shell_type == "F" ){ _nbf = 7;}
            if ( shell_type == "G" ){ _nbf = 9;}
    } else {
        // for combined shells, go over all contributions and add functions
        _nbf = 0;
        
           for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string local_shell = string(shell_type, i, 1);
                    _nbf += this->FindnumofFunc(local_shell);
            }   
    }
    return _nbf;
        
    }
    
    
    
    double getScale() { return _scale; }
    
    int getSize() { return _gaussians.size(); }
    
    // iterator over pairs (decay constant; contraction coefficient)
    typedef vector< GaussianPrimitive* >::iterator GaussianIterator;
    GaussianIterator firstGaussian() { return _gaussians.begin(); }
    GaussianIterator lastGaussian(){ return _gaussians.end(); }
   
    // adds a Gaussian 
    GaussianPrimitive*  addGaussian( double decay, std::vector<double> contraction ) 
    {
        GaussianPrimitive* gaussian = new GaussianPrimitive(decay, contraction, this);
        _gaussians.push_back( gaussian );
        return gaussian;
    }

    // adds a Gaussian of a pseudopotential
    GaussianPrimitive*  addGaussian( int power, double decay, std::vector<double> contraction ) 
    {
        GaussianPrimitive* gaussian = new GaussianPrimitive(power, decay, contraction, this);
        _gaussians.push_back( gaussian );
        return gaussian;
    }    
    
private:   

    // only class Element can construct shells    
    Shell( std::string type, double scale, Element* element = NULL ) : _type(type), _scale(scale) { ; }
    
    // only class Element can destruct shells
   ~Shell() 
   { 
       for (vector< GaussianPrimitive* >::iterator it = _gaussians.begin(); it != _gaussians.end() ; it++ ) delete (*it); 
       _gaussians.clear();
   }
    
    // shell type (S, P, D))
    std::string _type;
    // scaling factor
    double _scale;
     
    //Element* _element;

    // vector of pairs of decay constants and contraction coefficients
    vector< GaussianPrimitive* > _gaussians;

};

/*
 * A collection of shells associated with a specific element  
 */
class Element 
{   
    friend class BasisSet;
public:
    
    typedef vector< Shell* >::iterator ShellIterator;
    ShellIterator firstShell() { return _shells.begin(); }
    ShellIterator lastShell(){ return _shells.end(); }

    std::string getType() { return _type; }
    
    int getLmax() { return _lmax; }
    
    int getNcore() { return _ncore; }
    
    Shell* getShell( ShellIterator it ) { return (*it); }
    
    Shell* addShell( std::string shellType, double shellScale ) 
    { 
        Shell* shell = new Shell( shellType, shellScale, this );
        _shells.push_back(shell); 
        return shell;
    }
    
    vector<Shell*> getShells() { return _shells; }
    
private:  
    
    // only class BasisSet can create Elements
    Element( std::string type ) : _type(type) { ; }

    // used for the pseudopotential
    Element( std::string type, int lmax, int ncore ) : _type(type), _lmax(lmax), _ncore(ncore)  { ; }
    
    // only class BasisSet can destruct Elements
   ~Element() 
   { 
       for (vector< Shell* >::iterator it = _shells.begin(); it != _shells.end() ; it++ ) delete (*it); 
       _shells.clear();
   }    
   
    std::string _type;    
    // lmax is used in the pseudopotentials only (applies to the highest angular momentum lmax)
    int _lmax;
    // ncore is used in the pseudopotentials only (replaces ncore electrons))
    int _ncore;
    
    vector<Shell*> _shells;    
};

/*
 * A collection of elements and shells forms the basis set 
 */
class BasisSet 
{
public:
    
    void LoadBasisSet ( std::string name );

    void LoadPseudopotentialSet ( std::string name );
    
    Element* addElement(std::string elementType );
    
    // used for pseudopotentials only
    Element* addElement(std::string elementType, int lmax, int ncore );
 
    Element* getElement( std::string element_type ) {
        
         map<std::string,Element*>::iterator itm = _elements.find( element_type );
         
         if ( itm == _elements.end() ) throw std::runtime_error( "Basis set does not have element of type " + element_type );
         
         Element* element = (*itm).second;
         return element; 
     }
    
    ~BasisSet();
    
private:    
    
    map<std::string,Element*> _elements;
};


inline void BasisSet::LoadBasisSet ( std::string name ) 
{    
    Property basis_property;
 
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    std::string xmlFile = std::string(getenv("VOTCASHARE")) + std::string("/xtp/basis_sets/") + name + std::string(".xml");
    
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    list<Property*> elementProps = basis_property.Select("basis.element");
        
    for (list<Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) 
    {       
        std::string elementName = (*ite)->getAttribute<std::string>("name");
        Element *element = addElement( elementName );
        //cout << "\nElement " << elementName;
        
        list<Property*> shellProps = (*ite)->Select("shell");
        for (list<Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) 
        {            
            std::string shellType = (*its)->getAttribute<std::string>("type");
            double shellScale = (*its)->getAttribute<double>("scale");
            
            Shell* shell = element->addShell( shellType, shellScale );
            //cout << "\n\tShell " << shellType;
            
            list<Property*> constProps = (*its)->Select("constant");
            for (list<Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) 
            {
                double decay = (*itc)->getAttribute<double>("decay");
                // << " decay "<<decay<<endl;
                std::vector<double> contraction;
                contraction.resize(shell->getLmax()+1); 
                list<Property*> contrProps = (*itc)->Select("contractions");
                for (list<Property*> ::iterator itcont = contrProps.begin(); itcont != contrProps.end(); ++itcont)
                {
                    std::string contrType = (*itcont)->getAttribute<std::string>("type");
                    double contrFactor = (*itcont)->getAttribute<double>("factor");
                    //cout << " factor " << contrFactor << endl;
                    if ( contrType == "S" ) contraction[0] = contrFactor;
                    if ( contrType == "P" ) contraction[1] = contrFactor;
                    if ( contrType == "D" ) contraction[2] = contrFactor;
                    if ( contrType == "F" ) contraction[3] = contrFactor;
                    if ( contrType == "G" ) contraction[4] = contrFactor;
                }    
                shell->addGaussian(decay, contraction);
            }
            
        }
       
    }
    
}


inline void BasisSet::LoadPseudopotentialSet ( std::string name ) 
{    
    Property basis_property;
 
    // get the path to the shared folders with xml files
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    std::string xmlFile = std::string(getenv("VOTCASHARE")) + std::string("/xtp/basis_sets/") + name + std::string(".xml");
    
    bool success = load_property_from_xml(basis_property, xmlFile);
    
    if ( !success ) {; }
    
    list<Property*> elementProps = basis_property.Select("pseudopotential.element");
        
    for (list<Property*> ::iterator  ite = elementProps.begin(); ite != elementProps.end(); ++ite) 
    {       
        std::string elementName = (*ite)->getAttribute<std::string>("name");
        int lmax = (*ite)->getAttribute<int>("lmax");
        int ncore = (*ite)->getAttribute<int>("ncore");
        
        Element *element = addElement( elementName, lmax, ncore );
        //cout << "\nElement " << elementName;
        
        list<Property*> shellProps = (*ite)->Select("shell");
        for (list<Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) 
        {            
            std::string shellType = (*its)->getAttribute<std::string>("type");
            double shellScale = 1.0;
            
            Shell* shell = element->addShell( shellType, shellScale );
            //cout << "\n\tShell " << shellType;
            
            list<Property*> constProps = (*its)->Select("constant");
            for (list<Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) 
            {
                int power = (*itc)->getAttribute<int>("power");
                double decay = (*itc)->getAttribute<double>("decay");
                //double contraction = (*itc)->getAttribute<double>("contraction");
                std::vector<double> contraction;
                // just testing here with single value 
                contraction.push_back((*itc)->getAttribute<double>("contraction"));
                //shell->addGaussian(decay, contraction);
                shell->addGaussian(power, decay, contraction);
                //cout << "\n\t\t" << decay << " " << contraction << endl;
            }
            
        }
       
    }

}


// adding an Element to a Basis Set
inline Element* BasisSet::addElement( std::string elementType ) {
    Element *element = new Element( elementType );
    _elements[elementType] = element;
    return element;
};

// adding an Element to a Pseudopotential Library
inline Element* BasisSet::addElement( std::string elementType, int lmax, int ncore ) {
    Element *element = new Element( elementType, lmax, ncore );
    _elements[elementType] = element;
    return element;
};

// cleanup the basis set
inline BasisSet::~BasisSet() {
    
    for ( map< std::string,Element* >::iterator it = _elements.begin(); it !=  _elements.end(); it++ ) {
         delete (*it).second;
     }
    
    _elements.clear();
};


}}

#endif	/* BASISSET_H */

