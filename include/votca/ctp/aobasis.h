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

#ifndef __CTP_AOBASIS__H
#define	__CTP_AOBASIS__H

#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <votca/ctp/segment.h>
#include <boost/numeric/ublas/matrix.hpp>
#include "basisset.h"

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {
namespace ub = boost::numeric::ublas;
class AOShell;  
class AOBasis;

// Gaussian function: contraction*exp(-decay*r^2)
class AOGaussianPrimitive 
{
    friend class AOShell;
public:
    double decay;
    double contraction;
    int power; // used in pseudopotenials only
    AOShell* aoshell;
private:
    // private constructor, only a shell can create a primitive
    AOGaussianPrimitive( double _decay, double _contraction, AOShell *_aoshell = NULL ) 
    : decay(_decay), contraction(_contraction), aoshell(_aoshell) { ; }

    AOGaussianPrimitive( int _power, double _decay, double _contraction, AOShell *_aoshell = NULL ) 
    : power(_power), decay(_decay), contraction(_contraction), aoshell(_aoshell) { ; }
};      
    
/*
 * S, P, or D functions in a Gaussian-basis expansion
 */
class AOShell 
{
    //friend class AOElement;
    friend class AOBasis;
public:

    string getType() { return _type; }
    int    getNumFunc() { return _numFunc ;}
    int    getStartIndex() { return _startIndex ;}
    int    getOffset() { return _offset ;}
    
    
    int getLmax(  ) {
        int _lmax;
        if ( _type == "S" ) _lmax = 0;
        if ( _type == "SP" ) _lmax = 1;
        if ( _type == "SPD" ) _lmax = 2;
        if ( _type == "P" ) _lmax = 1;
        if ( _type == "PD" ) _lmax = 2;
        if ( _type == "D" ) _lmax = 2;
        return _lmax;
    }; 
    
    vec getPos() { return _pos; }
    double getScale() { return _scale; }
    
    int getSize() { return _gaussians.size(); }
    
    // iterator over pairs (decay constant; contraction coefficient)
    typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
    GaussianIterator firstGaussian() { return _gaussians.begin(); }
    GaussianIterator lastGaussian(){ return _gaussians.end(); }
   
    // adds a Gaussian 
    AOGaussianPrimitive*  addGaussian( double decay, double contraction ) 
    {
        AOGaussianPrimitive* gaussian = new AOGaussianPrimitive(decay, contraction, this);
        _gaussians.push_back( gaussian );
        return gaussian;
    }

    
private:   

    // only class Element can construct shells    
    AOShell( string type, double scale, int numFunc, int startIndex, int offset, vec pos, AOBasis* aobasis = NULL ) : _type(type), _scale(scale), _numFunc(numFunc), _startIndex(startIndex), _offset(offset), _pos(pos) { ; }
    
    // only class Element can destruct shells
   ~AOShell() 
   { 
       for (vector< AOGaussianPrimitive* >::iterator it = _gaussians.begin(); it != _gaussians.end() ; it++ ) delete (*it); 
       _gaussians.clear();
   }
    
    // shell type (S, P, D))
    string _type;
    // scaling factor
    double _scale;
    // number of functions in shell
    int _numFunc;
    int _startIndex;
    vec _pos;
    int _offset;
     
    AOBasis* _aobasis;

    // vector of pairs of decay constants and contraction coefficients
    vector< AOGaussianPrimitive* > _gaussians;

};

/*
 * A collection of shells associated with a specific element  
 */
class AOBasis 
{   
public:
       
       template< class T >
       static void ReorderMOs(ub::matrix<T> &v, vector<int> const &order )  { 
          // Sanity check
          if ( v.size2() != order.size() ) {
              cerr << "Size mismatch in ReorderMOs" << v.size2() << ":" << order.size() << endl;
              throw std::runtime_error( "Abort!");
              
          }
           for ( int _i_orbital = 0; _i_orbital < v.size1(); _i_orbital++ ){
        
                for ( int s = 1, d; s < order.size(); ++ s ) {
                    for ( d = order[s]; d < s; d = order[d] ) ;
                          if ( d == s ) while ( d = order[d], d != s ) swap( v(_i_orbital,s), v(_i_orbital,d) );
                }
           }
       }
       
      template< class T >   
      static void MultiplyMOs(ub::matrix<T> &v, vector<int> const &multiplier )  { 
          // Sanity check
          if ( v.size2() != multiplier.size() ) {
              cerr << "Size mismatch in MultiplyMOs" << v.size2() << ":" << multiplier.size() << endl;
              throw std::runtime_error( "Abort!");
          }

          for ( int _i_orbital = 0; _i_orbital < v.size1(); _i_orbital++ ){

               for ( int _i_basis = 0; _i_basis < v.size2(); _i_basis++ ){
        
                   v(_i_orbital, _i_basis ) = multiplier[_i_basis] * v(_i_orbital, _i_basis );
                   
               }
               
               
           }
       } 
       
       
       
    
    void AOBasisFill( BasisSet* bs , vector<Segment* > segments);
    int NumFuncShell( string shell );
    int OffsetFuncShell( string shell );
    int getAOBasisSize() {return _AOBasisSize; }
    
    typedef vector< AOShell* >::iterator AOShellIterator;
    AOShellIterator firstShell() { return _aoshells.begin(); }
    AOShellIterator lastShell(){ return _aoshells.end(); }

    // string getType() { return _type; }
    // int    getNumFunc() { return _numFunc ;}
        
    AOShell* getShell( AOShellIterator it ) { return (*it); }
    
    AOShell* addShell( string shellType, double shellScale, int shellFunc, int startIndex, int offset, vec pos ) 
    { 
        AOShell* aoshell = new AOShell( shellType, shellScale, shellFunc, startIndex, offset, pos, this );
        _aoshells.push_back(aoshell); 
        return aoshell;
    }
    
    vector<AOShell*> getShells() { return _aoshells; }
    
    AOBasis( ) { ; }
   ~AOBasis() 
   { 
       for (vector< AOShell* >::iterator it = _aoshells.begin(); it != _aoshells.end() ; it++ ) delete (*it); 
       _aoshells.clear();
   }    
   
   int _AOBasisSize;
   
   bool _is_stable;
   
    vector<AOShell*> _aoshells;

    void getReorderVector( string& package, vector<int>& neworder );
    void addReorderShell( string& package, string& shell, vector<int>& neworder );
    void getMultiplierVector( string& package, vector<int>& multiplier );
    void addMultiplierShell( string& package, string& shell, vector<int>& multiplier );  
    
    int getMaxFunctions ( );
};


inline int AOBasis::getMaxFunctions () {
    
    int _maxfunc = 0;
    
    // go through basisset
    for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
        AOShell* _shell = this->getShell( _is );
        int _func_here = _shell->getNumFunc();
        if ( _func_here > _maxfunc ) _maxfunc = _func_here;
    }

    return _maxfunc;
    
}


inline void AOBasis::getMultiplierVector( string& package, vector<int>& multiplier){

    // go through basisset
    for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
        AOShell* _shell = this->getShell( _is );
        string _type =  _shell->getType();
        this->addMultiplierShell(  package, _type, multiplier );
    }
    
}




inline void AOBasis::addMultiplierShell( string& package,  string& shell_type, vector<int>& multiplier ) {
    
    // current length of vector
    int _cur_pos = multiplier.size() -1 ;
    
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ 
           multiplier.push_back( 1 );
       }
       
       if ( shell_type == "P" ){ 
           multiplier.push_back( 1 );
           multiplier.push_back( 1 );
           multiplier.push_back( 1 );
       }
       if ( shell_type == "D" ){ 
           if ( package == "nwchem") {
               multiplier.push_back( -1  ); 
               multiplier.push_back( 1 );
               multiplier.push_back( 1 );
               multiplier.push_back( 1 ); 
               multiplier.push_back( 1 );               
           } else {
               cerr << "Tried to get multipliers d-functions from package " << package << ".";
               throw std::runtime_error( "Multiplication not implemented yet!");
           }
       }
       if ( shell_type == "F" ){ 
           cerr << "Tried to get multipliers for f-functions . ";
           throw std::runtime_error( "Multiplication not implemented yet!");
       }
       if ( shell_type == "G" ){
           cerr << "Tried to get multipliers g-functions . ";
           throw std::runtime_error( "Multiplication not implemented yet!");
       } 
    } else {
        // for combined shells, iterate over all contributions
        //_nbf = 0;
        for( int i = 0; i < shell_type.length(); ++i) {
           string local_shell =    string( shell_type, i, 1 );
           this->addMultiplierShell( package, local_shell, multiplier  );
        }
    }
    
    
}


inline void AOBasis::getReorderVector( string& package, vector<int>& neworder){

    // go through basisset
    for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
        AOShell* _shell = this->getShell( _is );
        string _type =  _shell->getType();
        this->addReorderShell( package,  _type, neworder );
    }
    
}


inline void AOBasis::addReorderShell( string& package,  string& shell_type, vector<int>& neworder ) {
    
    // current length of vector
    int _cur_pos = neworder.size() -1 ;
    
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ 
           neworder.push_back( _cur_pos + 1 );
       }
       
       if ( shell_type == "P" ){ 
           neworder.push_back( _cur_pos + 1 );
           neworder.push_back( _cur_pos + 2 );
           neworder.push_back( _cur_pos + 3 );
       }
       if ( shell_type == "D" ){ 
           if ( package == "gaussian"){
               neworder.push_back( _cur_pos + 4 );
               neworder.push_back( _cur_pos + 1 );
               neworder.push_back( _cur_pos + 2 );
               neworder.push_back( _cur_pos + 5 );
               neworder.push_back( _cur_pos + 3 );
           } else if ( package == "nwchem") {
               neworder.push_back( _cur_pos + 3  ); 
               neworder.push_back( _cur_pos + 2 );
               neworder.push_back( _cur_pos + 4 );
               neworder.push_back( -(_cur_pos + 1) ); // bloody inverted sign
               neworder.push_back( _cur_pos + 5 );               
           } else {
               cerr << "Tried to reorder d-functions from package " << package << ".";
               throw std::runtime_error( "Reordering not implemented yet!");
           }
       }
       if ( shell_type == "F" ){ 
           cerr << "Tried to reorder f-functions . ";
           throw std::runtime_error( "Reordering not implemented yet!");
       }
       if ( shell_type == "G" ){
           cerr << "Tried to reorder g-functions . ";
           throw std::runtime_error( "Reordering not implemented yet!");
       } 
    } else {
        // for combined shells, iterate over all contributions
        //_nbf = 0;
        for( int i = 0; i < shell_type.length(); ++i) {
           string local_shell =    string( shell_type, i, 1 );
           this->addReorderShell( package, local_shell, neworder  );
        }
    }
    
    
}




inline void AOBasis::AOBasisFill(BasisSet* bs , vector<Segment* > segments) {
    
        // cout << "\n Filling AO basis:" << endl;
        
        vector< Atom* > _atoms;
        vector< Atom* > ::iterator ait;
        vector< Segment* >::iterator sit;

       _AOBasisSize = 0;
       _is_stable = true; // _is_stable = true corresponds to gwa_basis%S_ev_stable = .false. 
        
        // loop over segments
        for (sit = segments.begin() ; sit != segments.end(); ++sit) {
        
            _atoms = (*sit)-> Atoms();
            // loop over atoms in segment
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get coordinates of this atom and convert from nm to Bohr
                vec     pos = (*ait)->getQMPos() * 18.897259886;
                // get element type of the atom
                string  name = (*ait)->getElement();
                // get the basis set entry for this element
                Element* element = bs->getElement(name);
                // and loop over all shells
                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                    Shell* shell = (*its);
                    // we don't like contracted basis sets yet
                    if ( shell->getSize() > 1 ) {
                        cerr << "No contracted basis sets!" << flush;
                    } else {
                        AOShell* aoshell = addShell( shell->getType(), shell->getScale(), NumFuncShell( shell->getType() ), _AOBasisSize, OffsetFuncShell( shell->getType() ), pos );
                        _AOBasisSize += NumFuncShell( shell->getType() );
                        for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                           GaussianPrimitive* gaussian = *itg;
                           aoshell->addGaussian(gaussian->decay, gaussian->contraction);
                        }
                    }
                }
            }
        }
         //cout << "Atomic orbitals basis set size: " << _AOBasisSize << endl;
    
    
}

   
    inline int AOBasis::NumFuncShell( string shell_type ) {
    int _nbf;
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 1;}
       if ( shell_type == "P" ){ _nbf = 3;}
       if ( shell_type == "D" ){ _nbf = 5;}
       if ( shell_type == "F" ){ _nbf = 7;}
       if ( shell_type == "G" ){ _nbf = 9;}
    } else {
        // for combined shells, sum over all contributions
        _nbf = 0;
        for( int i = 0; i < shell_type.length(); ++i) {
            string local_shell =    string( shell_type, i, 1 );
            _nbf += this->NumFuncShell( local_shell  );
        }
    }

    return _nbf;
}

    inline int AOBasis::OffsetFuncShell( string shell_type ) {
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
        for(int i = 0; i < shell_type.length(); ++i) {
            string local_shell = string( shell_type, i, 1 );
            int _test = this->OffsetFuncShell( local_shell  );
            if ( _test < _nbf ) { _nbf = _test;} 
        }   
    }
    return _nbf;
}

/*
 * A collection of elements and shells forms the basis set 
 */
/*class AOBasis 
{
public:
    
    void AOLoadBasisSet ( string name );
    
    AOElement* addElement(string elementType );
 
    AOElement* getElement( string element_type ) {
        
         map<string,AOElement*>::iterator itm = _aoelements.find( element_type );
         
         if ( itm == _aoelements.end() ) throw std::runtime_error( "Basis set does not have element of type " + element_type );
         
         AOElement* aoelement = (*itm).second;
         return aoelement; 
     }
    
    ~AOBasis();
    
private:    
    
    map<string,AOElement*> _aoelements;
};


inline void AOBasis::AOLoadBasisSet ( string name ) 
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
        AOElement *aoelement = addElement( elementName );
        //cout << "\nElement " << elementName;
        
        list<Property*> shellProps = (*ite)->Select("shell");
        for (list<Property*> ::iterator  its = shellProps.begin(); its != shellProps.end(); ++its) 
        {            
            string shellType = (*its)->getAttribute<string>("type");
            double shellScale = (*its)->getAttribute<double>("scale");
            
            AOShell* aoshell = aoelement->addShell( shellType, shellScale );
            //cout << "\n\tShell " << shellType;
            
            list<Property*> constProps = (*its)->Select("constant");
            for (list<Property*> ::iterator  itc = constProps.begin(); itc != constProps.end(); ++itc) 
            {
                double decay = (*itc)->getAttribute<double>("decay");
                double contraction = (*itc)->getAttribute<double>("contraction");
                aoshell->addGaussian(decay, contraction);
                //cout << "\n\t\t" << decay << " " << contraction << endl;
            }
            
        }
       
    }
    
}

// adding an Element to a Basis Set
inline AOElement* AOBasis::addElement( string elementType ) {
    AOElement *aoelement = new AOElement( elementType );
    _aoelements[elementType] = aoelement;
    return aoelement;
};



// cleanup the basis set
inline AOBasis::~AOBasis() {
    
    for ( map< string,AOElement* >::iterator it = _aoelements.begin(); it !=  _aoelements.end(); it++ ) {
         delete (*it).second;
     }
    
    _aoelements.clear();
};

 * */

}}

#endif	/* AOBASIS_H */

