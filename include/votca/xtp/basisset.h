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



namespace votca { namespace xtp {
 // shell type (S, P, D))
    int FindLmax(const std::string& _type);

    int FindLmin(const std::string& _type);

    int OffsetFuncShell(const std::string& shell_type);

    int NumFuncShell(const std::string& shell_type);
    int NumFuncShell_cartesian(const std::string& shell_type);

    int OffsetFuncShell_cartesian(const std::string& shell_type);
    
    std::vector<int> NumFuncSubShell(const std::string& shell_type);
    
 

class Shell;  
class Element;
class BasisSet;

// Gaussian function: contraction*exp(-decay*r^2)
class GaussianPrimitive 
{
    friend class Shell;
public:
    int _power; // used in pseudopotenials only
    double _decay;
    std::vector<double> _contraction;
private:
    // private constructor, only a shell can create a primitive
    GaussianPrimitive( double decay, std::vector<double> contraction ) 
    : _decay(decay),
    _contraction(contraction){ ; }

    GaussianPrimitive( int power, double decay, std::vector<double> contraction ) 
    : _power(power),
    _decay(decay),
    _contraction(contraction){ ; }
};      
    

class Shell 
{
    friend class Element;   
public:

    const std::string& getType() const{ return _type; }
    
    bool combined()const{
        if (_type.length()>1){
            return true;
        }
        return false;
    }
    
    int getLmax() const{
        return FindLmax(_type);
    }
    
    int getLmin() const{
        return FindLmin(_type);
    }
    
    int getnumofFunc() const{
        return NumFuncShell(_type);
    }; 

    double getScale() const{ return _scale; }
    
    int getSize() const{ return _gaussians.size(); }
    
    // iterator over pairs (decay constant; contraction coefficient)
    typedef std::vector< GaussianPrimitive* >::iterator GaussianIterator;
    GaussianIterator firstGaussian() { return _gaussians.begin(); }
    GaussianIterator lastGaussian(){ return _gaussians.end(); }
   
    // adds a Gaussian 
    GaussianPrimitive*  addGaussian( double decay, std::vector<double> contraction );
   

    // adds a Gaussian of a pseudopotential
    GaussianPrimitive*  addGaussian( int power, double decay, std::vector<double> contraction );
  
    
    friend std::ostream &operator<<(std::ostream &out, const Shell& shell);
private:   

    // only class Element can construct shells    
    Shell( std::string type, double scale, Element* element = NULL ) : _type(type), _scale(scale) { ; }
    
    // only class Element can destruct shells
   ~Shell() 
   { 
       for (std::vector< GaussianPrimitive* >::iterator it = _gaussians.begin(); it != _gaussians.end() ; it++ ) delete (*it); 
       _gaussians.clear();
   }
    
   
    std::string _type;
    // scaling factor
    double _scale;
     

    // vector of pairs of decay constants and contraction coefficients
    std::vector< GaussianPrimitive* > _gaussians;

};

/*
 * A collection of shells associated with a specific element  
 */
class Element 
{   
    friend class BasisSet;
public:
    
    typedef std::vector< Shell* >::iterator ShellIterator;
     std::vector< Shell* >::iterator begin() { return _shells.begin(); }
     std::vector< Shell* >::iterator end(){ return _shells.end(); }

    std::string getType(){ return _type; }
    
    int getLmax() const{ return _lmax; }
    
    int getNcore() const{ return _ncore; }
    
    Shell* getShell( ShellIterator it ) { return (*it); }
    
    Shell* addShell( std::string shellType, double shellScale ) 
    { 
        Shell* shell = new Shell( shellType, shellScale, this );
        _shells.push_back(shell); 
        return shell;
    }
    
    std::vector<Shell*> getShells() { return _shells; }
    
private:  
    
    // only class BasisSet can create Elements
    Element( std::string type ) : _type(type) { ; }

    // used for the pseudopotential
    Element( std::string type, int lmax, int ncore ) : _type(type), _lmax(lmax), _ncore(ncore)  { ; }
    
    // only class BasisSet can destruct Elements
   ~Element() 
   { 
       for (std::vector< Shell* >::iterator it = _shells.begin(); it != _shells.end() ; it++ ) delete (*it); 
       _shells.clear();
   }    
   
    std::string _type;    
    // lmax is used in the pseudopotentials only (applies to the highest angular momentum lmax)
    int _lmax;
    // ncore is used in the pseudopotentials only (replaces ncore electrons))
    int _ncore;
    
    std::vector<Shell*> _shells;    
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
        
         std::map<std::string,Element*>::iterator itm = _elements.find( element_type );
         
         if ( itm == _elements.end() ) throw std::runtime_error( "Basis set "+_name+" does not have element of type " + element_type );
         
         Element* element = (*itm).second;
         return element; 
     }
    
    std::map< std::string,Element* >::iterator begin() { return _elements.begin(); }
    std::map< std::string,Element* >::iterator end(){ return _elements.end(); }

    
        
    ~BasisSet();
    
private:    
    std::string _name;
    std::map<std::string,Element*> _elements;
};


}}

#endif	/* BASISSET_H */

