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

#include <votca/tools/property.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/qmatom.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/ctp/basisset.h>


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
    int power; // used in pseudopotenials only
    double decay;
    std::vector<double> contraction;
    AOShell* aoshell;
private:
    // private constructor, only a shell can create a primitive
    AOGaussianPrimitive( double _decay, std::vector<double> _contraction, AOShell *_aoshell = NULL ) 
    : decay(_decay),
            contraction(_contraction),
            aoshell(_aoshell) { ; }

    AOGaussianPrimitive( int _power, double _decay, std::vector<double> _contraction, AOShell *_aoshell = NULL ) 
    : power(_power),
    decay(_decay),
    contraction(_contraction),
    aoshell(_aoshell) { ; }
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
    int    getIndex() { return _atomindex;}
    string getName() { return _atomname;}
    
    int getLmax(  ) { return detlmax( _type );}
    /*
        int _lmax;
        if ( _type == "S" ) _lmax = 0;
        if ( _type == "SP" ) _lmax = 1;
        if ( _type == "SPD" ) _lmax = 2;
        if ( _type == "P" ) _lmax = 1;
        if ( _type == "PD" ) _lmax = 2;
        if ( _type == "D" ) _lmax = 2;
        
        
        return _lmax;
    };*/ 
    
    vec getPos() { return _pos; }
    double getScale() { return _scale; }
    
    int getSize() { return _gaussians.size(); }
    
    
    //vector<double> evalAOspace( double x, double y, double z , string type = "");
    //void EvalAOspace( ub::matrix_range<ub::matrix<double> >& AOvalues, double x, double y, double z , string type = "");
    void EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues, double x, double y, double z );
    void EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues,ub::matrix_range<ub::matrix<double> >& AODervalues, double x, double y, double z );
    //void EvalAOspace(ub::matrix<double>& AOvalues, double x, double y, double z , string type = "");
    
    void EvalAOIntegral(ub::matrix_range<ub::matrix<double> >& AOvalues);
    //vector< vector<double> > evalAOGradspace( double x, double y, double z , string type = "");
    //void EvalAOGradspace( ub::matrix_range<ub::matrix<double> >& AODerXvalues,ub::matrix_range<ub::matrix<double> >& AODerYvalues,ub::matrix_range<ub::matrix<double> >& AODerZvalues, double x, double y, double z , string type = "");
    void EvalAOGradspace( ub::matrix_range<ub::matrix<double> >& AODervalues, double x, double y, double z , string type = "");
    //void EvalAOGradspace( ub::matrix<double>& AODervalues, double x, double y, double z , string type = "");
    // iterator over pairs (decay constant; contraction coefficient)
    typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
    GaussianIterator firstGaussian() { return _gaussians.begin(); }
    GaussianIterator lastGaussian(){ return _gaussians.end(); }
   
    // adds a Gaussian 
    AOGaussianPrimitive*  addGaussian( double decay, std::vector<double> contraction ) 
    {
        AOGaussianPrimitive* gaussian = new AOGaussianPrimitive(decay, contraction, this);
        _gaussians.push_back( gaussian );
        return gaussian;
    }

    
private:   

    // only class Element can construct shells    
    AOShell( string type, double scale, int numFunc, int startIndex, int offset, vec pos, string atomname, int atomindex, AOBasis* aobasis = NULL ) : _type(type), _scale(scale), _numFunc(numFunc), _startIndex(startIndex), _offset(offset), _pos(pos) , _atomname(atomname), _atomindex(atomindex) { ; }
    
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
    int _offset;
    vec _pos;
    string _atomname;
    int _atomindex;
     
    //AOBasis* _aobasis;
    int detlmax( string shell );
    // vector of pairs of decay constants and contraction coefficients
    vector< AOGaussianPrimitive* > _gaussians;
    
    


};

/*
 * A collection of shells associated with a specific element  
 */
class AOBasis 
{   
public:
       
       // template< class T >
       //static void ReorderMOs(ub::matrix<T> &v, vector<int> const &order )  {
       //static void ReorderMOs(ub::matrix<T> &v, string start, string target )  { 
       void ReorderMOs(ub::matrix<double> &v, string start, string target )  { 
       
          // cout << " Reordering MOs from " << start << " to " << target << endl;
           
          // get reordering vector _start -> target 
          vector<int> order;
          this->getReorderVector( start, target, order);
           
          // Sanity check
          if ( v.size2() != order.size() ) {
              cerr << "Size mismatch in ReorderMOs" << v.size2() << ":" << order.size() << endl;
              throw std::runtime_error( "Abort!");
          }
           
          // actual swapping of coefficients
          for ( unsigned _i_orbital = 0; _i_orbital < v.size1(); _i_orbital++ ){
                for ( unsigned s = 1, d; s < order.size(); ++ s ) {
                    for ( d = order[s]; d < s; d = order[d] ){
                        ;
                    }
                          if ( d == s ) while ( d = order[d], d != s ) swap( v(_i_orbital,s), v(_i_orbital,d) );
                }
          }
           
          // NWChem has some strange minus in d-functions
          if ( start == "nwchem" || target == "nwchem" ){
              
              // get vector with multipliers, e.g. NWChem -> Votca (bloody sign for d_xz)
              vector<int> multiplier;
              this->getMultiplierVector(start, target, multiplier);
              // and reorder rows of _orbitals->_mo_coefficients() accordingly
              this->MultiplyMOs( v , multiplier);
              
          }
           
           
           
       }
       
      //template< class T >   
      void MultiplyMOs(ub::matrix<double> &v, vector<int> const &multiplier )  { 
          // Sanity check
          if ( v.size2() != multiplier.size() ) {
              cerr << "Size mismatch in MultiplyMOs" << v.size2() << ":" << multiplier.size() << endl;
              throw std::runtime_error( "Abort!");
          }

          for ( unsigned _i_orbital = 0; _i_orbital < v.size1(); _i_orbital++ ){

               for ( unsigned _i_basis = 0; _i_basis < v.size2(); _i_basis++ ){
        
                   v(_i_orbital, _i_basis ) = multiplier[_i_basis] * v(_i_orbital, _i_basis );
                   
               }
               
               
           }
       } 

      
      
    // void AOBasisFill( BasisSet* bs , vector<Segment* > segments);
    void AOBasisFill( BasisSet* bs , vector<QMAtom* > segments, int fragbreak = -1);
    void ECPFill( BasisSet* bs , vector<QMAtom* > segments); 
    
    int NumFuncShell( string shell );
    int NumFuncShell_cartesian( string shell );
    int OffsetFuncShell( string shell );
    int OffsetFuncShell_cartesian( string shell );
    int AOBasisSize() {return _AOBasisSize; }
    
    typedef vector< AOShell* >::iterator AOShellIterator;
    AOShellIterator firstShell() { return _aoshells.begin(); }
    AOShellIterator lastShell(){ return _aoshells.end(); }

    // string getType() { return _type; }
    // int    getNumFunc() { return _numFunc ;}
        
    AOShell* getShell( AOShellIterator it ) { return (*it); }
    
    AOShell* getShell( int idx ){ return _aoshells[idx] ;}
    
    AOShell* addShell( string shellType, double shellScale, int shellFunc, int startIndex, int offset, vec pos, string name, int index ) 
    { 
        AOShell* aoshell = new AOShell( shellType, shellScale, shellFunc, startIndex, offset, pos, name, index, this );
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
   int _AOBasisFragA;
   int _AOBasisFragB;
   
   bool _is_stable;
   
    vector<AOShell*> _aoshells;

    // void getReorderVector( string& package, vector<int>& neworder );
    void getReorderVector( string& start, string& target, vector<int>& neworder );
    //void addReorderShell( string& package, string& shell, vector<int>& neworder );
    void addReorderShell( string& start, string& target, string& shell, vector<int>& neworder );
    //void getMultiplierVector( string& package, vector<int>& multiplier );
    void getMultiplierVector( string& start, string& target, vector<int>& multiplier );
    //void addMultiplierShell( string& package, string& shell, vector<int>& multiplier );  
    void addMultiplierShell( string& start, string& target, string& shell, vector<int>& multiplier );  
    
    
    void getTransformationCartToSpherical( string& package, ub::matrix<double>& _trafomatrix );
    void addTrafoCartShell(  AOShell* shell , ub::matrix_range< ub::matrix<double> >& _submatrix );
    
    int getMaxFunctions ( );
};


        
inline void AOBasis::getTransformationCartToSpherical( string& package, ub::matrix<double>& _trafomatrix ){

    if ( package != "gaussian" ){
        cout << " I should not have been called, will do nothing! " << endl;
    } else {
        // go through basisset, determine function sizes
        int _dim_sph = 0;
        int _dim_cart = 0;
        for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
            AOShell* _shell = this->getShell( _is );
            string _type =  _shell->getType();
            
            _dim_sph  += this->NumFuncShell( _type );
            _dim_cart += this->NumFuncShell_cartesian( _type );
            
            
            //this->addTrafoShell( package,  _type, _trafomatrix );
        }   
        
        
        // cout << " Filling trafo matrix of size : " << _dim_sph << " : " << _dim_cart << endl;
        
        // initialize _trafomatrix
        _trafomatrix = ub::zero_matrix<double>( _dim_sph , _dim_cart );
        
        // now fill it
        int _row_start = 0;
        int _col_start = 0;
        for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
            AOShell* _shell = this->getShell( _is );
            string _type =  _shell->getType();
            int _row_end = _row_start + this->NumFuncShell( _type );
            int _col_end = _col_start + this->NumFuncShell_cartesian( _type );
            
            ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange( _trafomatrix, _row_start, _row_end, _col_start, _col_end);
            this->addTrafoCartShell(  _shell, _submatrix  );
            
            _row_start = _row_end;
            _col_start = _col_end;
            
        } 
        
        
    }

}


inline int AOShell::detlmax( string shell_type ) {
    int _lmax;
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _lmax = 0;}
       if ( shell_type == "P" ){ _lmax = 1;}
       if ( shell_type == "D" ){ _lmax = 2;}
       if ( shell_type == "F" ){ _lmax = 3;}
       if ( shell_type == "G" ){ _lmax = 4;}
    } else {
        // for combined shells check all contributions
        _lmax = 0;
        for( unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell =    string( shell_type, i, 1 );
            int _test = this->detlmax( local_shell  );
            if ( _test > _lmax ) { _lmax = _test;} 
        }
    }

    return _lmax;
}





inline void AOBasis::addTrafoCartShell( AOShell* shell , ub::matrix_range< ub::matrix<double> >& _trafo ){
    
    // cout << "getting trafo of shell type :" << shell->getType() << endl;
    // fill _local according to _lmax;
    int _lmax = shell->getLmax();
    string _type = shell->getType();
    
    int _sph_size = this->NumFuncShell( _type ) + this->OffsetFuncShell( _type );
    int _cart_size = this->NumFuncShell_cartesian( _type ) + this->OffsetFuncShell_cartesian( _type )  ;
    
    // cout << "    local size : " << _sph_size << " : " << _cart_size << endl;
    
    ub::matrix<double> _local =  ub::zero_matrix<double>(_sph_size,_cart_size);

    // s-functions
    _local(0,0) = 1.0; // s
    
    // p-functions
    if ( _lmax > 0 ){
        _local(1,1) = 1.0; 
        _local(2,2) = 1.0;
        _local(3,3) = 1.0;
    }

    // d-functions
    if ( _lmax > 1 ){
        _local(4,4) = -0.5;             // d3z2-r2 (dxx)
        _local(4,5) = -0.5;             // d3z2-r2 (dyy)
        _local(4,6) =  1.0;             // d3z2-r2 (dzz)
        _local(5,8) =  1.0;             // dxz
        _local(6,9) =  1.0;             // dyz
        _local(7,4) = 0.5*sqrt(3.0);    // dx2-y2 (dxx)
        _local(7,5) = -_local(7,4);      // dx2-y2 (dyy)
        _local(8,7) = 1.0;              // dxy
     }
  
    if ( _lmax > 2 ){
        cerr << " Gaussian input with f- functions or higher not yet supported!" << endl;
        exit(1);
    }

    // now copy to _trafo
    for ( int _i_sph = 0 ; _i_sph < this->NumFuncShell( _type ) ; _i_sph++ ){
        for  ( int _i_cart = 0 ; _i_cart < this->NumFuncShell_cartesian( _type ) ; _i_cart++ ){
            
            
            _trafo( _i_sph , _i_cart ) = _local( _i_sph + this->OffsetFuncShell( _type ) , _i_cart +  this->OffsetFuncShell_cartesian( _type ) );
            
        }
    }

    
}






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


inline void AOBasis::getMultiplierVector( string& start, string& target, vector<int>& multiplier){

    // go through basisset
    for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
        AOShell* _shell = this->getShell( _is );
        string _type =  _shell->getType();
        this->addMultiplierShell(  start, target, _type, multiplier );
    }
    
        }

        inline void AOBasis::addMultiplierShell(string& start, string& target, string& shell_type, vector<int>& multiplier) {


            if (target == "votca") {
                // current length of vector
                //int _cur_pos = multiplier.size() - 1;

                // single type shells defined here
                if (shell_type.length() == 1) {
                    if (shell_type == "S") {
                        multiplier.push_back(1);
                    }

                    if (shell_type == "P") {
                        multiplier.push_back(1);
                        multiplier.push_back(1);
                        multiplier.push_back(1);
                    }
                    if (shell_type == "D") {
                        if (start == "nwchem") {
                            multiplier.push_back(-1);
                            multiplier.push_back(1);
                            multiplier.push_back(1);
                            multiplier.push_back(1);
                            multiplier.push_back(1);
                        } else {
                            cerr << "Tried to get multipliers d-functions from package " << start << ".";
                            throw std::runtime_error("Multiplication not implemented yet!");
                        }
                    }
                    if (shell_type == "F") {
                        cerr << "Tried to get multipliers for f-functions . ";
                        throw std::runtime_error("Multiplication not implemented yet!");
                    }
                    if (shell_type == "G") {
                        cerr << "Tried to get multipliers g-functions . ";
                        throw std::runtime_error("Multiplication not implemented yet!");
                    }
                } else {
                    // for combined shells, iterate over all contributions
                    //_nbf = 0;
                    for (unsigned i = 0; i < shell_type.length(); ++i) {
                        string local_shell = string(shell_type, i, 1);
                        this->addMultiplierShell(start, target, local_shell, multiplier);
                    }
                }
            } else {

                cerr << "Tried to reorder functions (multiplier)  from " << start << " to " << target << endl;
                throw std::runtime_error("Reordering not implemented yet!");


            }

        }


inline void AOBasis::getReorderVector( string& start, string& target, vector<int>& neworder){

    // go through basisset
    for (vector< AOShell* >::iterator _is = this->firstShell(); _is != this->lastShell() ; _is++ ) {
        AOShell* _shell = this->getShell( _is );
        string _type =  _shell->getType();
        this->addReorderShell( start, target, _type, neworder );
    }
    
}


inline void AOBasis::addReorderShell( string& start, string& target,  string& shell_type, vector<int>& neworder ) {
    
    // current length of vector
    int _cur_pos = neworder.size() -1 ;
    
    if ( target == "votca" ){
    
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
           if ( start == "gaussian"){
               neworder.push_back( _cur_pos + 4 );
               neworder.push_back( _cur_pos + 1 );
               neworder.push_back( _cur_pos + 2 );
               neworder.push_back( _cur_pos + 5 );
               neworder.push_back( _cur_pos + 3 );
           } else if ( start == "nwchem") {
               neworder.push_back( _cur_pos + 3  ); 
               neworder.push_back( _cur_pos + 2 );
               neworder.push_back( _cur_pos + 4 );
               //neworder.push_back( -(_cur_pos + 1) ); // bloody inverted sign // BUG!!!!!!!
               neworder.push_back( _cur_pos + 1 ); 
               neworder.push_back( _cur_pos + 5 );               
           } else {
               cerr << "Tried to reorder d-functions from package " << start << ".";
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
        for( unsigned i = 0; i < shell_type.length(); ++i) {
           string local_shell =    string( shell_type, i, 1 );
           this->addReorderShell( start, target, local_shell, neworder  );
        }
    }
    } else {
        
        cerr << "Tried to reorder functions (neworder) from " << start << " to " << target << endl;
        throw std::runtime_error( "Reordering not implemented yet!");
        
    } 
    
}




inline void AOBasis::AOBasisFill(BasisSet* bs , vector<QMAtom* > _atoms, int _fragbreak  ) {
    
        vector< QMAtom* > :: iterator ait;
        std::vector < QMAtom* > :: iterator atom;

       _AOBasisSize = 0;
       _is_stable = true; // _is_stable = true corresponds to gwa_basis%S_ev_stable = .false. 
       
       int _atomidx = 0;
       
       // loop over atoms
       for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
          // get coordinates of this atom and convert from Angstrom to Bohr
          vec pos;
          pos.setX( (*ait)->x * 1.8897259886  );
          pos.setY( (*ait)->y * 1.8897259886  );
          pos.setZ( (*ait)->z * 1.8897259886  );
          // get element type of the atom
          string  name = (*ait)->type;
          // get the basis set entry for this element
          Element* element = bs->getElement(name);
                    // and loop over all shells
          for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
               Shell* shell = (*its);
               // we don't like contracted basis sets yet
               //if ( shell->getSize() > 1 ) {
               //    cerr << "We have a contracted basis set!" << flush;
               //} else {
                   AOShell* aoshell = addShell( shell->getType(), shell->getScale(), NumFuncShell( shell->getType() ), _AOBasisSize, OffsetFuncShell( shell->getType() ), pos, name, _atomidx );
                   _AOBasisSize += NumFuncShell( shell->getType() );
                   for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                      GaussianPrimitive* gaussian = *itg;
                      aoshell->addGaussian(gaussian->decay, gaussian->contraction);
                //   }
               }
          }
          
          if ( _atomidx < _fragbreak ) _AOBasisFragA = _AOBasisSize;
          
          _atomidx++;
      }
       
       if ( _fragbreak < 0 ) {
           _AOBasisFragA = _AOBasisSize;
           _AOBasisFragB = 0;
       } else {
           _AOBasisFragB = _AOBasisSize - _AOBasisFragA;
       }
           
    
}



inline void AOBasis::ECPFill(BasisSet* bs , vector<QMAtom* > _atoms  ) {
    
        vector< QMAtom* > :: iterator ait;
        std::vector < QMAtom* > :: iterator atom;

       _AOBasisSize = 0;
       _is_stable = true; // _is_stable = true corresponds to gwa_basis%S_ev_stable = .false. 
       
       int _atomidx = 0;
       
       // loop over atoms
       for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
          // get coordinates of this atom and convert from Angstrom to Bohr
          vec pos;
          pos.setX( (*ait)->x * 1.8897259886  );
          pos.setY( (*ait)->y * 1.8897259886  );
          pos.setZ( (*ait)->z * 1.8897259886  );
          // get element type of the atom
          string  name = (*ait)->type;
          // get the basis set entry for this element
          Element* element = bs->getElement(name);
          // cout << " Name " << name << endl;
          // and loop over all shells
          
          int lmax;
          for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
               Shell* shell = (*its);
               //cout << " Shell " << shell->getType() << endl;
               // we don't like contracted basis sets yet
               //if ( shell->getSize() > 1 ) {
               //    cerr << "We have a contracted basis set!" << flush;
               //} else {
               string local_shell =    string( shell->getType(), 0, 1 );
               int l;
               if ( local_shell == "S" ) l =0;
               if ( local_shell == "P" ) l =1;
               if ( local_shell == "D" ) l =2;
               if ( local_shell == "F" ) l =3;
               if (its == element->firstShell()) lmax = l;
               // first shell is local component, identification my negative angular momentum
               
               
                   AOShell* aoshell = addShell( shell->getType(), shell->getScale(), lmax, l, l, pos, name, _atomidx );
                   _AOBasisSize += NumFuncShell( shell->getType() );
                   for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                      GaussianPrimitive* gaussian = *itg;
                      aoshell->addGaussian( gaussian->decay, gaussian->contraction);
                //   }
               }
          }
          
        
          
          _atomidx++;
      }
       
           
    
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
        for( unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell =    string( shell_type, i, 1 );
            _nbf += this->NumFuncShell( local_shell  );
        }
    }

    return _nbf;
}
    
    
 inline int AOBasis::NumFuncShell_cartesian( string shell_type ) {
    int _nbf;
    // single type shells defined here
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 1;}
       if ( shell_type == "P" ){ _nbf = 3;}
       if ( shell_type == "D" ){ _nbf = 6;}
       if ( shell_type == "F" ){ _nbf = 10;}
       if ( shell_type == "G" ){ _nbf = 15;}
    } else {
        // for combined shells, sum over all contributions
        _nbf = 0;
        for( unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell =    string( shell_type, i, 1 );
            _nbf += this->NumFuncShell_cartesian( local_shell  );
        }
    }

    return _nbf;
}
    
    
    
       

    inline int AOBasis::OffsetFuncShell_cartesian( string shell_type ) {
    int _nbf;
    // single type shells
    if ( shell_type.length() == 1 ){
       if ( shell_type == "S" ){ _nbf = 0;}
       if ( shell_type == "P" ){ _nbf = 1;}
       if ( shell_type == "D" ){ _nbf = 4;}
       if ( shell_type == "F" ){ _nbf = 10;}
       if ( shell_type == "G" ){ _nbf = 20;}
    } else {
        // for combined shells, go over all contributions and find minimal offset
        _nbf = 1000;
        for(unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell = string( shell_type, i, 1 );
            int _test = this->OffsetFuncShell_cartesian( local_shell  );
            if ( _test < _nbf ) { _nbf = _test;} 
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
        for(unsigned i = 0; i < shell_type.length(); ++i) {
            string local_shell = string( shell_type, i, 1 );
            int _test = this->OffsetFuncShell( local_shell  );
            if ( _test < _nbf ) { _nbf = _test;} 
        }   
    }
    return _nbf;
        }

   
    
       inline void AOShell::EvalAOGradspace(ub::matrix_range<ub::matrix<double> >& gradAOvalues, double x, double y, double z, string type ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();
            
            
  typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

            double _expofactor = pow(2.0*alpha/pi,0.75) * exp( -alpha * distsq );
            
           // std::vector< std::vector<double> > _AOGradevaluated(3);
         
            // split combined shells
            int _i_func = -1;
            int i_act;
            for (unsigned i = 0; i < shell_type.length(); ++i) {
                string single_shell = string(shell_type, i, 1);
                // single type shells
                if ( single_shell == "S") {
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) += _contractions[0] * -2.0*alpha*center_x*_expofactor; // x gradient of s-function
                    gradAOvalues(i_act,1) += _contractions[0] * -2.0*alpha*center_y*_expofactor; // y gradient of s-function
                    gradAOvalues(i_act,2) += _contractions[0] * -2.0*alpha*center_z*_expofactor; // z gradient of s-function
                    _i_func = i_act;
                }
                if ( single_shell == "P") {
                    
                    // px -functions
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) += _contractions[1] * 2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_x*center_x)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // z gradient 
                    
                    // py -functions
                    i_act = _i_func +2;
                    gradAOvalues(i_act,0) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[1] * 2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_y*center_y) *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z*_expofactor; // z gradient 
                    
                     // pz -functions
                    i_act = _i_func +3;
                    gradAOvalues(i_act,0) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[1] * -2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[1] * 2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_z*center_z)*_expofactor; // z gradient                    
                    _i_func = i_act;

                }
                if ( single_shell == "D") {
                             
                    // dxz function
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) += _contractions[2] * 4.0*alpha * (center_z - 2.0*alpha*center_x*center_x*center_z)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[2] * -8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[2] * 4.0*alpha * (center_x - 2.0*alpha*center_x*center_z*center_z)*_expofactor; // z gradient 

                    // dyz function
                    i_act = _i_func+2;                    
                    gradAOvalues(i_act,0) += _contractions[2] * -8.0*alpha*alpha * center_x*center_y*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[2] * 4.0*alpha * (center_z - 2.0*alpha*center_y*center_y*center_z)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) += _contractions[2] * 4.0*alpha * (center_y - 2.0*alpha*center_y*center_z*center_z)*_expofactor; // z gradient 

                    // dxy function
                    i_act = _i_func+3;                    
                    gradAOvalues(i_act,0) += _contractions[2] * 4.0*alpha * (center_y - 2.0*alpha*center_x*center_x*center_y)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) += _contractions[2] * 4.0*alpha * (center_x - 2.0*alpha*center_x*center_y*center_y)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) +=_contractions[2] * -8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // z gradient 

                    // d3z2-r2-function
                    i_act = _i_func+4;                    
                    gradAOvalues(i_act,0) +=_contractions[2] * -4.0*alpha/sqrt(3.0) * center_x * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,1) +=_contractions[2] * -4.0*alpha/sqrt(3.0) * center_y * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,2) += _contractions[2] * 4.0*alpha/sqrt(3.0) * center_z * (2.0 - alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    
                    // dx2-y2-function
                    i_act = _i_func+5;                    
                    gradAOvalues(i_act,0) += _contractions[2] * 4.0*alpha * center_x * (1.0 - alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,1) +=_contractions[2] * -4.0*alpha * center_y * (1.0 + alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,2) += _contractions[2] * -4.0*alpha*alpha * center_z * (center_x*center_x - center_y*center_y) * _expofactor;


                }
                if ( single_shell == "F") {
                    cerr << " F functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
                if ( single_shell == "G") {
                    cerr << " G functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
            }
            }// contractions


        }
    /*
       inline void AOShell::EvalAOGradspace(ub::matrix<double>& gradAOvalues, double x, double y, double z, string type ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            double alpha = this->_gaussians[0]->decay; // only uncontracted for testing
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            double _expofactor = pow(2.0*alpha/pi,0.75) * exp( -alpha * distsq );
            
            std::vector< std::vector<double> > _AOGradevaluated(3);
         
            // split combined shells
            int _i_func = -1;
            int i_act;
            for (int i = 0; i < shell_type.length(); ++i) {
                string single_shell = string(shell_type, i, 1);
                // single type shells
                if ( single_shell == "S") {
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) =-2.0*alpha*center_x*_expofactor; // x gradient of s-function
                    gradAOvalues(i_act,1) =-2.0*alpha*center_y*_expofactor; // y gradient of s-function
                    gradAOvalues(i_act,2) =-2.0*alpha*center_z*_expofactor; // z gradient of s-function
                    _i_func = i_act;
                }
                if ( single_shell == "P") {
                    
                    // px -functions
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) =2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_x*center_x)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // z gradient 
                    
                    // py -functions
                    i_act = _i_func +2;
                    gradAOvalues(i_act,0) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_y*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_y*center_y) *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =-2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z*_expofactor; // z gradient 
                    
                     // pz -functions
                    i_act = _i_func +3;
                    gradAOvalues(i_act,0) =-2.0*sqrt(alpha) * 2.0*alpha*center_x*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =-2.0*sqrt(alpha) * 2.0*alpha*center_y*center_z *_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =2.0*sqrt(alpha) * (1.0 - 2.0*alpha*center_z*center_z)*_expofactor; // z gradient                    
                    _i_func = i_act;

                }
                if ( single_shell == "D") {
                             // dxz, dyz, dxy, d3z2-r2, dx2-y2
                    // dxz function
                    i_act = _i_func+1;
                    gradAOvalues(i_act,0) = 4.0*alpha * (center_z - 2.0*alpha*center_x*center_x*center_z)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) =-8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) = 4.0*alpha * (center_x - 2.0*alpha*center_x*center_z*center_z)*_expofactor; // z gradient 

                    // dyz function
                    i_act = _i_func+2;                    
                    gradAOvalues(i_act,0) =-8.0*alpha*alpha * center_x*center_y*center_z*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) = 4.0*alpha * (center_z - 2.0*alpha*center_y*center_y*center_z)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) = 4.0*alpha * (center_y - 2.0*alpha*center_y*center_z*center_z)*_expofactor; // z gradient 

                    // dxy function
                    i_act = _i_func+3;                    
                    gradAOvalues(i_act,0) = 4.0*alpha * (center_y - 2.0*alpha*center_x*center_x*center_y)*_expofactor; // x gradient 
                    gradAOvalues(i_act,1) = 4.0*alpha * (center_x - 2.0*alpha*center_x*center_y*center_y)*_expofactor; // y gradient 
                    gradAOvalues(i_act,2) =-8.0*alpha*alpha *center_x*center_y*center_z*_expofactor; // z gradient 

                    // d3z2-r2-function
                    i_act = _i_func+4;                    
                    gradAOvalues(i_act,0) =-4.0*alpha/sqrt(3.0) * center_x * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,1) =-4.0*alpha/sqrt(3.0) * center_y * (1.0 + alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    gradAOvalues(i_act,2) = 4.0*alpha/sqrt(3.0) * center_z * (2.0 - alpha *(3.0*center_z*center_z - distsq) ) * _expofactor;
                    
                    // dx2-y2-function
                    i_act = _i_func+5;                    
                    gradAOvalues(i_act,0) = 4.0*alpha * center_x * (1.0 - alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,1) =-4.0*alpha * center_y * (1.0 + alpha*(center_x*center_x - center_y*center_y)) * _expofactor;
                    gradAOvalues(i_act,2) =-4.0*alpha*alpha * center_z * (center_x*center_x - center_y*center_y) * _expofactor;


                }
                if ( single_shell == "F") {
                    cerr << " F functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
                if ( single_shell == "G") {
                    cerr << " G functions not implemented in AOGradeval at the moment!" << endl;
                    exit(1);
                }
            }



        } */
    
        
       
       inline void AOShell::EvalAOIntegral(ub::matrix_range<ub::matrix<double> >& AOvalues){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
      //      double center_x = x - this->_pos.getX();
    //        double center_y = y - this->_pos.getY();
  //          double center_z = z - this->_pos.getZ();
            // need decay constant
            
//            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            
            
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

                double _factor = pow(2.0 * pi / alpha, 0.75) ;

                // split combined shells
                int _i_func = -1;
                //int i_act;
                for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string single_shell = string(shell_type, i, 1);
                    // single type shells
                    if (single_shell == "S") {
                        AOvalues(0, _i_func + 1) += _contractions[0] * _factor; // s-function
                        _i_func++;
                    }
                }
            } // contractions

        }
       
       
       
           inline void AOShell::EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues, ub::matrix_range<ub::matrix<double> >& gradAOvalues, double x, double y, double z ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            
            
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

                double _expofactor = pow(2.0 * alpha / pi, 0.75) * exp(-alpha * distsq);

                // split combined shells
                int _i_func = -1;
                int i_act;
                for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string single_shell = string(shell_type, i, 1);
                    // single type shells
                    if (single_shell == "S") {
                        AOvalues(0, _i_func + 1) += _contractions[0] * _expofactor; // s-function

                            i_act = _i_func+1;
                            gradAOvalues(0, i_act) += _contractions[0] * -2.0 * alpha * center_x*_expofactor; // x gradient of s-function
                            gradAOvalues(1, i_act) += _contractions[0] * -2.0 * alpha * center_y*_expofactor; // y gradient of s-function
                            gradAOvalues(2, i_act) += _contractions[0] * -2.0 * alpha * center_z*_expofactor; // z gradient of s-function
                        
                        _i_func++;
                    }
                    if (single_shell == "P") {
                        AOvalues(0, _i_func + 1) += _contractions[1] * 2.0 * sqrt(alpha) * center_x*_expofactor; // px-function
                        AOvalues(0, _i_func + 2) += _contractions[1] * 2.0 * sqrt(alpha) * center_y*_expofactor; // py-function
                        AOvalues(0, _i_func + 3) += _contractions[1] * 2.0 * sqrt(alpha) * center_z*_expofactor; // pz-function
                            i_act = _i_func+1;
                            gradAOvalues(0, i_act) += _contractions[1] * 2.0 * sqrt(alpha) * (1.0 - 2.0 * alpha * center_x * center_x) * _expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_y*_expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_z*_expofactor; // z gradient 

                            // py -functions
                            i_act = _i_func + 2;
                            gradAOvalues(0, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_y*_expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[1] * 2.0 * sqrt(alpha) * (1.0 - 2.0 * alpha * center_y * center_y) * _expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_y * center_z*_expofactor; // z gradient 

                            // pz -functions
                            i_act = _i_func + 3;
                            gradAOvalues(0, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_x * center_z*_expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[1] * -2.0 * sqrt(alpha) * 2.0 * alpha * center_y * center_z *_expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[1] * 2.0 * sqrt(alpha) * (1.0 - 2.0 * alpha * center_z * center_z) * _expofactor; // z gradient      



                        _i_func += 3;
                    }
                    if (single_shell == "D") {
                        // dxz, dyz, dxy, d3z2-r2, dx2-y2
                        AOvalues(0, _i_func + 1) += _contractions[2] * 4.0 * alpha * center_x * center_z*_expofactor; // dxz-function
                        AOvalues(0, _i_func + 2) += _contractions[2] * 4.0 * alpha * center_y * center_z*_expofactor; // dyz-function
                        AOvalues(0, _i_func + 3) += _contractions[2] * 4.0 * alpha * center_x * center_y*_expofactor; // dxy-function
                        AOvalues(0, _i_func + 4) += _contractions[2] * 2.0 * alpha / sqrt(3.0)*(3.0 * center_z * center_z - distsq) * _expofactor; // d3z2-r2-function
                        AOvalues(0, _i_func + 5) += _contractions[2] * 2.0 * alpha * (center_x * center_x - center_y * center_y) * _expofactor; // dx2-y2-function


                            i_act = _i_func+1;
                            gradAOvalues(0, i_act) += _contractions[2] * 4.0 * alpha * (center_z - 2.0 * alpha * center_x * center_x * center_z) * _expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[2] * -8.0 * alpha * alpha * center_x * center_y * center_z*_expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[2] * 4.0 * alpha * (center_x - 2.0 * alpha * center_x * center_z * center_z) * _expofactor; // z gradient 

                            // dyz function
                            i_act = _i_func + 2;
                            gradAOvalues(0, i_act) += _contractions[2] * -8.0 * alpha * alpha * center_x * center_y * center_z*_expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[2] * 4.0 * alpha * (center_z - 2.0 * alpha * center_y * center_y * center_z) * _expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[2] * 4.0 * alpha * (center_y - 2.0 * alpha * center_y * center_z * center_z) * _expofactor; // z gradient 

                            // dxy function
                            i_act = _i_func + 3;
                            gradAOvalues(0, i_act) += _contractions[2] * 4.0 * alpha * (center_y - 2.0 * alpha * center_x * center_x * center_y) * _expofactor; // x gradient 
                            gradAOvalues(1, i_act) += _contractions[2] * 4.0 * alpha * (center_x - 2.0 * alpha * center_x * center_y * center_y) * _expofactor; // y gradient 
                            gradAOvalues(2, i_act) += _contractions[2] * -8.0 * alpha * alpha * center_x * center_y * center_z*_expofactor; // z gradient 

                            // d3z2-r2-function
                            i_act = _i_func + 4;
                            gradAOvalues(0, i_act) += _contractions[2] * -4.0 * alpha / sqrt(3.0) * center_x * (1.0 + alpha * (3.0 * center_z * center_z - distsq)) * _expofactor;
                            gradAOvalues(1, i_act) += _contractions[2] * -4.0 * alpha / sqrt(3.0) * center_y * (1.0 + alpha * (3.0 * center_z * center_z - distsq)) * _expofactor;
                            gradAOvalues(2, i_act) += _contractions[2] * 4.0 * alpha / sqrt(3.0) * center_z * (2.0 - alpha * (3.0 * center_z * center_z - distsq)) * _expofactor;

                            // dx2-y2-function
                            i_act = _i_func + 5;
                            gradAOvalues(0, i_act) += _contractions[2] * 4.0 * alpha * center_x * (1.0 - alpha * (center_x * center_x - center_y * center_y)) * _expofactor;
                            gradAOvalues(1, i_act) += _contractions[2] * -4.0 * alpha * center_y * (1.0 + alpha * (center_x * center_x - center_y * center_y)) * _expofactor;
                            gradAOvalues(2, i_act) += _contractions[2] * -4.0 * alpha * alpha * center_z * (center_x * center_x - center_y * center_y) * _expofactor;



                        
                        
                        
                        _i_func += 5;
                    }
                    if (single_shell == "F") {
                        cerr << " F functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                    if (single_shell == "G") {
                        cerr << " G functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                }
            } // contractions

        }
           
           
           
            inline void AOShell::EvalAOspace(ub::matrix_range<ub::matrix<double> >& AOvalues, double x, double y, double z ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();

            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
            // iterate over Gaussians in this shell
            for (GaussianIterator itr = firstGaussian(); itr != lastGaussian(); ++itr) {

                const double& alpha = (*itr)->decay;
                std::vector<double> _contractions = (*itr)->contraction;

                double _expofactor = pow(2.0 * alpha / pi, 0.75) * exp(-alpha * distsq);

                // split combined shells
                int _i_func = -1;
                
                for (unsigned i = 0; i < shell_type.length(); ++i) {
                    string single_shell = string(shell_type, i, 1);
                    // single type shells
                    if (single_shell == "S") {
                        AOvalues(0, _i_func + 1) += _contractions[0] * _expofactor; // s-function        
                        _i_func++;
                    }
                    if (single_shell == "P") {
                        AOvalues(0, _i_func + 1) += _contractions[1] * 2.0 * sqrt(alpha) * center_x*_expofactor; // px-function
                        AOvalues(0, _i_func + 2) += _contractions[1] * 2.0 * sqrt(alpha) * center_y*_expofactor; // py-function
                        AOvalues(0, _i_func + 3) += _contractions[1] * 2.0 * sqrt(alpha) * center_z*_expofactor; // pz-function

                        _i_func += 3;
                    }
                    if (single_shell == "D") {
                        // dxz, dyz, dxy, d3z2-r2, dx2-y2
                        AOvalues(0, _i_func + 1) += _contractions[2] * 4.0 * alpha * center_x * center_z*_expofactor; // dxz-function
                        AOvalues(0, _i_func + 2) += _contractions[2] * 4.0 * alpha * center_y * center_z*_expofactor; // dyz-function
                        AOvalues(0, _i_func + 3) += _contractions[2] * 4.0 * alpha * center_x * center_y*_expofactor; // dxy-function
                        AOvalues(0, _i_func + 4) += _contractions[2] * 2.0 * alpha / sqrt(3.0)*(3.0 * center_z * center_z - distsq) * _expofactor; // d3z2-r2-function
                        AOvalues(0, _i_func + 5) += _contractions[2] * 2.0 * alpha * (center_x * center_x - center_y * center_y) * _expofactor; // dx2-y2-function                             
                        _i_func += 5;
                    }
                    if (single_shell == "F") {
                        cerr << " F functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                    if (single_shell == "G") {
                        cerr << " G functions not implemented in AOeval at the moment!" << endl;
                        exit(1);
                    }
                }
            } // contractions

        }
        
           
           
           
        
           
        
               
  /*         inline void AOShell::EvalAOspace(ub::matrix<double>& AOvalues, double x, double y, double z, string type ){

            // need type of shell
            string shell_type = this->_type;
            // need position of shell
            double center_x = x - this->_pos.getX();
            double center_y = y - this->_pos.getY();
            double center_z = z - this->_pos.getZ();
            // need decay constant
            double alpha = this->_gaussians[0]->decay; // only uncontracted for testing
            double distsq = center_x*center_x + center_y*center_y +  center_z*center_z;
            const double pi = boost::math::constants::pi<double>();


            double _expofactor = pow(2.0*alpha/pi,0.75) * exp( -alpha * distsq );
            
            // split combined shells
            int _i_func = -1;
            for (int i = 0; i < shell_type.length(); ++i) {
                string single_shell = string(shell_type, i, 1);
                // single type shells
                if ( single_shell == "S") {
                    AOvalues(_i_func +1 ,0) = _expofactor; // s-function
                    _i_func++;
                }
                if ( single_shell == "P") {
                    AOvalues(_i_func +1,0) = 2.0*sqrt(alpha)*center_x*_expofactor; // px-function
                    AOvalues(_i_func +2,0) = 2.0*sqrt(alpha)*center_y*_expofactor; // py-function
                    AOvalues(_i_func +3,0) = 2.0*sqrt(alpha)*center_z*_expofactor; // pz-function
                    _i_func += 3;
                }
                if ( single_shell == "D") {
                    // dxz, dyz, dxy, d3z2-r2, dx2-y2
                    AOvalues(_i_func +1,0) = 4.0*alpha*center_x*center_z*_expofactor; // dxz-function
                    AOvalues(_i_func +2,0) = 4.0*alpha*center_y*center_z*_expofactor; // dyz-function
                    AOvalues(_i_func +3,0) = 4.0*alpha*center_x*center_y*_expofactor; // dxy-function
                    AOvalues(_i_func +4,0) = 2.0*alpha/sqrt(3.0)*(3.0*center_z*center_z - distsq)*_expofactor; // d3z2-r2-function
                    AOvalues(_i_func +5,0) = 2.0*alpha*(center_x*center_x - center_y*center_y)*_expofactor; // dx2-y2-function
                    _i_func += 5;
                }
                if ( single_shell == "F") {
                    cerr << " F functions not implemented in AOeval at the moment!" << endl;
                    exit(1);
                }
                if ( single_shell == "G") {
                    cerr << " G functions not implemented in AOeval at the moment!" << endl;
                    exit(1);
                }
            }   


        } */
        
        
    
    
}}

#endif	/* AOBASIS_H */

