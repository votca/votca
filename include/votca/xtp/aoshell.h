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

#ifndef __XTP_AOSHELL__H
#define	__XTP_AOSHELL__H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/basisset.h>


using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {
namespace ub = boost::numeric::ublas;
class AOBasis;
class AOShell;  




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

    
}}

#endif	/* AOSHELL_H */

