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

#ifndef __XTP_AOSHELL__H
#define	__XTP_AOSHELL__H


#include <boost/math/constants/constants.hpp>
#include <votca/xtp/eigen.h>
#include <votca/tools/constants.h>

#include <votca/xtp/basisset.h>
#include "qmatom.h"

namespace votca { namespace xtp {

class AOBasis;
class AOShell;  


// Gaussian function: contraction*exp(-decay*r^2)
class AOGaussianPrimitive 
{
    friend class AOShell;
public:
    
    

    double getPowfactor()const {return powfactor;}
    int    getPower()const{return power;}
    double getDecay()const {return decay;}
    const std::vector<double>& getContraction()const {return contraction;}
    const AOShell* getShell() const{return aoshell;}
private:
     
    int power; // used in pseudopotenials only
    double decay;
    std::vector<double> contraction;
    AOShell* aoshell;
    double powfactor;//used in evalspace to speed up DFT
    // private constructor, only a shell can create a primitive
    AOGaussianPrimitive( double _decay, std::vector<double> _contraction, AOShell *_aoshell = NULL ) 
    : power(-1),decay(_decay),
            contraction(_contraction),
            aoshell(_aoshell) {powfactor=pow(2.0 * decay / boost::math::constants::pi<double>(), 0.75) ; }

    AOGaussianPrimitive( int _power, double _decay, std::vector<double> _contraction, AOShell *_aoshell = NULL ) 
    : power(_power),
    decay(_decay),
    contraction(_contraction),
    aoshell(_aoshell) {powfactor=pow(2.0 * decay / boost::math::constants::pi<double>(), 0.75) ; }
};      
    
/*
 * S, P, or D functions in a Gaussian-basis expansion
 */
class AOShell 
{
    //friend class AOElement;
    friend class AOBasis;
public:

    const std::string& getType() const{ return _type; }
    int    getNumFunc() const{ return _numFunc ;}
    int    getStartIndex() const{ return _startIndex ;}
    int    getOffset() const{ return _offset ;}
    int    getIndex() const{ return _atomindex;}
    const std::string& getName() const{ return _atomname;}
    
    int getLmax(  ) const{ return _Lmax;}
    int getLmin(  ) const{ return _Lmin;}
    
    bool isNonLocal(  ) const{ return _nonlocal;}
    
    const tools::vec& getPos() const{ return _pos; }
    double getScale() const{ return _scale; }
    
    int getSize() const{ return _gaussians.size(); }
    
    void CalcMinDecay(){
     _mindecay=std::numeric_limits<double>::max();
     for(auto& gaussian:_gaussians){
         if(gaussian.getDecay()<_mindecay){
             _mindecay=gaussian.getDecay();
         }
     }
     return;
    }
    
    double getMinDecay() const{return _mindecay;}
    
    
  void EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>&  AOvalues, const tools::vec& grid_pos ) const;
  void EvalAOspace(Eigen::VectorBlock<Eigen::VectorXd>&  AOvalues,Eigen::Block< Eigen::MatrixX3d >& AODervalues, const tools::vec& grid_pos ) const;

    // iterator over pairs (decay constant; contraction coefficient)
    typedef std::vector< AOGaussianPrimitive >::const_iterator GaussianIterator;
    GaussianIterator begin() const{ return _gaussians.begin(); }
    GaussianIterator end()const{ return _gaussians.end(); }
   
    // adds a Gaussian 
    void  addGaussian( double decay, std::vector<double> contraction ) 
    {
        AOGaussianPrimitive gaussian = AOGaussianPrimitive(decay, contraction, this);
        _gaussians.push_back( gaussian );
        return;
    }
    // used for ecps
    void  addGaussian( int power, double decay, std::vector<double> contraction ) 
    {                                                                                
        AOGaussianPrimitive gaussian = AOGaussianPrimitive(power, decay, contraction, this); 
        _gaussians.push_back( gaussian );                                                  
        return;                                                                 
    }                                                                                

    void normalizeContraction();
    
private:   

    // only class aobasis can construct shells    
    AOShell( const Shell& shell, const QMAtom & atom, int startIndex)
            : _type(shell.getType()),_Lmax(shell.getLmax()),_Lmin(shell.getLmin()),
                    _scale(shell.getScale()), _numFunc(shell.getnumofFunc()),
                    _startIndex(startIndex), _offset(shell.getOffset()), _pos(atom.getPos()) , 
                    _atomname(atom.getType()), _atomindex(atom.getAtomID()) { ; }
    // for ECPs
    AOShell( const Shell& shell, const QMAtom & atom, int startIndex, bool nonlocal)
            : _type(shell.getType()),_Lmax(shell.getLmax()),_Lmin(shell.getLmin()),
                    _scale(shell.getScale()), _numFunc(shell.getnumofFunc()),
                    _startIndex(startIndex), _offset(shell.getOffset()), _pos(atom.getPos()) , 
                    _atomname(atom.getType()), _atomindex(atom.getAtomID()),_nonlocal(nonlocal) { ; }
            
    
    // only class aobasis can destruct shells
    ~AOShell(){};
    
    // shell type (S, P, D))
    std::string _type;
    int _Lmax;
    int _Lmin;
    // scaling factor
    double _scale;
    // number of functions in shell
    int _numFunc;
    double _mindecay;
    int _startIndex;
    int _offset;
    tools::vec _pos;
    std::string _atomname;
    int _atomindex;
    
    //used for ecp calculations
    bool _nonlocal;
     
    // vector of pairs of decay constants and contraction coefficients
    std::vector< AOGaussianPrimitive > _gaussians;
    
};

    
}}

#endif	/* AOSHELL_H */

