/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef __XTP_AOMATRIX__H
#define	__XTP_AOMATRIX__H

#include <votca/xtp/aobasis.h>
#include <votca/xtp/aoshell.h>
#include <votca/ctp/apolarsite.h>
#include <votca/ctp/polarseg.h>
#include <votca/xtp/votca_config.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#define BOOST_DISABLE_ASSERTS //could be used to slighlty speed up calculation but the compile time simply goes boom
#include <boost/multi_array.hpp>
#include "basisset.h"
//#include "linalg_tools.h"




namespace Cart {
        enum Index {
                    s,                                                                                                                                                     // s
                    x, y, z,                                                                                                                                               // p
                    xx, xy, xz, yy, yz, zz,                                                                                                                                // d
                    xxx, xxy, xxz, xyy, xyz, xzz, yyy, yyz, yzz, zzz,                                                                                                      // f
                    xxxx, xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyy, yyyz, yyzz, yzzz, zzzz,                                                              // g
                    xxxxx, xxxxy, xxxxz, xxxyy, xxxyz, xxxzz, xxyyy, xxyyz, xxyzz, xxzzz, xyyyy, xyyyz, xyyzz, xyzzz, xzzzz, yyyyy, yyyyz, yyyzz, yyzzz, yzzzz, zzzzz,     // h

                    xxxxxx, xxxxxy, xxxxxz, xxxxyy, xxxxyz, xxxxzz, xxxyyy, xxxyyz, xxxyzz, xxxzzz, xxyyyy, xxyyyz, xxyyzz, xxyzzz,                                        // i
                    xxzzzz, xyyyyy, xyyyyz, xyyyzz, xyyzzz, xyzzzz, xzzzzz, yyyyyy, yyyyyz, yyyyzz, yyyzzz, yyzzzz, yzzzzz, zzzzzz,

                    xxxxxxx, xxxxxxy, xxxxxxz, xxxxxyy, xxxxxyz, xxxxxzz, xxxxyyy, xxxxyyz, xxxxyzz, xxxxzzz, xxxyyyy, xxxyyyz,                                            // j
                    xxxyyzz, xxxyzzz, xxxzzzz, xxyyyyy, xxyyyyz, xxyyyzz, xxyyzzz, xxyzzzz, xxzzzzz, xyyyyyy, xyyyyyz, xyyyyzz,
                    xyyyzzz, xyyzzzz, xyzzzzz, xzzzzzz, yyyyyyy, yyyyyyz, yyyyyzz, yyyyzzz, yyyzzzz, yyzzzzz, yzzzzzz, zzzzzzz,

                    xxxxxxxx, xxxxxxxy, xxxxxxxz, xxxxxxyy, xxxxxxyz, xxxxxxzz, xxxxxyyy, xxxxxyyz, xxxxxyzz, xxxxxzzz, xxxxyyyy, xxxxyyyz, xxxxyyzz, xxxxyzzz, xxxxzzzz,  // k
                    xxxyyyyy, xxxyyyyz, xxxyyyzz, xxxyyzzz, xxxyzzzz, xxxzzzzz, xxyyyyyy, xxyyyyyz, xxyyyyzz, xxyyyzzz, xxyyzzzz, xxyzzzzz, xxzzzzzz, xyyyyyyy, xyyyyyyz,
                    xyyyyyzz, xyyyyzzz, xyyyzzzz, xyyzzzzz, xyzzzzzz, xzzzzzzz, yyyyyyyy, yyyyyyyz, yyyyyyzz, yyyyyzzz, yyyyzzzz, yyyzzzzz, yyzzzzzz, yzzzzzzz, zzzzzzzz,
                };
}

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;

    
    
    
    /* "superclass" AOSuperMatrix contains all common functionality for
     * atomic orbital matrix types
     */
        class AOSuperMatrix{
    public:
        
        static int getBlockSize( int size );
        
        
        static ub::matrix<double> getTrafo( const AOGaussianPrimitive* gaussian);
        
        void PrintIndexToFunction(const AOBasis& aobasis);
        
        
    };
    
    
    // base class for 1D atomic orbital matrix types (overlap, Coulomb, ESP)
    class AOMatrix : public AOSuperMatrix {
    public:
        

	// Access functions
	int Dimension(){ return  _aomatrix.size1();};
	ub::matrix<double> &Matrix(){ return _aomatrix ;};

        const ub::matrix<double> &Matrix() const{ return _aomatrix ;};
        
        void Initialize( int size ) {
            this->_aomatrix = ub::zero_matrix<double>(size,size);
        }
        
        void Fill(const AOBasis& aobasis, vec r = vec(0,0,0) , AOBasis* ecp = NULL );
        
        // matrix print 
        void Print( std::string _ident);
        // integrate F
        static std::vector<double> XIntegrate( int _n, double _T );
        // block fill prototype
        virtual void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix,const  AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp = NULL) {} ;

        // ~AOMatrix(){};
    protected:
        ub::matrix<double> _aomatrix; 
        vec _gridpoint;

    };
    
    
    
    /* base class class for 3D atomic orbital matrix types 
     * (in principle, we could make nD and 1D and 3D are just special types)
     */
    class AOMatrix3D : public AOSuperMatrix {
    public:
        
        
        void Initialize( unsigned size ) {
            _aomatrix.resize(3);
            for (int i = 0; i < 3 ; i++){
              _aomatrix[ i ] = ub::zero_matrix<double>(size,size);
            }
        }

        const std::vector<ub::matrix<double> > &Matrix() const{ return _aomatrix ;};

        // matrix print 
        void Print( std::string _ident);

        
        void Fill(const AOBasis& aobasis );

        // block fill prototype
        virtual void FillBlock(std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp = NULL) {} ;

        
        void Cleanup();
    protected:
        std::vector<ub::matrix<double> > _aomatrix; 
        
      //  ~AOMatrix3D();
        
    };
    
    
    
    /* derived class for atomic orbital gradient matrices, required for
     * momentum transition dipoles
     */
    class AOMomentum : public AOMatrix3D { 
        
        //block fill for gradient/momentum operator, implementation in aomomentum.cc
        void FillBlock( std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        
        
    };
    
    
    
    /* derived class for atomic orbital electrical dipole matrices, required for
     * electical transition dipoles
     */
    class AODipole : public AOMatrix3D { 
        
        //block fill for gradient/momentum operator, implementation in aomomentum.cc
        void FillBlock( std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        
        
    };
    
    
    // derived class for atomic orbital nuclear potential
    class AOESP : public AOMatrix{
    public:
        //block fill for overlap, implementation in aoesp.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        //void Print();
        void Fillnucpotential(const AOBasis& aobasis, std::vector<ctp::QMAtom*>& _atoms,bool _with_ecp=false );
        void Fillextpotential(const AOBasis& aobasis, const std::vector<ctp::PolarSeg*>& _sites);
        ub::matrix<double> &getNuclearpotential(){ return _nuclearpotential;}
        const ub::matrix<double> &getNuclearpotential()const{ return _nuclearpotential;}
        ub::matrix<double> &getExternalpotential(){ return _externalpotential;}
        const ub::matrix<double> &getExternalpotential()const{ return _externalpotential;}
        // ~AOESP();
    private:    
        ub::matrix<double> _nuclearpotential;
        ub::matrix<double> _externalpotential;
    };
    
    
    
    // derived class for Effective Core Potentials
    class AOECP : public AOMatrix{
    public:
        //block fill for overlap, implementation in aoesp.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);

        
        typedef boost::multi_array<double, 3> type_3D;
        
        
        ub::matrix<double> calcVNLmatrix(int _lmax_ecp,const vec& posC, const AOGaussianPrimitive* _g_row,const AOGaussianPrimitive* _g_col,const  ub::matrix<int>& _power_ecp,const ub::matrix<double>& _gamma_ecp,const ub::matrix<double>& _pref_ecp   );
        
        
        
        void getBLMCOF(int _lmax_ecp, int _lmax_dft, const vec& pos, type_3D& BLC, type_3D& C  );
        ub::vector<double> CalcNorms( double decay,int size);
        ub::vector<double> CalcInt_r_exp( int nmax, double decay );
    };
    

    
    // derived class for kinetic energy
    class AOKinetic : public AOMatrix{
    public:
        //block fill for overlap, implementation in aokinetic.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, const AOShell* _shell_row, const AOShell* _shell_col , AOBasis* ecp);
 
        
    };
    
    
    // derived class for atomic orbital overlap
    class AOOverlap : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        //void Print();
        
	//        ~AOOverlap();
        
    };
    
    class AODipole_Potential : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        
        void Fillextpotential(const AOBasis& aobasis, const std::vector<ctp::PolarSeg*>& _sites);
        ub::matrix<double> &getExternalpotential(){ return _externalpotential;}
        const ub::matrix<double> &getExternalpotential()const{ return _externalpotential;}
        
        //void Print();
        
	//        ~AOOverlap();
    private: 
        void setAPolarSite(ctp::APolarSite* site){
            apolarsite=site;
        };
        ctp::APolarSite* apolarsite;
        ub::matrix<double> _externalpotential;
    };
    
    class AOQuadrupole_Potential : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        //void Print();
        
        void Fillextpotential(const AOBasis& aobasis, const std::vector<ctp::PolarSeg*>& _sites);
        ub::matrix<double> &getExternalpotential(){ return _externalpotential;}
        const ub::matrix<double> &getExternalpotential()const{ return _externalpotential;}
        //void Print();
        
	//        ~AOOverlap();
    private: 
        void setAPolarSite(ctp::APolarSite* site){
            apolarsite=site;
        };
        ctp::APolarSite* apolarsite;
        ub::matrix<double> _externalpotential;
    };
    
    
 
    
    //derived class for atomic orbital Coulomb interaction
    class AOCoulomb : public AOMatrix{
    public:
        void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        int Symmetrize(const ub::matrix<double>& _gwoverlap_cholesky);
       
        
    private:
        typedef boost::multi_array<double, 3> ma_type;
        typedef ma_type::index index;
        
        
        
  
    };
}}

#endif	/* AOMATRIX_H */

