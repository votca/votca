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
#include <votca/xtp/multiarray.h>


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


    /* "superclass" AOSuperMatrix contains all common functionality for
     * atomic orbital matrix types
     */
        class AOSuperMatrix{
    public:
        
        static int getBlockSize( int _lmax );
        static Eigen::MatrixXd getTrafo( const AOGaussianPrimitive& gaussian);
        void PrintIndexToFunction(const AOBasis& aobasis);
        
        
    };
    
    
    // base class for 1D atomic orbital matrix types (overlap, Coulomb, ESP)
    class AOMatrix : public AOSuperMatrix {
    public:
        
	// Access functions
	int Dimension(){ return  _aomatrix.rows();};
	Eigen::MatrixXd &Matrix(){ return _aomatrix ;};

        const Eigen::MatrixXd &Matrix() const{ return _aomatrix ;};
        
        void Fill(const AOBasis& aobasis, vec r = vec(0,0,0) , AOBasis* ecp = NULL );
        
        // matrix print 
        void Print( std::string _ident);
        // integrate F
        static std::vector<double> XIntegrate( int _n, double _T );
        // block fill prototype
   
        virtual void FillBlock(Eigen::Block<Eigen::MatrixXd>&_matrix,const  AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp = NULL) {} ;

    protected:
        Eigen::MatrixXd _aomatrix; 
        vec _gridpoint;

    };
    
    
    
    /* base class class for 3D atomic orbital matrix types 
     * (in principle, we could make nD and 1D and 3D are just special types)
     */
    class AOMatrix3D : public AOSuperMatrix {
    public:

        const std::vector<Eigen::MatrixXd > &Matrix() const{ return _aomatrix ;};

        // matrix print 
        void Print( std::string _ident);

        void Fill(const AOBasis& aobasis );

        // block fill prototype
      
        virtual void FillBlock(std::vector<Eigen::Block<Eigen::MatrixXd> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp = NULL) {} ;
        
        void Cleanup();
    protected:
        std::vector<Eigen::MatrixXd > _aomatrix; 
        
    };
    
    
    
    /* derived class for atomic orbital gradient matrices, required for
     * momentum transition dipoles
     */
    class AOMomentum : public AOMatrix3D { 
        
        //block fill for gradient/momentum operator, implementation in aomomentum.cc
       
        void FillBlock(std::vector< Eigen::Block<Eigen::MatrixXd> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp = NULL);
        
    };
    
    
    
    /* derived class for atomic orbital electrical dipole matrices, required for
     * electical transition dipoles
     */
    class AODipole : public AOMatrix3D { 
        
        void FillBlock(std::vector< Eigen::Block<Eigen::MatrixXd> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp = NULL);
       
    };
    
    
    // derived class for atomic orbital nuclear potential
    class AOESP : public AOMatrix{
    public:
        //block fill for overlap, implementation in aoesp.cc
        
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix ,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        //void Print();
        void Fillnucpotential(const AOBasis& aobasis, std::vector<QMAtom*>& _atoms);
        void Fillextpotential(const AOBasis& aobasis, const std::vector<ctp::PolarSeg*>& _sites);
        Eigen::MatrixXd &getNuclearpotential(){ return _nuclearpotential;}
        const Eigen::MatrixXd &getNuclearpotential()const{ return _nuclearpotential;}
        Eigen::MatrixXd &getExternalpotential(){ return _externalpotential;}
        const Eigen::MatrixXd &getExternalpotential()const{ return _externalpotential;}
    private:    
        Eigen::MatrixXd _nuclearpotential;
        Eigen::MatrixXd _externalpotential;
    };
    
    
    
    // derived class for Effective Core Potentials
    class AOECP : public AOMatrix{
    public:
        //block fill for overlap, implementation in aoesp.cc
        
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        
        Eigen::MatrixXd calcVNLmatrix(int _lmax_ecp,const vec& posC,
                const AOGaussianPrimitive& _g_row,const AOGaussianPrimitive& _g_col,
                const  Eigen::Matrix<int,4,5>& _power_ecp,const Eigen::Matrix<double,4,5>& _gamma_ecp,
                const Eigen::Matrix<double,4,5>& _pref_ecp   );
        
        void getBLMCOF(int _lmax_ecp, int _lmax_dft, const vec& pos, tensor3d& BLC, tensor3d& C  );
        Eigen::VectorXd CalcNorms( double decay,int size);
        Eigen::VectorXd CalcInt_r_exp( int nmax, double decay );
    };
    

    
    // derived class for kinetic energy
    class AOKinetic : public AOMatrix{
    public:
        //block fill for overlap, implementation in aokinetic.cc
       
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix , const AOShell* _shell_row, const AOShell* _shell_col , AOBasis* ecp);
 
        
    };
    
    
    // derived class for atomic orbital overlap
    class AOOverlap : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
      
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        
    };
    
    class AODipole_Potential : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
       
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        
        void Fillextpotential(const AOBasis& aobasis, const std::vector<ctp::PolarSeg*>& _sites);
        Eigen::MatrixXd &getExternalpotential(){ return _externalpotential;}
        const Eigen::MatrixXd &getExternalpotential()const{ return _externalpotential;}
        
        //void Print();
        
	//        ~AOOverlap();
    private: 
        void setAPolarSite(ctp::APolarSite* site){
            apolarsite=site;
        };
        ctp::APolarSite* apolarsite;
        Eigen::MatrixXd _externalpotential;
    };
    
    class AOQuadrupole_Potential : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
   
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        //void Print();
        
        void Fillextpotential(const AOBasis& aobasis, const std::vector<ctp::PolarSeg*>& _sites);
        Eigen::MatrixXd &getExternalpotential(){ return _externalpotential;}
        const Eigen::MatrixXd &getExternalpotential()const{ return _externalpotential;}
        //void Print();
        
	//        ~AOOverlap();
    private: 
        void setAPolarSite(ctp::APolarSite* site){
            apolarsite=site;
        };
        ctp::APolarSite* apolarsite;
        Eigen::MatrixXd _externalpotential;
    };
    
    
 
    
    //derived class for atomic orbital Coulomb interaction
    class AOCoulomb : public AOMatrix{
    public:
    
        void FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp);
        Eigen::MatrixXd Pseudo_InvSqrt_GWBSE(const AOOverlap& _auxoverlap,double etol);
        Eigen::MatrixXd Pseudo_InvSqrt(double etol);
        int Removedfunctions(){return removedfunctions;}
        
    private:
        
        int removedfunctions;
       
        
        
        
  
    };
}}

#endif	/* AOMATRIX_H */

