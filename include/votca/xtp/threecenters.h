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

#ifndef __XTP_THREECENTERS__H
#define	__XTP_THREECENTERS__H

#include <boost/multi_array.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/logger.h>
#include <votca/tools/linalg.h>
#include <votca/xtp/votca_xtp_config.h>

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;
using namespace votca::tools;
/**
* \brief Calculates three electron overlap integrals for GW and DFT.
*
* 
* 
*/





namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    // due to different requirements for the data format for DFT and GW we have two different classes TCMatrix and TCMatrix_dft which inherit from TCrawMatrix
    
    class TCrawMatrix{    
        
    protected:
        
    bool FillThreeCenterOLBlock(  ub::matrix<double> & _subvector, AOShell* _shell, AOShell* _shell_row, AOShell* _shell_col);
    bool FillThreeCenterRepBlock(  ub::matrix<double> & _subvector, AOShell* _shell, AOShell* _shell_row, AOShell* _shell_col);
    //bool FillThreeCenterOLBlock(  ub::matrix<float> & _subvector, AOShell* _shell, AOShell* _shell_row, AOShell* _shell_col);
    void getTrafo(ub::matrix<double>& _trafo, int _lmax, const double& _decay,std::vector<double> contractions);
    
    int getBlockSize( int size );
    
    void XIntegrate( vector<double>& _FmT, const double& _T );
    
   
    };
    

    
    
    
    
    class TCMatrix_dft : public TCrawMatrix{
    public:
    
    void Fill( AOBasis& gwbasis, AOBasis& dftbasis);
    
    void Cleanup();
    
    int getSize(){return _matrix.size();}
    
    std::vector< ub::matrix<double> >& getData(){return  _matrix;}
    ub::matrix<double>& getDatamatrix( int i ){return  _matrix[i];}
    
    private:
        std::vector< ub::matrix<double> > _matrix;
    
        void FillBlock(AOShell* _shell, AOBasis& dftbasis) ; 
        
    };
    
    
    
    
    class TCMatrix : public TCrawMatrix {
    public:
    
        /// returns one level as a constant reference
        // const ub::matrix<double>& operator[](const int i) const { return _matrix[i]; }
        const ub::matrix<float>& operator[](const int i) const { return _matrix[i]; }
        //const ub::matrix<double>& M_bn( const int m) const { return _matrix[m]; }
        /// returns one level as a reference
        //ub::matrix<double>& operator[](const int i) { return _matrix[i]; }
        ub::matrix<float>& operator[](const int i) { return _matrix[i]; }
        //ub::matrix<double>& M_bn( const int m) { return _matrix[m]; }
        
        int size() {  return _matrix.size(); }
        
        //ub::vector< ub::matrix<double> >& matrix() { return this->_matrix ; }

        int get_mmin() { return this->mmin ;}
        int get_mmax() { return this->mmax ;}
        int get_nmin() { return this->nmin ;}
        int get_nmax() { return this->nmax ;}
        int get_mtot() { return this->mtotal ;}
        int get_ntot() { return this->ntotal ;}

        void set_mmin( int i ) { this->mmin = i ;}
        void set_mmax( int i ) { this->mmax = i ;}
        void set_nmin( int i ) { this->nmin = i ;}
        void set_nmax( int i ) { this->nmax = i ;}
        void set_mtot( int i ) { this->mtotal = i ;}
        void set_ntot( int i ) { this->ntotal = i ;}
        

        void Initialize ( int _basissize, int mmin, int mmax, int nmin, int nmax){

            // here as storage indices starting from zero
            set_mmin( mmin  );
            set_mmax( mmax  );
            set_nmin( nmin  );
            set_nmax( nmax  );
            set_mtot( mmax - mmin +1 );
            set_ntot( nmax - nmin +1 );

            
            // vector has mtotal elements
            _matrix.resize( this->get_mtot() );
            
            // each element is a gwabasis-by-n matrix, initialize to zero
            for ( int i = 0; i < this->get_mtot() ; i++){
                //_matrix[i] = ub::zero_matrix<double>(_basissize,ntotal);
                _matrix[i] = ub::zero_matrix<float>(_basissize,ntotal);
            }
        
        }
        
       
        
        void Prune ( int _basissize, int min, int max);
        
        void Print( string _ident);       
        
       
        void Fill( AOBasis& gwbasis, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals );
     

        // matrix print 
       

        void Symmetrize( const ub::matrix<double>& coulomb  );

        // ~TCMatrix();
        
        void Cleanup();
        
        // ub::matrix<double> matrixProd( int m, ub::matrix<double>& matrix);
        
    private:
        
        // store vector of matrices
        //std::vector< ub::matrix<double> > _matrix;
        std::vector< ub::matrix<float> > _matrix;
        
        // band summation indices
        int mmin;
        int mmax;
        int nmin;
        int nmax;
        int ntotal;
        int mtotal;
        
        
        
        
        
        
        // void FillBlock(ub::vector_range< ub::vector< ub::matrix<double> > >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        
        void FillBlock(std::vector< ub::matrix<double> >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        //void FillBlock(std::vector< ub::matrix<float> >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        
     
        
        
    };
    

}}

#endif	/* AOMATRIX_H */

