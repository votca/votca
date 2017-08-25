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

#ifndef __XTP_THREECENTERS__H
#define	__XTP_THREECENTERS__H
#define BOOST_DISABLE_ASSERTS 
#include <boost/multi_array.hpp>
#include <votca/xtp/aomatrix.h>
//matrix prod overload
#include <votca/tools/linalg.h>
//openmp 
#include <votca/xtp/votca_config.h>
#include <votca/xtp/orbitals.h>
#ifdef _OPENMP
#include <omp.h>
#endif



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
    typedef boost::multi_array<double, 3> ma_type;
    typedef boost::multi_array<ub::matrix<double>, 2> ma2_matrix_type; ////////////////////////////
            typedef boost::multi_array_types::extent_range range; //////////////////
            typedef ma_type::index index; /////////////////////
            ma_type::extent_gen extents; /////////////////////
            
            bool FillThreeCenterRepBlock(  ub::matrix<double> & _subvector, const AOShell* _shell, const AOShell* _shell_row,const AOShell* _shell_col);
    
    };
    
    
    class TCMatrix_dft : public TCrawMatrix{
    public:
    
    void Fill( AOBasis& gwbasis, AOBasis& dftbasis);
    
    void Cleanup();
    
    int getSize(){return _matrix.size();}
    
    std::vector< ub::symmetric_matrix<double> >& getData(){return  _matrix;}
    ub::symmetric_matrix<double>& getDatamatrix( int i ){return  _matrix[i];}
    const ub::symmetric_matrix<double>& getDatamatrix( int i )const{return  _matrix[i];}
    private:
        std::vector< ub::symmetric_matrix<double> > _matrix;
    
        void FillBlock(const AOShell* _shell,const AOBasis& dftbasis) ; 
        
    };


    class FCMatrix_dft : public TCrawMatrix{
    public:
    
    void Fill_4c_small_molecule(const AOBasis& dftbasis); ///////////////////
    
    const ub::vector<double>& get_4c_vector() { return _4c_vector;} ////////////////////
    
    private:
     bool FillFourCenterRepBlock(ub::matrix<double>& _subvector, const AOShell* _shell_1, const AOShell* _shell_2, const AOShell* _shell_3,const AOShell* _shell_4); ////////
    
        ub::vector<double> _4c_vector;
    };

    class TCMatrix : public TCrawMatrix {
    public:
    
        /// returns one level as a constant reference
        const ub::matrix<real_gwbse>& operator[](const int i) const { return _matrix[i]; }
     
        /// returns one level as a reference
        ub::matrix<real_gwbse>& operator[](const int i) { return _matrix[i]; }
        
        int size() {  return _matrix.size(); }
        
        //ub::vector< ub::matrix<double> >& matrix() { return this->_matrix ; }

        int get_mmin() { return mmin ;}
        int get_mmax() { return mmax ;}
        int get_nmin() { return nmin ;}
        int get_nmax() { return nmax ;}
        int get_mtot() { return mtotal ;}
        int get_ntot() { return ntotal ;}
        
        int get_beta()const{return basissize;}
        
        
         int get_mmin() const{ return mmin ;}
         int get_mmax() const{ return mmax ;}
         int get_nmin() const{ return nmin ;}
         int get_nmax() const{ return nmax ;}
         int get_mtot() const{ return mtotal ;}
         int get_ntot() const{ return ntotal ;}

        void set_mmin( int i ) { mmin = i ;}
        void set_mmax( int i ) { mmax = i ;}
        void set_nmin( int i ) { nmin = i ;}
        void set_nmax( int i ) { nmax = i ;}
        void set_mtot( int i ) { mtotal = i ;}
        void set_ntot( int i ) { ntotal = i ;}
        

        void Initialize ( int _basissize, int mmin, int mmax, int nmin, int nmax){

            // here as storage indices starting from zero
            set_mmin( mmin  );
            set_mmax( mmax  );
            set_nmin( nmin  );
            set_nmax( nmax  );
            set_mtot( mmax - mmin +1 );
            set_ntot( nmax - nmin +1 );
            basissize=_basissize;

            
            // vector has mtotal elements
            _matrix.resize( this->get_mtot() );
            
            // each element is a gwabasis-by-n matrix, initialize to zero
            for ( int i = 0; i < this->get_mtot() ; i++){
                _matrix[i] = ub::zero_matrix<real_gwbse>(basissize,ntotal);
            }
        
        }

        void Prune ( int _basissize, int min, int max);
        void Print( std::string _ident);       
        void Fill(const AOBasis& gwbasis,const AOBasis& dftbasis, const ub::matrix<double>& _dft_orbitals );
        
        void Symmetrize( const ub::matrix<double>& coulomb  );
   
        void Cleanup();
        
    private:
        
        // store vector of matrices
        std::vector< ub::matrix<real_gwbse> > _matrix;
        
        // band summation indices
        int mmin;
        int mmax;
        int nmin;
        int nmax;
        int ntotal;
        int mtotal;
        int basissize;
        bool FillThreeCenterOLBlock(  ub::matrix<double> & _subvector, const AOShell* _shell, const AOShell* _shell_row, const AOShell* _shell_col); 
        void FillBlock(std::vector< ub::matrix<double> >& _matrix,const  AOShell* _shell, const AOBasis& dftbasis,const ub::matrix<double>& _dft_orbitals ) ;
        
    };
    

}}

#endif	/* AOMATRIX_H */

