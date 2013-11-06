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

// UBLAS stops checking types and array bounds if this flag is defined
#define NDEBUG
#define BOOST_UBLAS_NDEBUG

#ifndef _CALC_GWBSE_TOOL_H
#define	_CALC_GWBSE_TOOL_H

#include <votca/ctp/segment.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>

#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/linalg.h>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
// #include <gsl/gsl_eigen.h>
// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_cblas.h>

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
/**
* \brief GWBSE implementation
*
* Requires a first-principles package, i.e. GAUSSIAN installation
*
* Callname: gwbse
*/

class GWBSE : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    GWBSE() {};
   ~GWBSE() {};

    string  Identify() { return "gwbse"; }
    void    Initialize( Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);

    void    CleanUp();

    // int getMlower(){ return mmin -1; };
    // int getMupper(){ return mmax -1; };
    
    
    

private:

    
    //bool   _maverick;
    
    // program tasks
    bool                                _do_qp_diag;
    bool                                _do_bse_singlets;
    bool                                _do_bse_triplets;
    
    // storage tasks
    bool                                _store_qp_pert;
    bool                                _store_qp_diag;
    bool                                _store_bse_singlets;
    bool                                _store_bse_triplets;
    
    
    
    string _outParent;
    string _outMonDir;
    
    string _package;
    Property _package_options;   
    
    string _gwpackage;
    Property _gwpackage_options; 
    
    // basis sets
    string                              _gwbasis_name;
    string                              _dftbasis_name;

    string                              _ranges;          // range types
    unsigned int                        _homo;            // HOMO index
    unsigned int                        _rpamin;
    unsigned int                        _rpamax;
    double                              _rpamaxfactor;    // RPA level range
    unsigned int                        _qpmin;
    unsigned int                        _qpmax;
    unsigned int                        _qptotal;
    double                              _qpminfactor;
    double                              _qpmaxfactor;     // QP level range
    double                              _bseminfactor;
    double                              _bsemaxfactor;
    unsigned int                        _bse_vmin;
    unsigned int                        _bse_vmax;
    unsigned int                        _bse_cmin;
    unsigned int                        _bse_cmax;
    unsigned int                        _bse_size;
    unsigned int                        _bse_vtotal;
    unsigned int                        _bse_ctotal;
         
    double                              _shift;  // pre-shift of DFT energies
    
    
    // RPA related variables and functions
    // container for the epsilon matrix
    ub::vector< ub::matrix<double> > _epsilon;
    // container for frequencies in screening (index 0: real part, index 1: imaginary part)
    ub::matrix<double> _screening_freq;
    
    void RPA_calculate_epsilon( TCMatrix& _Mmn_RPA , ub::matrix<double> _screening_freq , double _shift , ub::vector<double>& _dft_energies  );
    void RPA_prepare_threecenters( TCMatrix& _Mmn_RPA, TCMatrix& _Mmn_full, AOBasis& gwbasis, AOMatrix& gwoverlap, AOMatrix& gwoverlap_inverse     );

    
    // PPM related variables and functions
    ub::matrix<double> _ppm_phi;
    ub::vector<double> _ppm_freq;
    ub::vector<double> _ppm_weight;
    
    void PPM_construct_parameters( ub::matrix<double>& _overlap_cholesky_inverse   );
    
    // Sigma related variables and functions
    ub::matrix<double> _sigma_x; // exchange term
    ub::matrix<double> _sigma_c; // correlation term
    
    void sigma_prepare_threecenters( TCMatrix& _Mmn );
    void sigma_x_setup( TCMatrix& _Mmn );
    void sigma_c_setup( TCMatrix& _Mmn , ub::vector<double>& _edft );
    
    // QP variables and functions
    ub::vector<double> _qp_energies;
    ub::matrix<double> _vxc;
    ub::vector<double> _qp_diag_energies;     // those should be directly stored in 
    ub::matrix<double> _qp_diag_coefficients; // orbitals object, once the interface is set
    void FullQPHamiltonian();
    
    // BSE variables and functions
    ub::matrix<double> _eh_x;
    ub::matrix<double> _eh_d;
    ub::vector<double> _bse_singlet_energies;
    ub::matrix<double> _bse_singlet_coefficients;
    ub::vector<double> _bse_triplet_energies;
    ub::matrix<double> _bse_triplet_coefficients;
    void BSE_x_setup( TCMatrix& _Mmn );
    void BSE_d_setup( TCMatrix& _Mmn );
    void BSE_solve_triplets();
    void BSE_solve_singlets();
    
};


}}

#endif	/* _CALC_GWBSE_TOOL_H */
