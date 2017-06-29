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

#ifndef _VOTCA_XTP_GWBSE_H
#define _VOTCA_XTP_GWBSE_H

#include <votca/ctp/segment.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/threecenters.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/linalg.h>
#include <votca/xtp/votca_config.h>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;
       

/**
         * \brief Electronic excitations from GW-BSE
         *
         * Evaluates electronic excitations in molecular systems based on
         * many-body Green's functions theory within the GW approximation and the
         * Bethe-Salpeter equation. Requires molecular orbitals of the object
         * in GAUSSIAN, NWChem, or TURBOMOLE format.
         * 
         *  B. Baumeier, Y. Ma, D. Andrienko, M. Rohlfing
         *  J. Chem. Theory Comput. 8, 997-1002 (2012)
         * 
         *  B. Baumeier, D. Andrienko, M. Rohlfing
         *  J. Chem. Theory Comput. 8, 2790-2795 (2012) 
         * 
         */

        class GWBSE {
        public:

            GWBSE(Orbitals* orbitals) : _orbitals(orbitals),
            _qp_diag_energies(orbitals->QPdiagEnergies()),
            _qp_diag_coefficients(orbitals->QPdiagCoefficients()),
            _eh_x(orbitals->eh_x()),
            _eh_d(orbitals->eh_d()),
            _bse_singlet_energies(orbitals->BSESingletEnergies()),
            _bse_singlet_coefficients(orbitals->BSESingletCoefficients()),
            _bse_singlet_coefficients_AR(orbitals->BSESingletCoefficientsAR()),
            _bse_triplet_energies(orbitals->BSETripletEnergies()),
            _bse_triplet_coefficients(orbitals->BSETripletCoefficients()){};

            ~GWBSE() {};



            void Initialize(Property *options);

            std::string Identify() {
                return "gwbse";
            }

            void CleanUp();

            // int getMlower(){ return mmin -1; };
            // int getMupper(){ return mmax -1; };

            void setLogger(ctp::Logger* pLog) {
                _pLog = pLog;
            }

            bool Evaluate();

            // interfaces for options getting/setting


            std::string get_gwbasis_name() {
                return _gwbasis_name;
            }

            void set_gwbasis_name(std::string inp) {
                _gwbasis_name = inp;
            }

            std::string get_dftbasis_name() {
                return _dftbasis_name;
            }

            void set_dftbasis_name(std::string inp) {
                _dftbasis_name = inp;
            }

            std::vector< ub::matrix<double> > getExcitedStateDmat(std::string singletortriplet, int state);

            void addoutput(Property *_summary);

        private:

            ctp::Logger *_pLog;



            //bool   _maverick;

            // program tasks
            bool _do_qp_diag;
            bool _do_bse_diag;
            bool _do_bse_singlets;
            bool _do_bse_triplets;

            // storage tasks
            bool _store_qp_pert;
            bool _store_qp_diag;
            bool _store_bse_singlets;
            bool _store_bse_triplets;
            bool _store_eh_interaction;

            // automatic scissors iteration
            bool _iterate_qp;
            bool _qp_converged;

            // options for own Vxc calculation
            bool _doVxc;
            std::string _functional;
            std::string _grid;

            int _openmp_threads;

            // fragment definitions
            int _fragA;
            int _fragB;

            // BSE variant
            bool _do_full_BSE;


            std::string _outParent;
            std::string _outMonDir;

            std::string _package;
            Property _package_options;

            std::string _gwpackage;
            Property _gwpackage_options;

            // basis sets
            std::string _gwbasis_name;
            std::string _dftbasis_name;

            std::string _ranges; // range types
            unsigned int _homo; // HOMO index
            unsigned int _rpamin;
            unsigned int _rpamax;
            double _rpamaxfactor; // RPA level range
            unsigned int _qpmin;
            unsigned int _qpmax;
            unsigned int _qptotal;
            double _qpminfactor;
            double _qpmaxfactor; // QP level range
            double _bseminfactor;
            double _bsemaxfactor;
            double _ScaHFX;

            double _qp_limit; //convergence criteria for qp iteration [Hartree]]
            double _shift_limit; //convergence criteria for shift iteration [Hartree]]
            unsigned int _bse_vmin;
            unsigned int _bse_vmax;
            unsigned int _bse_cmin;
            unsigned int _bse_cmax;
            unsigned int _bse_size;
            unsigned int _bse_vtotal;
            unsigned int _bse_ctotal;
            int _bse_nmax;
            int _bse_nprint;

            double _shift; // pre-shift of DFT energies
            AOBasis _dftbasis;
            ub::matrix<double> _dft_orbitals;

            Orbitals* _orbitals;
            // RPA related variables and functions
            // container for the epsilon matrix
            std::vector< ub::matrix<double> > _epsilon;
            // container for frequencies in screening (index 0: real part, index 1: imaginary part)
            ub::matrix<double> _screening_freq;
            void symmetrize_threecenters(TCMatrix& _Mmn, ub::matrix<double>& _coulomb);
            void RPA_calculate_epsilon(const TCMatrix& _Mmn_RPA);

            void RPA_real(ub::matrix<double>& result, const TCMatrix& _Mmn_RPA, 
                   const double screening_freq);

            void RPA_imaginary(ub::matrix<double>& result, const TCMatrix& _Mmn_RPA, 
                    const double screening_freq);

            void RPA_prepare_threecenters(TCMatrix& _Mmn_RPA, const TCMatrix& _Mmn_full, AOBasis& gwbasis,
                    const AOMatrix& gwoverlap, const ub::matrix<double>& gwoverlap_inverse);


            // PPM related variables and functions
            ub::matrix<double> _ppm_phi;
            ub::vector<double> _ppm_freq;
            ub::vector<double> _ppm_weight;

            void PPM_construct_parameters(const ub::matrix<double>& _overlap_cholesky_inverse,const ub::matrix<double>& _overlap_cholesky_inverse_trans);

            // Sigma related variables and functions
            ub::symmetric_matrix<double> _sigma_x; // exchange term
            ub::symmetric_matrix<double> _sigma_c; // correlation term

            void sigma_prepare_threecenters(TCMatrix& _Mmn);
            void sigma_x_setup(const TCMatrix& _Mmn);
            void sigma_c_setup(const TCMatrix& _Mmn);

            // QP variables and functions
            ub::vector<double> _qp_energies;
            ub::matrix<double> _vxc;
            ub::vector<double>& _qp_diag_energies; // stored in orbitals object 
            ub::matrix<double>& _qp_diag_coefficients; // dito
            void FullQPHamiltonian();

            // BSE variables and functions
            ub::matrix<real_gwbse>& _eh_x; //stored in orbitals object
            ub::matrix<real_gwbse>& _eh_d; //stored in orbitals object
            ub::matrix<real_gwbse> _eh_d2; //because it is not stored in orbitals object
            ub::matrix<real_gwbse> _eh_qp; 

           
            ub::vector<real_gwbse>& _bse_singlet_energies; //stored in orbitals object
            ub::matrix<real_gwbse>& _bse_singlet_coefficients; //stored in orbitals object
            ub::matrix<real_gwbse>& _bse_singlet_coefficients_AR; //stored in orbitals object
            ub::vector<real_gwbse>& _bse_triplet_energies; //stored in orbitals object
            ub::matrix<real_gwbse>& _bse_triplet_coefficients; //stored in orbitals object

            std::vector< ub::matrix<double> > _interlevel_dipoles;
            std::vector< ub::matrix<double> > _interlevel_dipoles_electrical;
            void BSE_x_setup(TCMatrix& _Mmn);
            void BSE_d_setup(TCMatrix& _Mmn);
            void BSE_d2_setup(TCMatrix& _Mmn);
            void BSE_qp_setup();
            void BSE_Add_qp2H(ub::matrix<real_gwbse>& qp);
            void BSE_solve_triplets();
            void BSE_solve_singlets();
            void BSE_solve_singlets_BTDA();
            
            void Solve_nonhermitian(ub::matrix<double>& H, ub::matrix<double>& L);
            std::vector<int> _index2v;
            std::vector<int> _index2c;


            // some cleaner analysis
            void BSE_analyze_triplets();
            void BSE_analyze_singlets();
            void BSE_analyze_singlets_BTDA();

            void BSE_analyze_eh_interaction_Triplet(std::vector<real_gwbse>& _c_d, std::vector<real_gwbse>& _c_qp);
            void BSE_analyze_eh_interaction_Singlet(std::vector<real_gwbse>& _c_x,
                    std::vector<real_gwbse>& _c_d, std::vector<real_gwbse>& _c_qp);
            
            void BSE_analyze_eh_interaction_BTDA_singlet(std::vector<real_gwbse>& _c_x,
                    std::vector<real_gwbse>& _c_d, std::vector<real_gwbse>& _c_qp);


            void BSE_FragmentPopulations(const ub::matrix<real_gwbse>& _bse_coefficients,
            std::vector< ub::vector<double> >& popH, std::vector< ub::vector<double> >& popE,
            std::vector< ub::vector<double> >& Crgs);
            
            void BSE_FragmentPopulations_BTDA(const ub::matrix<real_gwbse>& _bse_coefficients,
            const ub::matrix<real_gwbse>& _bse_coefficients_AR,
            std::vector< ub::vector<double> >& popH, std::vector< ub::vector<double> >& popE,
            std::vector< ub::vector<double> >& Crgs);
            
            
            void BSE_FreeTransition_Dipoles();
            
            void BSE_CoupledTransition_Dipoles();
            
            void BSE_CoupledTransition_Dipoles_BTDA();

        };


    }
}

#endif /* _VOTCA_XTP_GWBSE_H */
