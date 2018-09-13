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
/// For an earlier history see ctp repo commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_EWALDND_H
#define VOTCA_XTP_EWALDND_H

#include <votca/xtp/ewaldactor.h>
#include <votca/xtp/xjob.h>
#include <votca/xtp/xinteractor.h>
#include <votca/xtp/xinductor.h>
#include <votca/xtp/qmthread.h>
#include <boost/multi_array.hpp>

namespace CSG = votca::csg;

namespace votca { namespace xtp {

class PolarTop;    
/// NOTE: This is not a conventional 3D Ewald summation, so use carefully
///       (tuned for the purpose of cluster energy calculations)
/// NOTE: PolarTop should be set-up with three containers: FGC, FGN, BGN
///       MGN is set-up in constructor using the real-space c/o (from input)
///       The topology is used to retrieve information on the sim. box (PB)
///       All polar segments should be positioned as nearest images of the
///       foreground charge density (FGC, FGN).
///       All polar segments should be appropriately charged (Q00, Q10, ...).
/// NOTE: The k-shell grouping algorithm can fail for strongly skewed boxes.
///        

class Ewald3DnD
{

public:

    Ewald3DnD(Topology *top, PolarTop *ptop, Property *opt, Logger *log);
    virtual ~Ewald3DnD();
    virtual std::string IdentifyMethod() = 0;
    
    // POLAR SYSTEM SET-UP
    void ExpandForegroundReduceBackground(double polar_R_co);
    void CoarseGrainDensities(bool cg_bg, bool cg_fg, double cg_radius);
    void SetupMidground(double R_co);
    void WriteDensitiesPDB(std::string pdbfile);
    void WriteDensitiesPtop(std::string fg, std::string mg, std::string bg);
    void WriteInductionStateTable();
    // K-VECTOR GENERATION
    virtual void GenerateKVectors(
        std::vector<PolarSeg*> &ps1, std::vector<PolarSeg*> &ps2) { ; }
    
    // THOLEWALD EVALUATION
    void Evaluate();
    void EvaluateFields(bool do_depolarize_fgc);
    void EvaluateInduction();    
    void EvaluateEnergy(std::vector<PolarSeg*> &target);
    void EvaluateRadialCorrection(std::vector<PolarSeg*> &target);
    void EvaluatePoisson();
    // APERIODIC EMBEDDING QMMM
    bool EvaluateInductionQMMM(bool, bool, bool, bool, bool);
    void EvaluateEnergyQMMM();
    
    // OUTPUT & ERROR COMMUNICATION
    bool Converged() { return _converged_R && _converged_K && _polar_converged; }
    Property GenerateOutputString();
    std::string GenerateErrorString();
    void ShowAgenda(Logger *log);
    void ShowFieldsTeaser(std::vector<PolarSeg*> &target, Logger *log);
    void ShowEnergySplitting(Logger *log);
    
    // ENERGY CALCULATOR METHODS
    virtual EWD::triple<> ConvergeRealSpaceSum(std::vector<PolarSeg*> &target);
    virtual EWD::triple<> ConvergeReciprocalSpaceSum(std::vector<PolarSeg*> &target) = 0;
    virtual EWD::triple<> CalculateForegroundCorrection(std::vector<PolarSeg*> &target);
    virtual EWD::triple<> CalculateHigherRankCorrection(std::vector<PolarSeg*> &target);
    virtual EWD::triple<> CalculateShapeCorrection(std::vector<PolarSeg*> &target)
        { return EWD::triple<>(0.0,0.0,0.0); }
    virtual EWD::triple<> CalculateK0Correction(std::vector<PolarSeg*> &target)
        { return EWD::triple<>(0.0,0.0,0.0); }
    
    // FIELD CALCULATOR METHODS
    virtual void Field_ConvergeRealSpaceSum() { ; }
    virtual void Field_ConvergeReciprocalSpaceSum() { ; }
    virtual void Field_CalculateForegroundCorrection() { ; }
    virtual void Field_CalculateShapeCorrection() { ; }
    
    // POTENTIAL CALCULATOR METHODS
    void EvaluatePotential(std::vector<PolarSeg*> &target, bool add_bg, bool add_mm1, bool add_qm0);
    virtual void Potential_ConvergeRealSpaceSum(std::vector<PolarSeg*> &target) { ; }
    virtual void Potential_ConvergeReciprocalSpaceSum(std::vector<PolarSeg*> &target) { ; }
    virtual void Potential_CalculateForegroundCorrection(std::vector<PolarSeg*> &target) { ; }
    virtual void Potential_CalculateShapeCorrection(std::vector<PolarSeg*> &target) { ; }
    
    // METHOD ANALYSIS
    virtual void ScanCutoff() { ; }    

    // FOREGROUND TRACKER
    class ForegroundTable
    {
    public:
        typedef boost::multi_array<bool,4> fgtable_t;
        typedef fgtable_t::index idx_t;

        ForegroundTable(int n_segs_cell, int na_max, int nb_max, int nc_max)
            : _na_max(na_max), _nb_max(nb_max), _nc_max(nc_max),
              _id_na_nb_nc__inFg(              
                fgtable_t(boost::extents [n_segs_cell]
                                         [2*na_max+1]
                                         [2*nb_max+1]
                                         [2*nc_max+1]))
        {
            for (idx_t i = 0; i < n_segs_cell; ++i) {
            for (idx_t a = 0; a < 2*na_max+1; ++a) {
            for (idx_t b = 0; b < 2*nb_max+1; ++b) {
            for (idx_t c = 0; c < 2*nc_max+1; ++c) {
                _id_na_nb_nc__inFg[i][a][b][c] = false;
            }}}}
        }

        void AddToForeground(int segid, int na, int nb, int nc) {
            idx_t i = segid-1;
            idx_t a = na + _na_max;
            idx_t b = nb + _nb_max;
            idx_t c = nc + _nc_max;
            assert (_id_na_nb_nc__inFg[i][a][b][c] == false);
            _id_na_nb_nc__inFg[i][a][b][c] = true;
        }

        bool IsInForeground(int segid, int na, int nb, int nc) {
            bool is_in_fg;
            if (std::abs(na) > _na_max 
             || std::abs(nb) > _nb_max 
             || std::abs(nc) > _nc_max) {
                is_in_fg = false;
            }
            else {
                idx_t i = segid-1;
                idx_t a = na + _na_max;
                idx_t b = nb + _nb_max;
                idx_t c = nc + _nc_max;
                is_in_fg = _id_na_nb_nc__inFg[i][a][b][c];
            }
            return is_in_fg;
        }

    private:
        //int _n_segs_cell;
        int _na_max;
        int _nb_max;
        int _nc_max;
        fgtable_t _id_na_nb_nc__inFg;
    };


protected:

    EwdInteractor _ewdactor;
    XInteractor _actor;
   
    // PERIODIC BOUNDARY
    Topology *_top;
    vec _center;
    
    // POLAR SEGMENTS
    // Part I - Ewald
    PolarTop *_ptop;

    // [-Wreorder] requires _log to be after _ptop
    Logger *_log;
    bool _started_from_archived_indu_state;
    
    std::vector< PolarSeg* > _bg_P;         // Period. density = BGN + FGN
    std::vector< PolarSeg* > _bg_N;         // Neutral background
    std::vector< PolarSeg* > _mg_N;         // Neutral midground
    std::vector< PolarSeg* > _fg_N;         // Neutral foreground
    std::vector< PolarSeg* > _fg_C;         // Charged foreground
    ForegroundTable *_fg_table;
    std::string _jobType;                   // Calculated from FGC charges
    bool _do_compensate_net_dipole;
    // Part II - Thole
    std::vector< PolarSeg* > _polar_qm0;
    std::vector< PolarSeg* > _polar_mm1;
    std::vector< PolarSeg* > _polar_mm2;    // Should not be used
    double _max_int_dist_qm0;

    // COARSE-GRAINING
    bool _coarse_do_cg_background;
    bool _coarse_do_cg_foreground;
    double _coarse_cg_radius;
    bool _coarse_cg_anisotropic;

    // TASKS
    bool _task_calculate_fields;
    bool _task_polarize_fg;
    bool _task_evaluate_energy;
    bool _task_apply_radial;
    bool _task_solve_poisson;
    bool _task_scan_cutoff;

    // CONVERGENCE
    // Part I - Ewald
    double _alpha;                     // _a = 1/(sqrt(2)*sigma)
    double _kfactor;
    double _rfactor;
    double _K_co;                      // k-space c/o
    double _R_co;                      // r-space c/o
    double _crit_dE;                   // Energy convergence criterion [eV]
    bool   _converged_R;               // Did R-space sum converge?
    bool   _converged_K;               // Did K-space sum converge?
    bool   _field_converged_R;
    bool   _field_converged_K;
    bool   _potential_converged_R;
    bool   _potential_converged_K;
    bool   _did_field_pin_R_shell;
    bool   _save_nblist;
    // Part II - Thole
    bool _polar_do_induce;
    double _polar_aDamp;
    double _polar_wSOR_N;
    double _polar_wSOR_C;
    double _polar_epstol;
    double _polar_cutoff;
    double _polar_converged;
    double _polar_radial_corr_epsilon;

    // LATTICE (REAL, RECIPROCAL)
    vec _a; vec _b; vec _c;            // Real-space lattice vectors
    int _na_max, _nb_max, _nc_max;     // Max. cell indices to sum over (R)
    vec _A; vec _B; vec _C;            // Reciprocal-space lattice vectors
    int _NA_max, _NB_max, _NC_max;     // Max. cell indices to sum over (K)
    double _LxLy;                      // |a^b|
    double _LxLyLz;                    // a*|b^c|
    std::string _shape;                     // Summation shape (for 3D corr. term)

    EWD::VectorSort<EWD::MaxNorm,vec> _maxsort;
    EWD::VectorSort<EWD::EucNorm,vec> _eucsort;
    EWD::VectorSort<EWD::KNorm,EWD::KVector> _kvecsort;

    std::vector<EWD::KVector> _kvecs_2_0;   // K-vectors with two components = 0
    std::vector<EWD::KVector> _kvecs_1_0;   // K-vectors with one component  = 0
    std::vector<EWD::KVector> _kvecs_0_0;   // K-vectors with no  components = 0
    double _kxyz_s1s2_norm;
    bool _did_generate_kvectors;

    // THOLEWALD ENERGIES
    // Part I - Aperiodic Embedding : FGC
    EWD::triple<> _ER;                 // R-space sum
    EWD::triple<> _EC;                 // R-space correction
    EWD::triple<> _EK;                 // K-space sum
    EWD::triple<> _E0;                 // K-space K=0 contribution
    EWD::triple<> _ET;                 // ER - EC + EK + E0
    EWD::triple<> _EDQ;                // Higher-Rank FGC->MGN correction
    EWD::triple<> _EJ;                 // Geometry-dependent correction
    double _outer; double _outer_epp; double _outer_eppu;
    double _inner; double _inner_epp; double _inner_eppu; double _inner_ework;
    // Part I' - Aperiodic Embedding : QM/MM Splitting
    EWD::triple<> _ER_MM1; EWD::triple<> _ER_QM0;
    EWD::triple<> _EC_MM1; EWD::triple<> _EC_QM0;
    EWD::triple<> _EK_MM1; EWD::triple<> _EK_QM0;
    EWD::triple<> _EJ_MM1; EWD::triple<> _EJ_QM0;
    EWD::triple<> _E0_MM1; EWD::triple<> _E0_QM0;
    EWD::triple<> _ET_MM1; EWD::triple<> _ET_QM0;
    // Part II - Thole
    double _polar_ETT;
    double _polar_EPP;  double _polar_EPU;  double _polar_EUU;
    double _polar_EF00; double _polar_EF01; double _polar_EF02;
    double _polar_EF11; double _polar_EF12;
    double _polar_EM0;  double _polar_EM1;  double _polar_EM2;
    // Part I + Part II
    double _Estat;      double _Eindu;      double _Eppuu;
    // Radial dielectric correction
    double _polar_ERC;

    // TIMING (WALL CLOCK)
    double _t_total;
    double _t_coarsegrain;
    double _t_fields;
    double _t_induction;
    double _t_energy;
    double _t_radial;
    

};


}}


#endif // VOTCA_XTP_EWALDND_H
