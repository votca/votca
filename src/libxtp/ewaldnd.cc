#include <votca/xtp/ewaldnd.h>
#include <votca/xtp/poissongrid.h>
#include <boost/format.hpp>
#include <algorithm>
#include <boost/math/special_functions/round.hpp>
#include <boost/timer/timer.hpp>


namespace votca { namespace xtp {

using boost::format;

// TODO Develop potential-based scheme k-vector grading
// TODO Deprecate _save_nblist (default = false should always be applied)
// TODO Deprecate K0 term (use shape-term for terminology instead)
// TODO Fields should also accept <target> densities
// TODO Fields should also accept <add_mm1> and <add_qm0> option, or not?
//      (Currently practiced in potentials, which is inconsistent)
// TODO Introduce damping "operator" in EwdActor


Ewald3DnD::~Ewald3DnD() {
    vector< PolarSeg* >::iterator sit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
    
    _polar_qm0.clear();
    _polar_mm1.clear();
    _polar_mm2.clear();
    
    delete _fg_table;
    _fg_table = 0;
}
    
    
Ewald3DnD::Ewald3DnD(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
    :_log(log), _top(top), _ptop(ptop), _fg_table(0) {
    
    // EVALUATE OPTIONS
    string pfx = "options.ewald";
    // Multipoles: started from archive?
    if (opt->exists(pfx+".multipoles.polar_bg")) {
        string ptop_file = opt->get(pfx+".multipoles.polar_bg").as<string>();
        if (ptop_file != "") _started_from_archived_indu_state = true;
        else _started_from_archived_indu_state = false;
    }
    else
        _started_from_archived_indu_state = false;
    // Ewald parameters
    string cmethod = opt->get(pfx+".coulombmethod.method").as<string>();
    assert(cmethod == "ewald" && "<::Ewald3DnD> CMETHOD NOT IMPLEMENTED");
    if (opt->exists(pfx+".coulombmethod.cutoff")) {
        _R_co = opt->get(pfx+".coulombmethod.cutoff").as<double>();
    }
    else
        _R_co = 6.;
    if (opt->exists(pfx+".coulombmethod.shape")) {
        _shape = opt->get(pfx+".coulombmethod.shape").as<string>();
    }
    else
        _shape = "xyslab";
    if (opt->exists(pfx+".coulombmethod.save_nblist")) {
        _save_nblist = opt->get(pfx+".coulombmethod.save_nblist").as<bool>();
    }
    else
        _save_nblist = false;
    // Preprocessing
    if (opt->exists(pfx+".coulombmethod.dipole_corr"))
        _do_compensate_net_dipole = 
            opt->get(pfx+".coulombmethod.dipole_corr").as<bool>();
    else
        _do_compensate_net_dipole = false;
    // Convergence
    if (opt->exists(pfx+".convergence.energy")) {
        _crit_dE = opt->get(pfx+".convergence.energy").as<double>();
    }
    else
        _crit_dE = 1e-5;
    if (opt->exists(pfx+".convergence.kfactor"))
        _kfactor = opt->get(pfx+".convergence.kfactor").as<double>();
    else
        _kfactor = 100.;
    if (opt->exists(pfx+".convergence.rfactor"))
        _rfactor = opt->get(pfx+".convergence.rfactor").as<double>();
    else
        _rfactor = 6.;
    // Polar parameters
    string pmethod = opt->get(pfx+".coulombmethod.method").as<string>();
    assert(pmethod == "ewald" && "<::Ewald3DnD> PMETHOD NOT IMPLEMENTED");
    if (opt->exists(pfx+".polarmethod.induce"))
        _polar_do_induce = opt->get(pfx+".polarmethod.induce").as<bool>();
    else
        _polar_do_induce = false;
    if (opt->exists(pfx+".polarmethod.cutoff")) 
        _polar_cutoff = opt->get(pfx+".polarmethod.cutoff").as<double>();
    else
        _polar_cutoff = 0.0;
    if (opt->exists(pfx+".polarmethod.wSOR_N"))
        _polar_wSOR_N = opt->get(pfx+".polarmethod.wSOR_N").as<double>();
    else
        _polar_wSOR_N = 0.35;
    if (opt->exists(pfx+".polarmethod.wSOR_C"))
        _polar_wSOR_C = opt->get(pfx+".polarmethod.wSOR_C").as<double>();
    else
        _polar_wSOR_C = 0.30;
    if (opt->exists(pfx+".polarmethod.aDamp"))
        _polar_aDamp = opt->get(pfx+".polarmethod.aDamp").as<double>();
    else
        _polar_aDamp = 0.390;
    if (opt->exists(pfx+".polarmethod.tolerance"))
        _polar_epstol = opt->get(pfx+".polarmethod.tolerance").as<double>();
    else
        _polar_epstol = 0.001;
    if (opt->exists(pfx+".polarmethod.radial_dielectric")) {
        _polar_radial_corr_epsilon 
            = opt->get(pfx+".polarmethod.radial_dielectric").as<double>();
    }
    else
        _polar_radial_corr_epsilon = 4.;
    // Coarse-graining
    if (opt->exists(pfx+".coarsegrain.cg_background")) {
        _coarse_do_cg_background = 
            opt->get(pfx+".coarsegrain.cg_background").as<bool>();
    }
    else
        _coarse_do_cg_background = false;
    if (opt->exists(pfx+".coarsegrain.cg_foreground")) {
        _coarse_do_cg_foreground =
            opt->get(pfx+".coarsegrain.cg_foreground").as<bool>();
    }
    else
        _coarse_do_cg_foreground = false;
    if (opt->exists(pfx+".coarsegrain.cg_radius")) {
        _coarse_cg_radius =
            opt->get(pfx+".coarsegrain.cg_radius").as<double>();
    }
    else
        _coarse_cg_radius = _polar_cutoff;
    if (opt->exists(pfx+".coarsegrain.cg_anisotropic")) {
        _coarse_cg_anisotropic =
            opt->get(pfx+".coarsegrain.cg_anisotropic").as<bool>();
    }
    else
        _coarse_cg_anisotropic = false;
    // Tasks to perform
    if (opt->exists(pfx+".tasks.calculate_fields")) {
        _task_calculate_fields 
            = opt->get(pfx+".tasks.calculate_fields").as<bool>();
    }
    else
        _task_calculate_fields = false;
    if (opt->exists(pfx+".tasks.polarize_fg")) {
        _task_polarize_fg = opt->get(pfx+".tasks.polarize_fg").as<bool>();
    }
    else
        _task_polarize_fg = false;
    if (opt->exists(pfx+".tasks.evaluate_energy")) {
        _task_evaluate_energy
            = opt->get(pfx+".tasks.evaluate_energy").as<bool>();
    }
    else
        _task_evaluate_energy = false;
    if (opt->exists(pfx+".tasks.apply_radial")) {
        _task_apply_radial
            = opt->get(pfx+".tasks.apply_radial").as<bool>();
    }
    else
        _task_apply_radial = false;
    if (opt->exists(pfx+".tasks.solve_poisson")) {
        _task_solve_poisson
            = opt->get(pfx+".tasks.solve_poisson").as<bool>();
    }
    else
        _task_solve_poisson = false;
    if (opt->exists(pfx+".tasks.scan_cutoff")) {
        _task_scan_cutoff
            = opt->get(pfx+".tasks.scan_cutoff").as<bool>();
    }
    else
        _task_scan_cutoff = false;
    
    // EWALD INTERACTION PARAMETERS (GUESS ONLY)
    _K_co = _kfactor/_R_co;
    _alpha = _rfactor/_R_co;
    _ewdactor = EwdInteractor(_alpha, _polar_aDamp);
    
    _did_field_pin_R_shell = false;
    _did_generate_kvectors = false;
    
    // SET-UP REAL & RECIPROCAL SPACE
    _a = _top->getBox().getCol(0);
    _b = _top->getBox().getCol(1);
    _c = _top->getBox().getCol(2);
    _LxLyLz = _a*(_b^_c);
    _LxLy = abs(_a ^ _b);
    
    _A = 2*M_PI/_LxLyLz * _b^_c;
    _B = 2*M_PI/_LxLyLz * _c^_a;
    _C = 2*M_PI/_LxLyLz * _a^_b;

    _na_max = ceil((_R_co+_polar_cutoff)/maxnorm(_a)-0.5)+1;
    _nb_max = ceil((_R_co+_polar_cutoff)/maxnorm(_b)-0.5)+1;
    _nc_max = ceil((_R_co+_polar_cutoff)/maxnorm(_c)-0.5)+1;

    _NA_max = ceil(_K_co/maxnorm(_A));
    _NB_max = ceil(_K_co/maxnorm(_B));
    _NC_max = ceil(_K_co/maxnorm(_C));
    
    // SET-UP POLAR GROUNDS (FORE-, MID-, BACK-)
    _center = ptop->getCenter();
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
    
    _fg_C = ptop->FGC();
    _fg_N = ptop->FGN();
    _bg_N = ptop->BGN();
    
    // Construct periodic neutral background
    _bg_P.insert(_bg_P.end(), _fg_N.begin(), _fg_N.end());
    _bg_P.insert(_bg_P.end(), _bg_N.begin(), _bg_N.end());
    // Apply system net dipole compensation if desired
    if (_do_compensate_net_dipole) {
        vector<PolarSeg*>::iterator sit; 
        vector<APolarSite*> ::iterator pit;
        for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
            (*sit)->CalcIsCharged();
        }
        vec system_dpl(0,0,0);
        int charged_count = 0;
        for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
            PolarSeg* pseg = *sit;
            if (!pseg->IsCharged()) continue;
            for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
                charged_count += 1;
                system_dpl += (*pit)->getPos() * (*pit)->getQ00();
                if ((*pit)->getRank() > 0) {
                    system_dpl += (*pit)->getQ1();
                }
            }
        }
        LOG(logINFO,*_log) << "  o System Q1 compensation: " << system_dpl 
            << "  (apply to " << charged_count << " polar sites)" << flush;
        vec atomic_compensation_dpl = - system_dpl/charged_count;
        int compensation_count = 0;
        for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
            PolarSeg* pseg = *sit;
            if (!pseg->IsCharged()) continue;
            for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
                compensation_count += 1;
                if ((*pit)->getRank() < 1) {
                    (*pit)->setQ1(atomic_compensation_dpl);
                    (*pit)->setRank(1);
                }
                else {
                    vec new_dpl = (*pit)->getQ1()+atomic_compensation_dpl;
                    (*pit)->setQ1(new_dpl);
                }
            }
        }
        assert(compensation_count == charged_count);
    }
    // Grow foreground according to induction cut-off
    this->ExpandForegroundReduceBackground(_polar_cutoff);
    // Coarse-grain as demanded by input
    boost::timer::cpu_timer cpu_t;
    cpu_t.start();
    boost::timer::cpu_times t0 = cpu_t.elapsed();
    this->CoarseGrainDensities(_coarse_do_cg_background, 
        _coarse_do_cg_foreground, _coarse_cg_radius);
    boost::timer::cpu_times t1 = cpu_t.elapsed();
    _t_coarsegrain = (t1.wall-t0.wall)/1e9/60.;
    
    // SET-UP MIDGROUND (INCLUDING PERIODIC IMAGES IF REQUIRED)
    LOG(logINFO,*_log) << flush;
    LOG(logINFO,*_log) << "Generate periodic images. ";
    this->SetupMidground(_R_co);
    
    // CALCULATE COG POSITIONS, NET CHARGE
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    double Q_fg_C = 0.0;
    double Q_fg_C_2nd = 0.0;
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;  
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
        double Qseg = (*sit)->CalcTotQ();
        Q_fg_C += Qseg;
        Q_fg_C_2nd += Qseg*Qseg / _fg_C.size();
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_N += (*sit)->CalcTotQ();
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_mg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_P += (*sit)->CalcTotQ();
    }
    
    // DETERMINE JOB TYPE
    int iQ1 = boost::math::iround(Q_fg_C);
    int iQ2 = boost::math::iround(Q_fg_C_2nd);
    if (iQ1 == +1)         _jobType = "hole-like";
    else if (iQ1 == -1)    _jobType = "electron-like";
    else if (iQ1 == 0 && iQ2 == 0) _jobType = "neutral";
    else if (iQ1 == 0 && iQ2 > 0) _jobType = "charge-transfer-like";
    else _jobType = "bipolaron-like";

    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size:")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+5d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+5d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| ~ %2$+5d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+5d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P % _bg_P.size()).str()
        << flush << (format("  o Job type '%3$s' (iQ1=%1$d, iQ2=%2$d)") % iQ1 % iQ2 % _jobType).str()
        << flush;
    
    if (std::abs(Q_bg_P) > 1e-2) {
        cout << endl;
        cout << endl << format("***************************** ERROR ******************************");
        cout << endl << format("       Background charge |Q(BGP)| is larger than 0.01e.");
        cout << endl << format("       Be more precise: e.g. rounding error?");
        cout << endl << format("       Or think again: e.g. erroneous parametrization?");
        cout << endl << format("******************************************************************");
        cout << endl;
    }
    
    // CALCULATE NET DIPOLE OF BGP & FGC
    vec netdpl_bgP = vec(0,0,0);
    vec netdpl_fgC = vec(0,0,0);
    double qzz_bgP = 0.0;
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_bgP += (*pit)->getPos() * (*pit)->getQ00();
            if ((*pit)->getRank() > 0)
                netdpl_bgP += (*pit)->getQ1();
            qzz_bgP += (*pit)->getQ00() * ((*pit)->getPos().getZ() * (*pit)->getPos().getZ());
        }
    }
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_fgC += (*pit)->getPos() * (*pit)->getQ00();
            if ((*pit)->getRank() > 0)
                netdpl_fgC += (*pit)->getQ1();
        }
    }
    
    LOG(logINFO,*_log)
        << (format("Net dipole moment of background density")).str()
        << flush << (format("  o D(BGP) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_bgP.getX() % netdpl_bgP.getY() % netdpl_bgP.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o D(FGC) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") 
        % netdpl_fgC.getX() % netdpl_fgC.getY() % netdpl_fgC.getZ()).str();
    LOG(logINFO,*_log)
        << flush << (format("  o Sigma q|z|**2 [e*nm**2] = %1$+1.7f   ")
        % qzz_bgP) << flush;
    
    
    // ZERO ENERGIES
    _ER  = EWD::triple<>(0, 0, 0);
    _EC  = EWD::triple<>(0, 0, 0);
    _EK  = EWD::triple<>(0, 0, 0);
    _E0  = EWD::triple<>(0, 0, 0);
    _ET  = EWD::triple<>(0, 0, 0);
    _EDQ = EWD::triple<>(0, 0, 0);
    _EJ  = EWD::triple<>(0, 0, 0);
    _polar_ETT = 0;
    _polar_EPP = 0;  _polar_EPU = 0;  _polar_EUU = 0;
    _polar_EF00 = 0; _polar_EF01 = 0; _polar_EF02 = 0;
    _polar_EF11 = 0; _polar_EF12 = 0;
    _polar_EM0 = 0;  _polar_EM1 = 0;  _polar_EM2 = 0;
    _Estat = 0;      _Eindu = 0;      _Eppuu = 0;
    _polar_ERC = 0;
    
    return;
}


void Ewald3DnD::ExpandForegroundReduceBackground(double polar_R_co) {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "Set-up polar grounds" << flush;
    
    vector<PolarSeg*>::iterator sit1;
    vector<PolarSeg*>::iterator sit2;
    vector<PolarSeg*>::iterator sit3;
    
    assert(_polar_mm1.size() == 0 
        && _polar_qm0.size() == 0 
        && _polar_mm2.size() == 0);
    
    // Target containers
    vector<PolarSeg*> exp_fg_C;
    vector<PolarSeg*> exp_fg_N;
    vector<PolarSeg*> red_bg_N;
    
    // Image boxes to consider, set-up boolean foreground table
    int polar_na_max = ceil(polar_R_co/maxnorm(_a)-0.5)+1;
    int polar_nb_max = ceil(polar_R_co/maxnorm(_b)-0.5)+1;
    int polar_nc_max = ceil(polar_R_co/maxnorm(_c)-0.5)+1;
    
    LOG(logDEBUG,*_log) << "  o Expanding cell space for neighbour search:"
            " +/-" << polar_na_max << " x +/-" << polar_nb_max << " x +/-" 
            << polar_nc_max << flush;
    
    _fg_table = new ForegroundTable(
        _bg_P.size(), polar_na_max, polar_nb_max, polar_nc_max);
    
    // Max. distance between any two segments in FGC before expansion
    _max_int_dist_qm0 = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = sit1+1; sit2 < _fg_C.end(); ++sit2) {
            double R = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
            _max_int_dist_qm0 = (R > _max_int_dist_qm0) ? R : _max_int_dist_qm0;
        }
    }
    
    LOG(logDEBUG,*_log) << "  o Max. int. distance (QM0): " << 
        _max_int_dist_qm0 << "nm" << flush;
    
    // Foreground remains in foreground + contributes to QM0
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        exp_fg_C.push_back(*sit1);
        _polar_qm0.push_back(*sit1);
        _fg_table->AddToForeground((*sit1)->getId(), 0, 0, 0);
    }
    for (sit1 = _fg_N.begin(); sit1 < _fg_N.end(); ++sit1) {
        exp_fg_N.push_back(*sit1);
    }    
    
    // Background expands and migrates to foreground OR remains background
    int allocated_count_n = 0;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        PolarSeg *seg_bg = *sit1;
        for (int na = -polar_na_max; na < polar_na_max+1; ++na) {
        for (int nb = -polar_nb_max; nb < polar_nb_max+1; ++nb) {
        for (int nc = -polar_nc_max; nc < polar_nc_max+1; ++nc) {            
            vec L = na*_a + nb*_b + nc*_c;
            
            bool in_central_cell = (na == 0 && nb == 0 && nc == 0);
            bool identical = false;
            bool within_range = false;
            
            for (sit2 = _fg_C.begin(), sit3 = _fg_N.begin(); 
                 sit2 < _fg_C.end(); 
                 ++sit2, ++sit3) {
                PolarSeg *seg_fg = *sit2;
                // Identical ?
                if (in_central_cell && seg_bg->getId() == seg_fg->getId()) {
                    assert(identical==false);
                    identical = true;
                }
                // Within range ?
                // NOTE We have to truncate the decimal places for the radius
                // and the cut-off before we draw the comparison - otherwise
                // numerics may play a prank on us and produce different
                // foregrounds for different charge states, making energy
                // differences for those states pretty much meaningless.
                // Note that abs(dR_L-polar_R_co) < 1e-xy won't do here,
                // since we then to a large degree still rely on machine 
                // precision.
                // NOTE Calculate distance with respect to NEUTRAL segment
                // to achieve consistency between threads in case neutral
                // and charged geometries differ slightly
                double dR_L = votca::tools::abs(
                    seg_bg->getPos() + L - (*sit3)->getPos());
                // Compare up to 3 decimal places (= 1e-3 nm)
                if (int(dR_L*1e3+0.5) <= int(polar_R_co*1e3+0.5)) {
                    within_range = true;
                }
            }
            
            if (!identical && within_range) {
                _fg_table->AddToForeground(seg_bg->getId(), na, nb, nc);
                // Add new shifted clone to fgC and MM1, depolarize = true
                PolarSeg *seg_bg_clone_fgc = new PolarSeg(seg_bg, true);
                seg_bg_clone_fgc->Translate(L);
                exp_fg_C.push_back(seg_bg_clone_fgc);
                _polar_mm1.push_back(seg_bg_clone_fgc);
                // Add original to fgN OR a shifted clone if image box != 0,
                // depolarize = false
                if (in_central_cell)
                    exp_fg_N.push_back(seg_bg);
                else {
                    allocated_count_n += 1;
                    PolarSeg *seg_bg_clone_fgn = new PolarSeg(seg_bg, false);
                    seg_bg_clone_fgn->Translate(L);
                    exp_fg_N.push_back(seg_bg_clone_fgn);
                }
            }
            else if (!identical && !within_range) {
                // Add original to bgN OR a shifted clone if image box != 0
                if (in_central_cell)
                    red_bg_N.push_back(seg_bg);
                else {
                    ;
                    //PolarSeg *seg_bg_clone_bgn = new PolarSeg(seg_bg);
                    //seg_bg_clone_bgn->Translate(L);
                    //red_bg_N.push_back(seg_bg_clone_bgn);
                }
            }
            else {
                ;
            }
        }}} // Sum over image boxes
    } // Sum over periodic neutral background
    
    // Exchange new for old containers
    _fg_C.clear();
    _fg_N.clear();
    _bg_N.clear();
    _fg_C = exp_fg_C;
    _fg_N = exp_fg_N;
    _bg_N = red_bg_N;
    
    bool clean = true;
    _ptop->setFGC(_fg_C, true);
    _ptop->setFGN(_fg_N, true);
    _ptop->setBGN(_bg_N, true);
    clean = false; // Already deleted via fgC
    _ptop->setQM0(_polar_qm0, clean);
    _ptop->setMM1(_polar_mm1, clean);
    _ptop->setMM2(_polar_mm2, clean);    
    
    // Sanity checks
    assert(_polar_qm0.size()+_polar_mm1.size() == _fg_C.size());
    assert(_fg_C.size() == _fg_N.size());
    
    return;
}


void Ewald3DnD::CoarseGrainDensities(bool cg_bg, bool cg_fg, double cg_radius) {
    
    LOG(logDEBUG,*_log) << "Coarse-graining agenda" << flush;
    LOG(logDEBUG,*_log) << "  o Coarse-grain background:   " << ((_coarse_do_cg_background) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Coarse-grain foreground:   " << ((_coarse_do_cg_foreground) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Anisotropic P-tensor:      " << ((_coarse_cg_anisotropic) ? "yes" : "no") << flush;
    LOG(logDEBUG,*_log) << "  o Coarse-graining radius:    " << _coarse_cg_radius << " nm" << flush;
    
    // COARSE-GRAIN BACKGROUND
    if (cg_bg) {
        LOG(logDEBUG,*_log) << "Coarse-grain background" << flush;
        
        int count_bgp = 0;
        int count_fgn = 0;
        //int count_fgc = 0;
        int count_bgp_id_in_fg = 0;

        for (vector<PolarSeg*>::iterator sit = _bg_P.begin();
            sit != _bg_P.end(); ++sit) {
            if (_fg_table->IsInForeground((*sit)->getId(), 0, 0, 0)) {
                ++count_bgp_id_in_fg;
            }
            // By checking whether there are more polar sites than polar
            // fragments in the segment, one avoids coarse-graining twice
            assert((*sit)->size() >= (*sit)->PolarFrags().size()
               && "<::CoarseGrainDensities> BGN: FEWER POLAR SITES THAN FRAGS");
            if ((*sit)->size() > (*sit)->PolarFrags().size()) {
                //cout << "\rMST DBG ...   o BGP ID = " << (*sit)->getId() 
                //     << "   " << flush;
                (*sit)->Coarsegrain(_coarse_cg_anisotropic);
                ++count_bgp;
            }
        }

        for (vector<PolarSeg*>::iterator sit = _fg_N.begin();
            sit != _fg_N.end(); ++sit) {
            // By checking whether there are more polar sites than polar
            // fragments in the segment, one avoids coarse-graining twice
            assert((*sit)->size() >= (*sit)->PolarFrags().size()
               && "<::CoarseGrainDensities> FGN: FEWER POLAR SITES THAN FRAGS");
            if ((*sit)->size() > (*sit)->PolarFrags().size()) {
                //cout << "\rMST DBG ...   o FGN ID = " << (*sit)->getId() 
                //     << "   " << flush;
                (*sit)->Coarsegrain(_coarse_cg_anisotropic);
                ++count_fgn;
            }
        }
        
        LOG(logDEBUG,*_log) << "  o Coarse-grained "
             << count_bgp << "/" << _bg_P.size() << " (BGP) "
             << count_fgn << "/" << _fg_N.size() << " (FGN) "
             << "with " << count_bgp_id_in_fg << " FG-table counts" << flush;
    }
    // COARSE-GRAIN FOREGROUND
    if (cg_fg) {
        LOG(logDEBUG,*_log) << "Coarse-grain foreground" << flush;
        
        int count_fgc = 0;
        
        for (vector<PolarSeg*>::iterator sit = _fg_C.begin();
            sit != _fg_C.end(); ++sit) {
            assert((*sit)->size() >= (*sit)->PolarFrags().size()
               && "<::CoarseGrainDensities> FGC: FEWER POLAR SITES THAN FRAGS");            
            // Only coarse-grain the site if the minimum distance of approach
            // to any of the sites in QM0 is larger than cg_radius
            // Negative radius: Coarse-grain all sites in FGC
            double min_R = std::abs(2*cg_radius);
            for (vector<PolarSeg*>::iterator qit = _polar_qm0.begin();
                qit != _polar_qm0.end(); ++qit) {
                double R = votca::tools::abs((*sit)->getPos()-(*qit)->getPos());
                if (R < min_R) min_R = R;
            }

			if (!(min_R <= _polar_cutoff+1e-3)) {
				cout << endl << "ASSERTION IMMINENT" << endl;
			}
            assert(min_R <= _polar_cutoff+1e-3
                && "<::CoarseGrainDensities> FGC: INCONSISTENT WITH EXPANSION");
            // Different from above, here there should be no danger of 
            // coarse-graining twice
            if (min_R > cg_radius) {
                //cout << "\rMST DBG ...   o FGC ID = " << (*sit)->getId() 
                //     << "   " << flush;
                (*sit)->Coarsegrain(_coarse_cg_anisotropic);
                ++count_fgc;
            }
        }
        
        LOG(logDEBUG,*_log) << "  o Coarse-grained "
             << count_fgc << "/" << _fg_C.size() << " (FGC) " << flush;
    }    
    return;
}


void Ewald3DnD::SetupMidground(double R_co) {
    // SET-UP MIDGROUND
    // NOTE No periodic-boundary correction here: We require that all
    //      polar-segment CoM coords. be folded with respect to the central
    //      segment
    // NOTE Excludes interaction of fg polar segments with neutral selfs in 
    //      real-space sum
    // NOTE Includes periodic images if within cut-off    

    // CLEAR ANY EXTANT MIDGROUND
    vector<PolarSeg*>::iterator sit;
    vector<PolarSeg*>::iterator sit2;
    vector<APolarSite*> ::iterator pit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    _mg_N.clear();
    
    // SAMPLE MIDGROUND FROM BGP EXCLUDING CENTRAL SEG.
    assert(_fg_N.size() == _fg_C.size());
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        // Periodic images
        for (int na = -_na_max; na < _na_max+1; ++na) {
        for (int nb = -_nb_max; nb < _nb_max+1; ++nb) {
        for (int nc = -_nc_max; nc < _nc_max+1; ++nc) {
            vec L = na*_a + nb*_b + nc*_c;
            // In foreground ?
            bool is_in_fg = _fg_table->IsInForeground(pseg->getId(),na,nb,nc);            
            if (!is_in_fg) {
                bool is_within_range = false;
                // Within range ?
                // NOTE Calculate distance with respect to NEUTRAL foreground
                // to achieve consistency between threads in case neutral
                // and charged geometries differ slightly
                for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
                    vec pos_L = pseg->getPos() + L;
                    double dR_L = abs((*sit2)->getPos()-pos_L);
                    if (dR_L <= R_co) {
                        is_within_range = true;
                        break;
                    }
                }
                // Add if appropriate, depolarize = false
                if (is_within_range) {
                    PolarSeg *newSeg = new PolarSeg(pseg, false);
                    newSeg->Translate(L);
                    _mg_N.push_back(newSeg);
                }
            }
            else ;
        }}} // Loop over na, nb, nc
    } // Loop over BGP
    
    return;
}


void Ewald3DnD::WriteDensitiesPDB(string pdbfile) {
    // COORDINATE OUTPUT FOR VISUAL CHECK
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;    
    FILE *out;
    out = fopen(pdbfile.c_str(),"w");
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGC");
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "FGN");
        }
    }
    for (sit = _polar_qm0.begin(); sit < _polar_qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "QM0");
        }
    }
    for (sit = _polar_mm1.begin(); sit < _polar_mm1.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MM1");
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "MGN");
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGN");
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->WritePdbLine(out, "BGP");
        }
    }
    fclose(out);
    return;
}


void Ewald3DnD::WriteDensitiesPtop(string fg, string mg, string bg) {
//    // FGC, FGN, BGN, QM0, MM1, MM2
//    _ptop->SaveToDrive(fg);
//    // MGN
//    PolarTop mg_ptop;
//    mg_ptop.setBGN(_mg_N, false);
//    mg_ptop.SaveToDrive(mg);
//    // BGP
//    PolarTop bg_ptop;
//    bg_ptop.setBGN(_bg_P, false);
//    bg_ptop.SaveToDrive(bg);
//    
//    string fg_pdb = fg + ".pdb";
//    string mg_pdb = mg + ".pdb";
//    string bg_pdb = bg + ".pdb";
//    vector<PolarSeg*>::iterator sit; 
//    vector<APolarSite*> ::iterator pit;    
//    FILE *out;
//    out = fopen(fg_pdb.c_str(),"w");
//    for (sit = _polar_qm0.begin(); sit < _polar_qm0.end(); ++sit) {        
//        PolarSeg* pseg = *sit;        
//        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
//            (*pit)->WritePdbLine(out, "QM0");
//        }
//    }
//    for (sit = _polar_mm1.begin(); sit < _polar_mm1.end(); ++sit) {        
//        PolarSeg* pseg = *sit;        
//        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
//            (*pit)->WritePdbLine(out, "MM1");
//        }
//    }
//    fclose(out);
//    out = fopen(mg_pdb.c_str(),"w");
//    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
//        PolarSeg* pseg = *sit;        
//        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
//            (*pit)->WritePdbLine(out, "MGN");
//        }
//    }
//    fclose(out);
//    out = fopen(bg_pdb.c_str(),"w");
//    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
//        PolarSeg* pseg = *sit;        
//        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
//            (*pit)->WritePdbLine(out, "BGP");
//        }
//    }
//    fclose(out);
//    return;
    

    vector<PolarSeg*>::iterator sit;
    vector<PolarSeg*>::iterator sit2;
    vector<APolarSite*> ::iterator pit;
    assert(_fg_N.size() == _fg_C.size());
    
    std::ofstream ofs;
    ofs.open((fg+bg).c_str(), ofstream::out);
    ofs << "<ptop>\n";

    // Multiplied system: box info
    vec na_a = (2*_na_max+1)*_a;
    vec nb_b = (2*_nb_max+1)*_b;
    vec nc_c = (2*_nc_max+1)*_c;
    ofs << "\t<box>\n";
    ofs << (format("\t\t<a>%1$1.4f %2$1.4f %3$1.4f</a>\n") 
        % na_a.getX() % na_a.getY() % na_a.getZ());
    ofs << (format("\t\t<b>%1$1.4f %2$1.4f %3$1.4f</b>\n") 
        % nb_b.getX() % nb_b.getY() % nb_b.getZ());
    ofs << (format("\t\t<c>%1$1.4f %2$1.4f %3$1.4f</c>\n") 
        % nc_c.getX() % nc_c.getY() % nc_c.getZ());
    ofs << "\t</box>\n";    
    
    // Foreground
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        PolarSeg *pseg = *sit;
        // Position
        vec pos = pseg->getPos();
        // Net induced dipole
        vec u1_tot = vec(0,0,0);                
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit != (*sit)->end(); ++pit) {            
            u1_tot += (*pit)->getU1();     
        }
        // Output: polar segment
        ofs << "\t<pseg>\n";
        ofs << (format("\t\t<id>%1$d</id>\n") % (*sit)->getId());
        ofs << (format("\t\t<size>%1$d</size>\n") % (*sit)->size());
        ofs << (format("\t\t<name>%1$s</name>\n") % _top->getSegment((*sit)->getId())->getName());
        ofs << (format("\t\t<region>fgc</region>\n"));
        ofs << (format("\t\t<pos>%1$1.4f %2$1.4f %3$1.4f</pos>\n") 
            % pos.getX() % pos.getY() % pos.getZ());
        ofs << (format("\t\t<dpl>%1$1.7e %2$1.7e %3$1.7e</dpl>\n") 
            % u1_tot.getX() % u1_tot.getY() % u1_tot.getZ());
        // Output: polar sites in segment
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit != (*sit)->end(); ++pit) {
            vec pos = (*pit)->getPos();
            vec u1 = (*pit)->getU1();
            ofs << "\t\t<psit>\n";
            ofs << (format("\t\t\t<pos>%1$1.4f %2$1.4f %3$1.4f</pos>\n") 
                % pos.getX() % pos.getY() % pos.getZ());
            ofs << (format("\t\t\t<dpl>%1$1.7e %2$1.7e %3$1.7e</dpl>\n") 
                % u1.getX() % u1.getY() % u1.getZ());            
            ofs << "\t\t</psit>\n";
        }
        ofs << "\t</pseg>\n";
    }    
    
    // Background
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        // Periodic images
        for (int na = -_na_max; na < _na_max+1; ++na) {
        for (int nb = -_nb_max; nb < _nb_max+1; ++nb) {
        for (int nc = -_nc_max; nc < _nc_max+1; ++nc) {
            vec L = na*_a + nb*_b + nc*_c;
            // In foreground ?
            bool is_in_fg = _fg_table->IsInForeground(pseg->getId(),na,nb,nc);
            if (!is_in_fg) {
                // Position
                vec pos_L = pseg->getPos() + L;
                // Net induced dipole
                vec u1_tot = vec(0,0,0);                
                for (PolarSeg::iterator pit = (*sit)->begin();
                    pit != (*sit)->end(); ++pit) {            
                    u1_tot += (*pit)->getU1();     
                }
                // Output: polar segment
                ofs << "\t<pseg>\n";
                ofs << (format("\t\t<id>%1$d</id>\n") % (*sit)->getId());
                ofs << (format("\t\t<size>%1$d</size>\n") % (*sit)->size());
                ofs << (format("\t\t<name>%1$s</name>\n") % _top->getSegment((*sit)->getId())->getName());
                ofs << (format("\t\t<region>bgp*\\fgc</region>\n"));
                ofs << (format("\t\t<pos>%1$1.4f %2$1.4f %3$1.4f</pos>\n") 
                    % pos_L.getX() % pos_L.getY() % pos_L.getZ());
                ofs << (format("\t\t<dpl>%1$1.7e %2$1.7e %3$1.7e</dpl>\n") 
                    % u1_tot.getX() % u1_tot.getY() % u1_tot.getZ());
                // Output: polar sites in segment
                for (PolarSeg::iterator pit = (*sit)->begin();
                    pit != (*sit)->end(); ++pit) {
                    vec pos = (*pit)->getPos() + L;
                    vec u1 = (*pit)->getU1();
                    ofs << "\t\t<psit>\n";
                    ofs << (format("\t\t\t<pos>%1$1.4f %2$1.4f %3$1.4f</pos>\n") 
                        % pos.getX() % pos.getY() % pos.getZ());
                    ofs << (format("\t\t\t<dpl>%1$1.7e %2$1.7e %3$1.7e</dpl>\n") 
                        % u1.getX() % u1.getY() % u1.getZ());            
                    ofs << "\t\t</psit>\n";
                }
                ofs << "\t</pseg>\n";
            }
            else ;
        }}} // Loop over na, nb, nc
    } // Loop over BGP
    
    ofs.close();
    
    return;
}


void Ewald3DnD::WriteInductionStateTable() {
    // Field effect
    if (false && tools::globals::verbose) {
        string tabfile = "polarized_"
            + boost::lexical_cast<string>(_polar_qm0[0]->getId())
            + "_" + _jobType + ".tab";
        std::ofstream ofs;
        ofs.open(tabfile.c_str(), ofstream::out);
        vector<PolarSeg*>::iterator sit1, sit2;
        PolarSeg::iterator pit1, pit2;
        for (sit1 = _fg_C.begin(), sit2 = _fg_N.begin(); 
            sit1 < _fg_C.end();
            ++sit1, ++sit2) {
            PolarSeg *pseg_h = *sit1;
            PolarSeg *pseg_n = *sit2;

            vec dr = pseg_h->getPos() - _polar_qm0[0]->getPos();
            double dR = votca::tools::abs(dr);

            assert(pseg_h->getId() == pseg_n->getId());
            assert(pseg_h->size() == pseg_n->size());

            for (pit1 = pseg_h->begin(), pit2 = pseg_n->begin(); 
                pit1 < pseg_h->end(); 
                ++pit1, ++pit2) {

                vec fp_h = (*pit1)->getFieldP();
                vec fu_h = (*pit1)->getFieldU();
                vec ft_h = fp_h+fu_h;

                vec fp_n = (*pit2)->getFieldP();
                vec fu_n = (*pit2)->getFieldU();
                vec ft_n = fp_n+fu_n;

                vec dfp = fp_h - fp_n;
                vec dfu = fu_h - fu_n;
                vec dft = ft_h - ft_n;

                double frac_prj_dfp = dfp*dr/dR / votca::tools::abs(dfp);
                double frac_prj_dfu = dfu*dr/dR / votca::tools::abs(dfu);
                double frac_prj_dft = dft*dr/dR / votca::tools::abs(dft);

                ofs << (format("dR2 %1$+1.7f fput456 %2$+1.7e %3$+1.7e %4$+1.7e "
                        "fprjput8910 %5$+1.7e %6$+1.7e %7$+1.7e\n")
                    % dR % (votca::tools::abs(dfp)*EWD::int2V_m) % (votca::tools::abs(dfu)*EWD::int2V_m) % (votca::tools::abs(dft)*EWD::int2V_m)
                    % frac_prj_dfp % frac_prj_dfu % frac_prj_dft);
            }

        }
        ofs.close();
    }
    
    if (true && tools::globals::verbose) {
        std::ofstream ofs;
        string tabfile = this->IdentifyMethod() + ".indu_state_"
            + boost::lexical_cast<string>(_polar_qm0[0]->getId())
            + "_" + _jobType + ".tab";
        ofs.open(tabfile.c_str(), ofstream::out);
        for (vector<PolarSeg*>::iterator sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
            PolarSeg *pseg = *sit1;
            //Segment *seg = _top->getSegment(pseg->getId());
            for (PolarSeg::iterator pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
                vec fp = (*pit1)->getFieldP();
                vec fu = (*pit1)->getFieldU();
                vec u1 = (*pit1)->getU1();
                vec pos = (*pit1)->getPos();
                ofs << (format("SEGID2 %1$4d   ") % (pseg->getId()));
                ofs << (format("XYZ6 %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (pos.getX())
                        % (pos.getY()) 
                        % (pos.getZ())).str();
                ofs << (format("FP10 %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (fp.getX()*EWD::int2V_m)
                        % (fp.getY()*EWD::int2V_m) 
                        % (fp.getZ()*EWD::int2V_m)).str();
                ofs << (format("FU14 %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (fu.getX()*EWD::int2V_m)
                        % (fu.getY()*EWD::int2V_m) 
                        % (fu.getZ()*EWD::int2V_m)).str();
                ofs << (format("U118 %1$+1.7e %2$+1.7e %3$+1.7e   ") 
                        % (u1.getX())
                        % (u1.getY()) 
                        % (u1.getZ())).str() << endl;
            }
        }
        for (vector<PolarSeg*>::iterator sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            PolarSeg *pseg = *sit1;
            //Segment *seg = _top->getSegment(pseg->getId());
            for (PolarSeg::iterator pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
                vec fp = (*pit1)->getFieldP();
                vec fu = (*pit1)->getFieldU();
                vec u1 = (*pit1)->getU1();
                vec pos = (*pit1)->getPos();
                ofs << (format("SEGID2 %1$4d   ") % (pseg->getId()));
                ofs << (format("XYZ6 %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (pos.getX())
                        % (pos.getY()) 
                        % (pos.getZ())).str();
                ofs << (format("FP10 %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (fp.getX()*EWD::int2V_m)
                        % (fp.getY()*EWD::int2V_m) 
                        % (fp.getZ()*EWD::int2V_m)).str();
                ofs << (format("FU14 %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (fu.getX()*EWD::int2V_m)
                        % (fu.getY()*EWD::int2V_m) 
                        % (fu.getZ()*EWD::int2V_m)).str();
                ofs << (format("U118 %1$+1.7e %2$+1.7e %3$+1.7e   ") 
                        % (u1.getX())
                        % (u1.getY()) 
                        % (u1.getZ())).str() << endl;
            }
        }
        ofs.close();
    }
    
    if (false && tools::globals::verbose) {
        string tabfile = "polarized_"
            + boost::lexical_cast<string>(_polar_qm0[0]->getId())
            + "_" + _jobType + ".tab";
        std::ofstream ofs;
        ofs.open(tabfile.c_str(), ofstream::out);
        vector<PolarSeg*>::iterator sit, sit1, sit2;
        PolarSeg::iterator pit, pit1, pit2;
        for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
            PolarSeg *pseg = *sit;
            // Periodic images
            for (int na = -_na_max; na < _na_max+1; ++na) {
            for (int nb = -_nb_max; nb < _nb_max+1; ++nb) {
            for (int nc = -_nc_max; nc < _nc_max+1; ++nc) {
                vec L = na*_a + nb*_b + nc*_c;
                // In foreground ?
                bool is_in_fg = _fg_table->IsInForeground(pseg->getId(),na,nb,nc);            
                if (!is_in_fg) {
                    for (PolarSeg::iterator pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
                        vec pos = (*pit1)->getPos() + L;
                        vec u1_c = (*pit1)->getU1();
                        vec u1_n = (*pit1)->getU1();
                        ofs << (format("ID %1$4d   ") % (pseg->getId()));
                        ofs << (format("XYZ %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                                % (pos.getX())
                                % (pos.getY()) 
                                % (pos.getZ())).str();
                        ofs << (format("U1n %1$+1.7e %2$+1.7e %3$+1.7e   ") 
                                % (u1_n.getX())
                                % (u1_n.getY()) 
                                % (u1_n.getZ())).str();
                        ofs << (format("U1c %1$+1.7e %2$+1.7e %3$+1.7e   ") 
                                % (u1_c.getX())
                                % (u1_c.getY()) 
                                % (u1_c.getZ())).str();
                        ofs << "BG" << endl;
                    }
                }
                else ;
            }}} // Loop over na, nb, nc
        } // Loop over BGP
        for (sit1 = _fg_C.begin(), sit2 = _fg_N.begin(); 
            sit1 < _fg_C.end();
            ++sit1, ++sit2) {
            PolarSeg *pseg_h = *sit1;
            PolarSeg *pseg_n = *sit2;

            //vec dr = pseg_h->getPos() - _polar_qm0[0]->getPos();
            //double dR = votca::tools::abs(dr);

            assert(pseg_h->getId() == pseg_n->getId());
            assert(pseg_h->size() == pseg_n->size());

            for (pit1 = pseg_h->begin(), pit2 = pseg_n->begin(); 
                pit1 < pseg_h->end(); 
                ++pit1, ++pit2) {
                
                vec pos = (*pit2)->getPos();
                vec u1_n = (*pit2)->getU1();
                vec u1_c = (*pit1)->getU1();
                
                ofs << (format("ID %1$4d   ") % (pseg_n->getId()));
                ofs << (format("XYZ %1$+1.7e %2$+1.7e %3$+1.7e    ") 
                        % (pos.getX())
                        % (pos.getY()) 
                        % (pos.getZ())).str();
                ofs << (format("U1n %1$+1.7e %2$+1.7e %3$+1.7e   ") 
                        % (u1_n.getX())
                        % (u1_n.getY()) 
                        % (u1_n.getZ())).str();
                ofs << (format("U1c %1$+1.7e %2$+1.7e %3$+1.7e   ") 
                        % (u1_c.getX())
                        % (u1_c.getY()) 
                        % (u1_c.getZ())).str();
                ofs << "FG" << endl;
            }

        } // Loop over (FGC, FGN)
        ofs.close();
    }
    return;
}


void Ewald3DnD::ShowAgenda(Logger *log) {
    
    LOG(logDEBUG,*log) << flush;
    
    LOG(logDEBUG,*log) << "System & Ewald parameters (" << IdentifyMethod() << ")" << flush;
    LOG(logDEBUG,*log) << "  o Real-space unit cell:      " << _a << " x " << _b << " x " << _c << flush;
    LOG(logDEBUG,*log) << "  o Real-space c/o (guess):    " << _R_co << " nm" << flush;
    LOG(logDEBUG,*log) << "  o na(max), nb(max), nc(max): " << _na_max << ", " << _nb_max << ", " << _nc_max << flush;
    LOG(logDEBUG,*log) << "  o 1st Brillouin zone:        " << _A << " x " << _B << " x " << _C << flush;
    LOG(logDEBUG,*log) << "  o Reciprocal-space c/o:      " << _K_co << " 1/nm" << flush;
    LOG(logDEBUG,*log) << "  o R-K switching param.       " << _alpha << " 1/nm" << flush;
    LOG(logDEBUG,*log) << "  o Unit-cell volume:          " << _LxLyLz << " nm**3" << flush;
    LOG(logDEBUG,*log) << "  o LxLy (for 3D2D EW):        " << _LxLy << " nm**2" << flush;
    LOG(logDEBUG,*log) << "  o kx(max), ky(max), kz(max): " << _NA_max << ", " << _NB_max << ", " << _NC_max << flush;
    
    LOG(logDEBUG,*log) << "Tasks to perform (" << IdentifyMethod() << ")" << flush;
    LOG(logDEBUG,*log) << "  o Scan interaction range:    " << ((_task_scan_cutoff) ? "yes" : "no") << flush;
    LOG(logDEBUG,*log) << "  o Calculate fg fields:       " << ((_task_calculate_fields) ? "yes" : "no") << flush;
    LOG(logDEBUG,*log) << "  o Polarize foreground:       " << ((_task_polarize_fg) ? "yes" : "no") << flush;
    LOG(logDEBUG,*log) << "  o Evaluate energy:           " << ((_task_evaluate_energy) ? "yes" : "no") << flush;
    LOG(logDEBUG,*log) << "  o Apply radial correction:   " << ((_task_apply_radial) ? "yes" : "no") << flush;
    
    return;
}


void Ewald3DnD::ShowFieldsTeaser(vector<PolarSeg*> &target, Logger *log) {

    int fieldCount = 0;
    for (vector<PolarSeg*>::iterator sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        PolarSeg *pseg = *sit1;
        Segment *seg = _top->getSegment(pseg->getId());
        LOG(logDEBUG,*log) << "ID = " << pseg->getId() << " (" << seg->getName() << ") " << flush;
        for (PolarSeg::iterator pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
            vec fp = (*pit1)->getFieldP();
            vec fu = (*pit1)->getFieldU();
            vec u1 = (*pit1)->getU1();
            LOG(logDEBUG,*log)
               << (format("FPU = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m    ") 
                    % (fp.getX()*EWD::int2V_m+fu.getX()*EWD::int2V_m)
                    % (fp.getY()*EWD::int2V_m+fu.getY()*EWD::int2V_m) 
                    % (fp.getZ()*EWD::int2V_m+fu.getZ()*EWD::int2V_m)).str();
            LOG(logDEBUG,*log)
               << (format("U1* = (%1$+1.7e %2$+1.7e %3$+1.7e) e*nm") 
                    % (u1.getX())
                    % (u1.getY()) 
                    % (u1.getZ())).str() << flush;
            fieldCount += 1;
            if (fieldCount > 10) {
                LOG(logDEBUG,*log)
                    << "FPU = ... ... ..." << flush;
                break;
            }
        }
        if (fieldCount > 10) break;
    }
    return;
}


void Ewald3DnD::ShowEnergySplitting(Logger *log) {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  + EPP(FGC->MGN)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _ER.Sum() % _ER._pp % _ER._pu % _ER._uu).str()
        << flush << (format("  + EKK(FGC->BGP)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EK.Sum() % _EK._pp % _EK._pu % _EK._uu).str()       
        << flush << (format("  - EPP(FGC->FGN)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EC.Sum() % _EC._pp % _EC._pu % _EC._uu).str()
        << flush << (format("  + EK0(FGC->BGP)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _E0.Sum() % _E0._pp % _E0._pu % _E0._uu).str() 
        << flush << (format("  + EDQ(FGC->MGN)  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EDQ.Sum() % _EDQ._pp % _EDQ._pu % _EDQ._uu).str()
        << flush << (format("  + EJ(shape-dep.) = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _EJ.Sum() % _EJ._pp % _EJ._pu % _EJ._uu).str()
        << flush << (format("  = -----------------------------------------------------------------------------------")).str()
        << flush << (format("  + SUM(E) (0,Q,J) = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % _ET.Sum() % _ET._pp % _ET._pu % _ET._uu).str()
        << flush;
    
    LOG(logINFO,*_log)
        << (format("Interaction FGC <> FGC")).str()
        << flush << (format("  + EF [00 01 11]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % (_polar_EF00+_polar_EF01+_polar_EF11) % _polar_EF00 % _polar_EF01 % _polar_EF11).str()
        << flush << (format("  + EF [02 12 22]  = %1$+1.7e = *****ZERO*****  *****ZERO*****  *****ZERO***** eV") 
            % 0.0).str()
        << flush << (format("  + EM [0  1  2 ]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  *****ZERO***** eV") 
            % (_polar_EM0+_polar_EM1) % _polar_EM0 % _polar_EM1).str()
        << flush << (format("  o E  [PP PU UU]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV") 
            % (_polar_EPP+_polar_EPU+_polar_EUU) % _polar_EPP % _polar_EPU % _polar_EUU).str()
        << flush << (format("  = -----------------------------------------------------------------------------------")).str()
        << flush << (format("  + SUM(E)         = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV")
            % _inner % _inner_epp % _inner_eppu % _inner_ework).str()
        << flush;
    
    LOG(logINFO,*_log)            
        << (format("Interaction FGC <> FGC(i) u ***(o)")).str()
        << flush << (format("  + Ei [pp+pu+iw]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  %4$+1.7e eV")
            % _inner % _inner_epp % _inner_eppu % _inner_ework).str()
        << flush << (format("  + Eo [pp+pu+iw]  = %1$+1.7e = %2$+1.7e  %3$+1.7e  ************** eV")
            % _outer % _outer_epp % _outer_eppu).str()
        << flush << (format("  = ===================================================================================")).str()
        << flush << (format("  + E  [stat+ind]  = %1$+1.7e = %2$+1.7e  %3$+1.7e eV")
            % _Eppuu % _Estat % _Eindu).str()
        << flush;
    
    return;
}


void Ewald3DnD::Evaluate() {
    
    this->ShowAgenda(_log);

    // TEASER OUTPUT PERMANENT FIELDS
    LOG(logDEBUG,*_log) << flush << "Background fields (FGN):" << flush;
    this->ShowFieldsTeaser(_fg_N, _log);

    // PERFORM TASKS
    if (_task_scan_cutoff) ScanCutoff();
    boost::timer::cpu_timer cpu_t;
    cpu_t.start();
    boost::timer::cpu_times t0 = cpu_t.elapsed();
    if (_task_calculate_fields) EvaluateFields(true);
    boost::timer::cpu_times t1 = cpu_t.elapsed();
    if (_task_polarize_fg) EvaluateInduction();
    else _polar_converged = true;
    boost::timer::cpu_times t2 = cpu_t.elapsed();
    if (_task_evaluate_energy) EvaluateEnergy(_fg_C);
    boost::timer::cpu_times t3 = cpu_t.elapsed();
    if (_task_apply_radial) EvaluateRadialCorrection(_fg_C);
    boost::timer::cpu_times t4 = cpu_t.elapsed();
    if (_task_solve_poisson) EvaluatePoisson();
    
    _t_fields    = (t1.wall-t0.wall)/1e9/60.;
    _t_induction = (t2.wall-t1.wall)/1e9/60.;
    _t_energy    = (t3.wall-t2.wall)/1e9/60.;
    _t_radial    = (t4.wall-t3.wall)/1e9/60.;
    
//    EvaluateInductionQMMM(true, true, true, true, true);
//    EvaluateInductionQMMM(false, false, false, false, false);
//    EvaluateEnergyQMMM();
    /* 
    EvaluateInductionQMMM(true, true, true, true, true);
    EvaluateEnergyQMMM();
    bool add_bg = true;
    bool add_mm1 = true;
    bool add_qm0 = false;
    EvaluatePotential(_polar_qm0, add_bg, add_mm1, add_qm0);
    EvaluateEnergyQMMM();
    */
    
    // TEASER OUTPUT PERMANENT FIELDS
    LOG(logDEBUG,*_log) << flush << "Foreground fields (FGC):" << flush;
    this->ShowFieldsTeaser(_fg_C, _log);
    
    // COMPUTE ENERGY SPLITTING
    _outer_epp = _ET._pp;
    _outer_eppu = _ET._pu + _ET._uu;
    _outer = _outer_epp + _outer_eppu;
    
    _inner_epp = _polar_EPP;
    _inner_eppu = _polar_EF00+_polar_EF01+_polar_EF02+_polar_EF11+_polar_EF12 - _polar_EPP;
    _inner_ework = _polar_EM0+_polar_EM1;
    _inner = _inner_epp+_inner_eppu+_inner_ework;
    
    _Estat = _outer_epp + _inner_epp;
    _Eindu = _outer_eppu + _inner_eppu + _inner_ework;
    _Eppuu = _Estat + _Eindu;
    
    this->ShowEnergySplitting(_log);
    
    // ADDITIONAL OUTPUT (IF VERBOSE)
    WriteInductionStateTable();
    
    // CLEAN-UP
    for (vector<PolarSeg*>::iterator sit1 = _fg_C.begin(); 
        sit1 != _fg_C.end(); ++sit1) {
        (*sit1)->ClearPolarNbs();
    }
    
    // TIMING
    _t_total = _t_coarsegrain+_t_fields+_t_induction+_t_energy;    
    LOG(logDEBUG,*_log) << flush << (format("Timing (T = %1$1.2f min)") % (_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Coarsegrain> = %1$2.2f%%") % (100*_t_coarsegrain/_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Fields>      = %1$2.2f%%") % (100*_t_fields/_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Induction>   = %1$2.2f%%") % (100*_t_induction/_t_total)) << flush;
    LOG(logDEBUG,*_log) << (format("  o Usage <Energy>      = %1$2.2f%%") % (100*_t_energy/_t_total)) << flush;    
    LOG(logDEBUG,*_log) << flush;
    
    return;
}


void Ewald3DnD::EvaluateFields(bool do_depolarize_fgc) {

	vector<PolarSeg*>::iterator sit;
	vector<APolarSite*> ::iterator pit;

    // RESET FIELDS IF REQUESTED
	if (do_depolarize_fgc) {
		for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
			PolarSeg* pseg = *sit;
			for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
				(*pit)->Depolarize();
			}
		}
	}

    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    Field_ConvergeRealSpaceSum();    

    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    Field_ConvergeReciprocalSpaceSum();

    // SHAPE-CORRECTION (3D3D)/ K0-CORRECTION (3D2D)
    Field_CalculateShapeCorrection();

    // FOREGROUND CORRECTION (3D2D && 3D3D)
    Field_CalculateForegroundCorrection();
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    LOG(logDEBUG,*_log) << flush << "Foreground fields:" << flush;
    int fieldCount = 0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {        
        PolarSeg* pseg = *sit1;        
        for (pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
            vec fp = (*pit1)->getFieldP();
            LOG(logDEBUG,*_log)
               << (format("F = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m") 
                    % (fp.getX()*EWD::int2V_m) 
                    % (fp.getY()*EWD::int2V_m) 
                    % (fp.getZ()*EWD::int2V_m)).str() << flush;
            fieldCount += 1;
            if (fieldCount > 10) {
                LOG(logDEBUG,*_log)
                    << "F = ... ... ..." << flush;
                break;
            }
        }
        if (fieldCount > 10) break;
    }
    
    return;
}


void Ewald3DnD::EvaluateInduction() {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << format("Call inductor on FGC = QM0 u MM1") << flush;
    LOG(logDEBUG,*_log) << (format("  o |QM0|, |MM1|, |MM2|        %1$d %2$d %3$d") 
            % _ptop->QM0().size() % _ptop->MM1().size() % _ptop->MM2().size()).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Polarization cut-off:      ")).str() << _polar_cutoff << " nm " << flush;
    LOG(logDEBUG,*_log) << (format("  o With induction:            %1$s") % ((_polar_do_induce) ? "yes" : "no")) << flush;
    LOG(logDEBUG,*_log) << (format("  o Thole sharpness parameter: ")).str() << _polar_aDamp << flush;
    LOG(logDEBUG,*_log) << (format("  o SOR mixing factor:         ")).str() << _polar_wSOR_N << " (N) " << _polar_wSOR_C << " (C) "  << flush;
    LOG(logDEBUG,*_log) << (format("  o Iterations (max):          512")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Tolerance (dU/U):          %1$1.3e") % _polar_epstol).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Induce within QM0:         yes")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Subthreads:                single")).str() << flush;
    LOG(logDEBUG,*_log) << (format("  o Started from archive:      %1$s") % ((_started_from_archived_indu_state) ? "yes" : "no")) << flush;
    
    
    if (_started_from_archived_indu_state) {
        LOG(logDEBUG,*_log) << "Reusing induction state from archive." << flush;
        vector<PolarSeg*>::iterator sit1, sit2;
        PolarSeg::iterator pit1, pit2;
        for (sit1 = _fg_C.begin(), sit2 = _fg_N.begin(); 
            sit1 < _fg_C.end();
            ++sit1, ++sit2) {
            PolarSeg *pseg_c = *sit1;
            PolarSeg *pseg_n = *sit2;
            assert(pseg_c->getId() == pseg_n->getId());
            assert(pseg_c->size() == pseg_n->size());
            for (pit1 = pseg_c->begin(), pit2 = pseg_n->begin(); 
                pit1 < pseg_c->end(); 
                ++pit1, ++pit2) {
                vec U1_n = (*pit2)->getU1();
                (*pit1)->setU1(U1_n);
            }
        }
    }
    
    
    // Forge XJob object to comply with XInductor interface
    bool polar_has_permanent_fields = true;
    XJob polar_xjob = XJob(_ptop, polar_has_permanent_fields);
    
    // INITIALIZE XINDUCTOR
    bool    polar_induce_intra_pair = true;
    int     polar_subthreads = 1;
    int     polar_maxIter = 512;
    bool    polar_maverick = _log->isMaverick(); // TODO Extract from _log
    
    XInductor polar_xind = XInductor(_polar_do_induce, 
                                     polar_induce_intra_pair, 
                                     polar_subthreads,
                                     _polar_wSOR_N,
                                     _polar_wSOR_C,
                                     _polar_epstol,
                                     polar_maxIter,
                                     _polar_aDamp,
                                     polar_maverick,
                                     _top);
    polar_xind.setLog(_log);
    polar_xind.Evaluate(&polar_xjob);
    
    // SAVE CONVERGENCE
    _polar_converged = polar_xind.hasConverged();
    
    // SAVE RESULTS
    _polar_ETT = polar_xjob.getETOT();
    _polar_EPP = polar_xjob.getEPP();
    _polar_EPU = polar_xjob.getEPU();
    _polar_EUU = polar_xjob.getEUU();
    
    _polar_EF00 = polar_xjob.getEF00();
    _polar_EF01 = polar_xjob.getEF01();
    _polar_EF02 = polar_xjob.getEF02();
    _polar_EF11 = polar_xjob.getEF11();
    _polar_EF12 = polar_xjob.getEF12();
    _polar_EM0 = polar_xjob.getEM0();
    _polar_EM1 = polar_xjob.getEM1();
    _polar_EM2 = polar_xjob.getEM2();
    return;
}


bool Ewald3DnD::EvaluateInductionQMMM(
    bool do_reset, 
    bool do_reuse_bgp_state,
    bool do_calc_perm_fields, 
    bool do_calc_perm_bg_fields, 
    bool do_calc_perm_fg_fields) {
    
    // TO BE TAKEN CARE OF OUTSIDE OF THIS METHOD
    // o Make sure to retain U1 from previous iteration (if existing)
    // o Set polarizabilities of QM0 to zero except for first QMMM iteration
    
    // DEFINITION OF MULTIPOLAR DENSITIES
    //   o  X   = QM0 (excitation)
    //   o  P   = MM1 (polarization cloud)
    //   o  A   = FGC = QM0 u MM1
    //   o  B*~ = BGP \ FGC
    
    // 1st QMMM iteration
    // do_reset = true;
    // do_reuse_bgp_state = true; (ptop-feed required)
    // do_calc_perm_fields = true;
    
    // nth > 1st QMMM iteration
    // do_reset = false;
    // do_reuse_bgp_state = false;
    // do_calc_perm_fields = true;    
    
    vector<PolarSeg*>::iterator sit;
    vector<PolarSeg*>::iterator sit1;
    vector<PolarSeg*>::iterator sit2;
    PolarSeg::iterator pit;
    PolarSeg::iterator pit1;
    PolarSeg::iterator pit2;
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "Started MM Polarization" << flush;
    
    // DEPOLARIZATION STAGE
    if (do_reset) {
        LOG(logDEBUG,*_log) << " o Carry out complete depolarization on FGC" << flush;
        for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->Depolarize();
            }
        }
    }
    else {
        LOG(logDEBUG,*_log) << " o Carry out partial depolarization on FGC" << flush;
        for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->ResetFieldU();
            }
        }
    }
    
    if (do_reuse_bgp_state && _started_from_archived_indu_state) {
        LOG(logDEBUG,*_log) << " o Reuse induction state from archive" << flush;
        vector<PolarSeg*>::iterator sit1, sit2;
        PolarSeg::iterator pit1, pit2;
        for (sit1 = _fg_C.begin(), sit2 = _fg_N.begin(); sit1 < _fg_C.end();
            ++sit1, ++sit2) {
            PolarSeg *pseg_c = *sit1;
            PolarSeg *pseg_n = *sit2;
            assert(pseg_c->getId() == pseg_n->getId());
            assert(pseg_c->size() == pseg_n->size());
            for (pit1 = pseg_c->begin(), pit2 = pseg_n->begin(); pit1 < pseg_c->end(); 
                ++pit1, ++pit2) {
                vec U1_n = (*pit2)->getU1();
                (*pit1)->setU1(U1_n);
            }
        }
    }
    
    // BACKGROUND FIELDS (CONSTANT DURING QMMM ITERATION, STORED AS FP)
    if (do_calc_perm_fields) {
        LOG(logDEBUG,*_log) << flush;
        LOG(logDEBUG,*_log) << "Compute fields for aperiodic embedding" << flush;
        
        if (do_calc_perm_bg_fields) {
			// FIELDS ON 'A' GENERATED BY APERIODIC BACKGROUND : FP[Ql,U1 e B*~] -> X u P
			LOG(logDEBUG,*_log) << " o FP[Ql,U1 e B*~] -> X u P" << flush;
			bool do_depolarize_fgc = false;
			this->EvaluateFields(do_depolarize_fgc);
		}

        if (do_calc_perm_fg_fields) {
            // FIELDS ON 'P' GENERATED BY FOREGROUND (WITHOUT X) : FP[Ql e P] -> P
            LOG(logDEBUG,*_log) << " o FP[Ql e P]      -> P" << flush;
            for (sit1 = _polar_mm1.begin(); sit1 != _polar_mm1.end(); ++sit1) {
            for (sit2 = sit1+1; sit2 != _polar_mm1.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 != (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    _actor.FieldPerm(*(*pit1), *(*pit2));
                }}
            }}

            // FIELDS ON 'X' GENERATED BY FOREGROUND (WITHOUT X) : FP[Ql e P] -> X
            LOG(logDEBUG,*_log) << " o FP[Ql e P]      -> X" << flush;
            for (sit1 = _polar_mm1.begin(); sit1 != _polar_mm1.end(); ++sit1) {
            for (sit2 = _polar_qm0.begin(); sit2 != _polar_qm0.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 != (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit2), *(*pit1));
                    _actor.FieldPermAsPerm_At_By(*(*pit2), *(*pit1));
                }}
            }}
        }

        // FIRST-ORDER INDUCTION
        if (!do_reuse_bgp_state) {
            for (sit = _fg_C.begin(); sit != _fg_C.end(); ++sit) {
                for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                    (*pit)->InduceDirect();
                }
            }
        }
    }
    
    
    // CLASSICAL SCF (ITERATE UNTIL CONVERGED)
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "Start classical SCF with background fields" << flush;
    LOG(logDEBUG,*_log) << flush;
    _log->DisablePreface();
    
    // CONVERGENCE PARAMETERS
    int n_iter   = 1;
    int max_iter = 512;
    bool converged = false;
    double wSOR = _polar_wSOR_N;
    double eTOL = _polar_epstol;
    if (_jobType == "neutral") wSOR = _polar_wSOR_N;
    else wSOR = _polar_wSOR_C;        
    
    for ( ; n_iter < max_iter+1; ++n_iter) {
        
        LOG(logDEBUG,*_log) << (boost::format("o %1$2d") % n_iter) << flush;
        if (n_iter % 10 == 0) {
            _log->EnablePreface();
            LOG(logDEBUG,*_log) << flush;
            _log->DisablePreface();
        }
        
        // Reset induction fields
        for (sit = _fg_C.begin(); sit != _fg_C.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->ResetFieldU();
            }
        }
        
        // QM contribution (changes during QMMM iteration) : FU(Q e X) -> P
        for (sit1 = _polar_mm1.begin(); sit1 != _polar_mm1.end(); ++sit1) {
        for (sit2 = _polar_qm0.begin(); sit2 != _polar_qm0.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
            for (pit2 = (*sit2)->begin(); pit2 != (*sit2)->end(); ++pit2) {
                _actor.BiasIndu(*(*pit1), *(*pit2));
                _actor.FieldPermAsIndu_At_By(*(*pit1), *(*pit2));
            }}
        }}
        
        // Intramolecular contribution : FU(U1 e Ai) -> Ai
        for (sit = _fg_C.begin(); sit != _fg_C.end(); ++sit) {
            for (pit1 = (*sit)->begin(); pit1 != (*sit)->end(); ++pit1) {
            for (pit2 = pit1+1; pit2 != (*sit)->end(); ++pit2) {
                _actor.BiasIndu(*(*pit1), *(*pit2));
                _actor.FieldIndu(*(*pit1), *(*pit2));
            }}
        }
        
        // Intermolecular contribution : FU(U1 e A) -> A
        for (sit1 = _fg_C.begin(); sit1 != _fg_C.end(); ++sit1) {
        for (sit2 = sit1+1; sit2 != _fg_C.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
            for (pit2 = (*sit2)->begin(); pit2 != (*sit2)->end(); ++pit2) {
                _actor.BiasIndu(*(*pit1), *(*pit2));
                _actor.FieldIndu(*(*pit1), *(*pit2));
            }}
        }}
        
        // Induce again
        for (sit = _fg_C.begin(); sit != _fg_C.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->Induce(wSOR);
            }
        }
        
        // Check for convergence
        converged       = true;
        double  maxdU_U = -1;
        double  avgdU_U = 0.0;
        double  rmsdU   = 0.0;
        int     baseN   = 0;
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
             for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                 double dU_U = (*pit1)->HistdU();
                 avgdU_U += dU_U;
                 double dU2 = (*pit1)->HistdU2();
                 rmsdU += dU2;
                 ++baseN;
                 if ( dU_U > maxdU_U ) { maxdU_U = dU_U; }
                 if ( dU_U > eTOL) { converged = false; }
             }
        }
        avgdU_U /= baseN;
        rmsdU /= baseN;
        rmsdU = sqrt(rmsdU);
        if (avgdU_U < eTOL/10.) { converged = true; }        
        
        // Break if converged
        if (converged) { 
            _polar_converged = true;
            break;
        }
        else if (n_iter == max_iter) {
            _polar_converged = false;
            break;
        }       
    }
    
    _log->EnablePreface();
    
    // CLEAN-UP STAGE
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
            (*pit1)->ResetU1Hist();
        }
    }
    
    // ENERGY EVALUATION VIA INDUCTOR    
    // Forge XJob object to comply with XInductor interface
    bool polar_has_permanent_fields = true;
    XJob polar_xjob = XJob(_ptop, polar_has_permanent_fields);
    
    // Initialize XInductor
    bool    polar_induce_intra_pair = true;
    int     polar_subthreads = 1;
    int     polar_maxIter = 512;
    bool    polar_maverick = _log->isMaverick();
    
    XInductor polar_xind = XInductor(_polar_do_induce, 
                                     polar_induce_intra_pair, 
                                     polar_subthreads,
                                     _polar_wSOR_N,
                                     _polar_wSOR_C,
                                     _polar_epstol,
                                     polar_maxIter,
                                     _polar_aDamp,
                                     polar_maverick,
                                     _top);
    polar_xind.setLog(_log);
    polar_xind.Configure(&polar_xjob);
    polar_xind.Energy(&polar_xjob);
    
    // SAVE RESULTS
    _polar_ETT = polar_xjob.getETOT();
    _polar_EPP = polar_xjob.getEPP();
    _polar_EPU = polar_xjob.getEPU();
    _polar_EUU = polar_xjob.getEUU();
    
    _polar_EF00 = polar_xjob.getEF00();
    _polar_EF01 = polar_xjob.getEF01();
    _polar_EF02 = polar_xjob.getEF02();
    _polar_EF11 = polar_xjob.getEF11();
    _polar_EF12 = polar_xjob.getEF12();
    _polar_EM0 = polar_xjob.getEM0();
    _polar_EM1 = polar_xjob.getEM1();
    _polar_EM2 = polar_xjob.getEM2();
    
    return _polar_converged;
}


void Ewald3DnD::EvaluateEnergy(vector<PolarSeg*> &target) {
    
    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    EWD::triple<> EPP_fgC_mgN = ConvergeRealSpaceSum(target);    
    
    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    EWD::triple<> EKK_fgC_bgP = ConvergeReciprocalSpaceSum(target);       
    
    // K=0 TERM (FOR 3D2D)
    EWD::triple<> EK0_fgC_bgP = CalculateK0Correction(target);
    
    // SHAPE-CORRECTION (FOR 3D3D)
    EWD::triple<> EJ_fgC_bgP = CalculateShapeCorrection(target);    
    
    // REAL-SPACE HIGHER-RANK CORRECTION (3D2D && 3D3D)
    EWD::triple<> EDQ_fgC_mgN = CalculateHigherRankCorrection(target); // ! OBSOLETE !
    
    // FOREGROUND CORRECTION (3D2D && 3D3D)
    EWD::triple<> EPP_fgC_fgN = CalculateForegroundCorrection(target);
    
    _ER  = EPP_fgC_mgN * EWD::int2eV;
    _EK  = EKK_fgC_bgP * EWD::int2eV;
    _E0  = EK0_fgC_bgP * EWD::int2eV;
    _EJ  = EJ_fgC_bgP  * EWD::int2eV;
    _EDQ = EDQ_fgC_mgN * EWD::int2eV;
    _EC  = EPP_fgC_fgN * EWD::int2eV;
    
    // NOTE THE (-) IN FRONT OF _EC (WHICH IS A COMPENSATION TERM)
    _ET  = _ER + _EK + _E0 + _EJ + _EDQ - _EC;

    return;
}


void Ewald3DnD::EvaluateEnergyQMMM() {    
    
    _log->setReportLevel(logERROR);
    
    // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
    EWD::triple<> EPP_mm1_mgN = ConvergeRealSpaceSum(_polar_mm1);
    EWD::triple<> EPP_qm0_mgN = ConvergeRealSpaceSum(_polar_qm0);
    
    if (tools::globals::verbose) {
        LOG(logINFO,*_log) << flush;
        LOG(logINFO,*_log) << "R(MM1) " << EPP_mm1_mgN*EWD::int2eV << flush;
        LOG(logINFO,*_log) << "R(MM1) " << EPP_qm0_mgN*EWD::int2eV << flush;
    }    
    
    // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)       
    EWD::triple<> EKK_mm1_bgP = ConvergeReciprocalSpaceSum(_polar_mm1);       
    EWD::triple<> EKK_qm0_bgP = ConvergeReciprocalSpaceSum(_polar_qm0);
    
    if (tools::globals::verbose) {
        LOG(logINFO,*_log) << flush;
        LOG(logINFO,*_log) << "K " << EKK_mm1_bgP*EWD::int2eV << flush;
        LOG(logINFO,*_log) << "K " << EKK_qm0_bgP*EWD::int2eV << flush;
    }
    
    // K=0 TERM (FOR 3D2D)
    EWD::triple<> EK0_mm1_bgP = CalculateK0Correction(_polar_mm1);
    EWD::triple<> EK0_qm0_bgP = CalculateK0Correction(_polar_qm0);
    
    if (tools::globals::verbose) {
        LOG(logINFO,*_log) << flush;
        LOG(logINFO,*_log) << "0 " << EK0_mm1_bgP*EWD::int2eV << flush;
        LOG(logINFO,*_log) << "0 " << EK0_qm0_bgP*EWD::int2eV << flush;
    }
    
    
    // SHAPE-CORRECTION (FOR 3D3D)
    EWD::triple<> EJ_mm1_bgP = CalculateShapeCorrection(_polar_mm1);
    EWD::triple<> EJ_qm0_bgP = CalculateShapeCorrection(_polar_qm0);
    
    if (tools::globals::verbose) {
        LOG(logINFO,*_log) << flush;
        LOG(logINFO,*_log) << "J " << EJ_mm1_bgP*EWD::int2eV << flush;
        LOG(logINFO,*_log) << "J " << EJ_qm0_bgP*EWD::int2eV << flush;
    }
    
    
    // FOREGROUND CORRECTION (3D2D && 3D3D)
    EWD::triple<> EPP_mm1_fgN = CalculateForegroundCorrection(_polar_mm1);
    EWD::triple<> EPP_qm0_fgN = CalculateForegroundCorrection(_polar_qm0);
    
    if (tools::globals::verbose) {
        LOG(logINFO,*_log) << flush;
        LOG(logINFO,*_log) << "C " << EPP_mm1_fgN*EWD::int2eV << flush;
        LOG(logINFO,*_log) << "C " << EPP_qm0_fgN*EWD::int2eV << flush;    
    }
    
    _log->setReportLevel(logDEBUG);
    
    // STORE ENERGIES RESOLVED ACCORDING TO MM1 AND QM0
    _ER_MM1  = EPP_mm1_mgN * EWD::int2eV;
    _ER_QM0  = EPP_qm0_mgN * EWD::int2eV;
    
    _EK_MM1  = EKK_mm1_bgP * EWD::int2eV;
    _EK_QM0  = EKK_qm0_bgP * EWD::int2eV;
    
    _E0_MM1  = EK0_mm1_bgP * EWD::int2eV;
    _E0_QM0  = EK0_qm0_bgP * EWD::int2eV;
    
    _EJ_MM1  = EJ_mm1_bgP  * EWD::int2eV;
    _EJ_QM0  = EJ_qm0_bgP  * EWD::int2eV;
    
    _EC_MM1  = EPP_mm1_fgN * EWD::int2eV;
    _EC_QM0  = EPP_qm0_fgN * EWD::int2eV;
    
    // NOTE THE (-) IN FRONT OF _EC_... (WHICH IS A COMPENSATION TERM)
    _ET_MM1  = _ER_MM1 + _EK_MM1 + _E0_MM1 + _EJ_MM1 - _EC_MM1;
    _ET_QM0  = _ER_QM0 + _EK_QM0 + _E0_QM0 + _EJ_QM0 - _EC_QM0;
    
    // COMPOUND ENERGIES
    _ER = _ER_MM1 + _ER_QM0;
    _EK = _EK_MM1 + _EK_QM0;
    _E0 = _E0_MM1 + _E0_QM0;
    _EJ = _EJ_MM1 + _EJ_QM0;
    _EC = _EC_MM1 + _EC_QM0;
    
    _ET = _ET_MM1 + _ET_QM0;
    
    if (true || tools::globals::verbose) {
        LOG(logINFO,*_log) << flush;
        LOG(logINFO,*_log) << "QM-MM splitting for background interaction (eV)" << flush;
        LOG(logINFO,*_log) << "  o E(MM1)       = " << _ET_MM1 << flush;
        LOG(logINFO,*_log) << "  o E(QM0)       = " << _ET_QM0 << flush;
        LOG(logINFO,*_log) << "  o E(QM0 u MM1) = " << _ET_MM1+_ET_QM0 << flush;
    }
    
    // COMPUTE ENERGY SPLITTING
	_outer_epp = _ET._pp;
	_outer_eppu = _ET._pu + _ET._uu;
	_outer = _outer_epp + _outer_eppu;

	_inner_epp = _polar_EPP;
	_inner_eppu = _polar_EF00+_polar_EF01+_polar_EF02+_polar_EF11+_polar_EF12 - _polar_EPP;
	_inner_ework = _polar_EM0+_polar_EM1;
	_inner = _inner_epp+_inner_eppu+_inner_ework;

	_Estat = _outer_epp + _inner_epp;
	_Eindu = _outer_eppu + _inner_eppu + _inner_ework;
	_Eppuu = _Estat + _Eindu;

    return;
}


void Ewald3DnD::EvaluateRadialCorrection(vector<PolarSeg*> &target) {
    
    // ATTENTION This method depolarizes the midground. Do not call prematurely.
    
    LOG(logINFO,*_log) << flush;
    LOG(logINFO,*_log) << "Apply dielectric radial correction" << flush;
    LOG(logINFO,*_log) << "  o Radial screening constant: " 
        << _polar_radial_corr_epsilon << flush;
    LOG(logINFO,*_log) << "  o Radial extension: " << _polar_cutoff << "nm <> "
        << _polar_cutoff+_R_co << "nm" << flush;
    
    vector<PolarSeg*>::iterator sit1;
    vector<PolarSeg*>::iterator sit2;
    PolarSeg::iterator pit1;
    PolarSeg::iterator pit2;
    
    LOG(logDEBUG,*_log) << "Steps I-IV" << flush;
    // Depolarize midground
    LOG(logDEBUG,*_log) << "  o Depolarize midground" << flush;
    for (sit1 = _mg_N.begin(); sit1 != _mg_N.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
            (*pit1)->Depolarize();
        }
    }    
    // Add screened field contribution from QM0
    LOG(logDEBUG,*_log) << "  o Add screened field contribution (QM0)" << flush;
    for (sit1 = _mg_N.begin(); sit1 != _mg_N.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            for (sit2 = _polar_qm0.begin(); sit2 != _polar_qm0.end(); ++sit2) {
                for (pit2 = (*sit2)->begin(); pit2 != (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    _actor.FieldPerm_At_By(*(*pit1), *(*pit2), _polar_radial_corr_epsilon);
                }
            }
        }
    }    
    // Induce
    LOG(logDEBUG,*_log) << "  o Induce (direct)" << flush;
    for (sit1 = _mg_N.begin(); sit1 != _mg_N.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 != (*sit1)->end(); ++pit1) {
            (*pit1)->InduceDirect();
        }
    }    
    // Evaluate Energy
    LOG(logDEBUG,*_log) << "  o Evaluate stabilization energy" << flush;
    _polar_ERC = 0.0;
    for (sit1 = _mg_N.begin(); sit1 != _mg_N.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            for (sit2 = _polar_qm0.begin(); sit2 != _polar_qm0.end(); ++sit2) {
                for (pit2 = (*sit2)->begin(); pit2 != (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    _polar_ERC += _actor.EnergyInter_Indu_Perm(*(*pit1), *(*pit2));
                }
            }
        }
    }
    _polar_ERC *= EWD::int2eV;
    LOG(logINFO,*_log) << "=> ERC = " << _polar_ERC << "eV" << flush;
    return;
}


void Ewald3DnD::EvaluatePoisson() {
    
    POI::PoissonGrid poisson_grid(_top, _fg_C, _bg_P, _log);
}


void Ewald3DnD::EvaluatePotential(vector<PolarSeg*> &target, bool add_bg, 
    bool add_mm1, bool add_qm0) {
    
    // RESET POTENTIALS
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = target.begin(); sit < target.end(); ++sit) {        
        PolarSeg* pseg = *sit;
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->ResetPhi(true, true);
        }
    }
    
    // APERIODIC-PERIODIC BACKGROUND
    if (add_bg) {
        // REAL-SPACE CONTRIBUTION (3D2D && 3D3D)
        Potential_ConvergeRealSpaceSum(target);    

        // RECIPROCAL-SPACE CONTRIBUTION (3D2D && 3D3D)
        Potential_ConvergeReciprocalSpaceSum(target);

        // SHAPE-CORRECTION (3D3D)/ K0-CORRECTION (3D2D)
        Potential_CalculateShapeCorrection(target);
        
        // FOREGROUND CORRECTION (3D2D && 3D3D)
        Potential_CalculateForegroundCorrection(target);
    }

    // FOREGROUND EXCLUDING QM
    if (add_mm1) {    
        vector<PolarSeg*>::iterator sit1; 
        vector<APolarSite*> ::iterator pit1;
        vector<PolarSeg*>::iterator sit2; 
        vector<APolarSite*> ::iterator pit2;
        for (sit1 = _polar_mm1.begin(); sit1 != _polar_mm1.end(); ++sit1) {
            PolarSeg *pseg1 = *sit1;
            for (sit2 = target.begin(); sit2 != target.end(); ++sit2) {
                PolarSeg *pseg2 = *sit2;
                if (pseg1 == pseg2) assert(false);
                for (pit1 = pseg1->begin(); pit1 != pseg1->end(); ++pit1) {
                    for (pit2 = pseg2->begin(); pit2 != pseg2->end(); ++pit2) {
                        _actor.BiasIndu(*(*pit2), *(*pit1));
                        _actor.Potential_At_By(*(*pit2), *(*pit1));
                    }
                }
            }
        }
    }
    
    // FOREGROUND EXCLUDING MM
    if (add_qm0) {
        vector<PolarSeg*>::iterator sit1; 
        vector<APolarSite*> ::iterator pit1;
        vector<PolarSeg*>::iterator sit2; 
        vector<APolarSite*> ::iterator pit2;
        for (sit1 = _polar_qm0.begin(); sit1 != _polar_qm0.end(); ++sit1) {
            PolarSeg *pseg1 = *sit1;
            for (sit2 = target.begin(); sit2 != target.end(); ++sit2) {
                PolarSeg *pseg2 = *sit2;
                if (pseg1 == pseg2) assert(false);
                for (pit1 = pseg1->begin(); pit1 != pseg1->end(); ++pit1) {
                    for (pit2 = pseg2->begin(); pit2 != pseg2->end(); ++pit2) {
                        _actor.BiasIndu(*(*pit2), *(*pit1));
                        _actor.Potential_At_By(*(*pit2), *(*pit1));
                    }
                }
            }
        }
    }
    
    double q_phi = 0.0;
    for (sit = _polar_qm0.begin(); sit < _polar_qm0.end(); ++sit) {        
        PolarSeg* pseg = *sit;
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            q_phi += (*pit)->getQ00()*(*pit)->getPhi();
        }
    }
    
    LOG(logINFO,*_log) << flush << "Potential q*phi " << q_phi*EWD::int2eV << flush;
    
    return;
}


EWD::triple<> Ewald3DnD::ConvergeRealSpaceSum(vector<PolarSeg*> &target) {
    
    LOG(logDEBUG,*_log) << flush;

    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    // REAL-SPACE CONVERGENCE   
    double dR = 0.1; // radius increment [nm]
    _converged_R = false;
    double prev_ER = 0.0;
    double this_ER = 0.0;
    for (int i = 0; i < 1000; ++i) {
        double Rc = _R_co + (i-1)*dR;
        // Set-up midground
        this->SetupMidground(Rc);
        // Calculate interaction energy
        this_ER = 0.;
        for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
            for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        this_ER += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _alpha);
                    }
                }
            }
        }
        double dER_rms = sqrt((this_ER-prev_ER)*(this_ER-prev_ER))*EWD::int2eV;
        LOG(logDEBUG,*_log)
            << (format("Rc = %1$+1.7f   |MGN| = %3$5d nm   ER = %2$+1.7f eV   dER(rms) = %4$+1.7f") 
            % Rc % (this_ER*EWD::int2eV) % _mg_N.size() % dER_rms).str() << flush;
        if (i > 0 && dER_rms < _crit_dE) {
            _converged_R = true;
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision as of Rc = %1$+1.3f nm") 
                % Rc ) << flush;
            break;
        }
        prev_ER = this_ER;
    }
    return EWD::triple<>(this_ER,0,0);
}


EWD::triple<> Ewald3DnD::CalculateForegroundCorrection(vector<PolarSeg*> &target) {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EPP_fgC_fgN = 0.0;
    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_fgN += _actor.E_QQ_ERF(*(*pit1), *(*pit2), _alpha);
                }
            }
        }
    }
    return EWD::triple<>(EPP_fgC_fgN,0,0);
}


EWD::triple<> Ewald3DnD::CalculateHigherRankCorrection(vector<PolarSeg*> &target) {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EDQ_fgC_mgN = 0.0;
    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    EDQ_fgC_mgN += _actor.E_Q0_DQ(*(*pit1), *(*pit2));
                }
            }
        }
    }
    return EWD::triple<>(EDQ_fgC_mgN,0,0);
}


string Ewald3DnD::GenerateErrorString() {
    string rstr;
    rstr += (format("Converged R-sum = %1$s, converged K-sum = %2$s, ")
        % ((_converged_R) ? "true" : "false")
        % ((_converged_K) ? "true" : "false")).str();
    rstr += (format("converged induction = %1$s")
        % ((_polar_converged) ? "true" : "false")).str();
    return rstr;
}


Property Ewald3DnD::GenerateOutputString() {
    
    Property prop;
    Property &out = prop.add("output","");
    Property *next = NULL;    
    
    next = &out.add("summary", "");
    next->add("type", _jobType);
    next->add("xyz", (format("%1$+1.7f %2$+1.7f %3$+1.7f") 
        % _center.getX() % _center.getY() % _center.getZ()).str())
        .setAttribute("unit","nm");
    next->add("total", (format("%1$+1.7f") 
        % _Eppuu).str())
        .setAttribute("unit","eV");
    next->add("estat", (format("%1$+1.7f") 
        % _Estat).str())
        .setAttribute("unit","eV");
    next->add("eindu", (format("%1$+1.7f") 
        % _Eindu).str())
        .setAttribute("unit","eV");
    
//    next = &out.add("splitting", "");
//    next->add("R-term", (format("%1$+1.7f") % _ER.Sum()).str());
//    next->add("K-term", (format("%1$+1.7f") % _EK.Sum()).str());
//    next->add("O-term", (format("%1$+1.7f") % _E0.Sum()).str());
//    next->add("J-term", (format("%1$+1.7f") % _EJ.Sum()).str());
//    next->add("C-term", (format("%1$+1.7f") % _EC.Sum()).str());
//    next->add("Q-term", (format("%1$+1.7f") % _EDQ.Sum()).str());
    
    next = &out.add("terms_i", "");
    next->add("F-00-01-11", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _polar_EF00 % _polar_EF01 % _polar_EF11).str());
    next->add("M-00-11---", (format("%1$+1.5e %2$+1.5e") % _polar_EM0 % _polar_EM1).str());
    next->add("E-PP-PU-UU", (format("%1$+1.5e %2$+1.5e %3$+1.5e") % _polar_EPP % _polar_EPU % _polar_EUU).str());
    
    next = &out.add("terms_o", "");
    next->add("R-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _ER.Sum() % _ER._pp % _ER._pu % _ER._uu).str());
    next->add("K-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EK.Sum() % _EK._pp % _EK._pu % _EK._uu).str());
    next->add("O-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _E0.Sum() % _E0._pp % _E0._pu % _E0._uu).str());
    next->add("J-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EJ.Sum() % _EJ._pp % _EJ._pu % _EJ._uu).str());
    next->add("C-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EC.Sum() % _EC._pp % _EC._pu % _EC._uu).str());
    next->add("Q-pp-pu-uu", (format("%1$+1.5e = %2$+1.5e %3$+1.5e %4$+1.5e") % _EDQ.Sum() % _EDQ._pp % _EDQ._pu % _EDQ._uu).str());
    
    next = &out.add("terms_c", "");
    next->add("radial_corr", (format("%1$+1.7e") % _polar_ERC).str());
    
    next = &out.add("shells", "");
    next->add("FGC", (format("%1$d") % _fg_C.size()).str());
    next->add("FGN", (format("%1$d") % _fg_N.size()).str());
    next->add("MGN", (format("%1$d") % _mg_N.size()).str());
    next->add("BGN", (format("%1$d") % _bg_N.size()).str());
    next->add("BGP", (format("%1$d") % _bg_P.size()).str());
    next->add("QM0", (format("%1$d") % _polar_qm0.size()).str());
    next->add("MM1", (format("%1$d") % _polar_mm1.size()).str());
    next->add("MM2", (format("%1$d") % _polar_mm2.size()).str());
    
    next = &out.add("timing", "");
    next->add("t_total", (format("%1$1.2f") % _t_total).str())
        .setAttribute("unit","min");
    next->add("t_wload", (format("%1$1.2f %2$1.2f %3$1.2f %4$1.2f")
        % (_t_coarsegrain) % (_t_fields)
        % (_t_induction)   % (_t_energy)).str())
        .setAttribute("unit","min");
        
    return prop;
}
    
    
    
    
    
    
}}
