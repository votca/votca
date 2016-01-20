#include <votca/xtp/pewald3d.h>
#include <boost/format.hpp>
#include <algorithm>
#include <boost/date_time/posix_time/posix_time.hpp>
//#include <boost/timer/timer.hpp>





namespace votca { namespace xtp {

using boost::format;
    

PEwald3D3D::~PEwald3D3D() { ; }
    
    
PEwald3D3D::PEwald3D3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
  : Ewald3DnD(top, ptop, opt, log) {}


void PEwald3D3D::GenerateKVectors(vector<PolarSeg*> &ps1, vector<PolarSeg*> &ps2) {
    
    // Take care of norm for grading function
    // All three components non-zero
    //              S(kx)*S(ky)*S(kz)
    // G = A(k) * ---------------------
    //            (<S(kx)><S(ky)><S(kz)>)**(2/3)
    // Component i zero
    //                   S(kj)*S(kk)
    // G = A(k) * -------------------------
    //             (<S(kj)><S(kk)>)**(1/2)
    // Components i,j zero
    // => All S(k) calculated anyway, no need to grade
    // We can use the same grading function if we set
    //
    // S(ki=0) = <S(ki)>**(2/3) (<S(kj)><S(kk)>)**(1/6)
    
    assert(!_did_generate_kvectors);    
    
    LOG(logINFO,*_log) << flush
        << "Generating K-vectors (Sx*Sy*Sz grading system in place)" << flush;
    
    vector< EWD::KVector > kvecs_2_0; // 2 components zero
    vector< EWD::KVector > kvecs_1_0; // 1 component zero
    vector< EWD::KVector > kvecs_0_0; // 0 components zero
    
    // CONTAINERS FOR GRADING K-VECTORS
    vector< double > kx_s1s2;
    kx_s1s2.push_back(1);
    vector< double > ky_s1s2;
    ky_s1s2.push_back(1);
    vector< double > kz_s1s2;
    kz_s1s2.push_back(1);
    double avg_kx_s1s2 = 0.0;
    double avg_ky_s1s2 = 0.0;
    double avg_kz_s1s2 = 0.0;
    
    // TWO COMPONENTS ZERO, ONE NON-ZERO
    LOG(logINFO,*_log)
        << "  o K-lines through origin: Exploring K resonances" << flush;
    for (int i = 1; i < _NA_max+1; ++i) {
        vec k = +i*_A;
        EWD::triple<EWD::cmplx> ppuu_posk = _ewdactor.S1S2(k, ps1, ps2);        
        kx_s1s2.push_back(0.5*std::abs(ppuu_posk._pp._re));
        avg_kx_s1s2 += 0.5*std::abs(ppuu_posk._pp._re);
        EWD::KVector kvec_pos = EWD::KVector(+1*k,0.);
        EWD::KVector kvec_neg = EWD::KVector(-1*k,0.);
        kvecs_2_0.push_back(kvec_pos);
        kvecs_2_0.push_back(kvec_neg);
    }
    avg_kx_s1s2 /= _NA_max;
    
    for (int i = 1; i < _NB_max+1; ++i) {
        vec k = +i*_B;
        EWD::triple<EWD::cmplx> ppuu_posk = _ewdactor.S1S2(k, ps1, ps2);        
        ky_s1s2.push_back(0.5*std::abs(ppuu_posk._pp._re));
        avg_ky_s1s2 += 0.5*std::abs(ppuu_posk._pp._re);
        EWD::KVector kvec_pos = EWD::KVector(+1*k,0);
        EWD::KVector kvec_neg = EWD::KVector(-1*k,0);
        kvecs_2_0.push_back(kvec_pos);
        kvecs_2_0.push_back(kvec_neg);
    }
    avg_ky_s1s2 /= _NB_max;
    
    for (int i = 1; i < _NC_max+1; ++i) {
        vec k = +i*_C;
        EWD::triple<EWD::cmplx> ppuu_posk = _ewdactor.S1S2(k, ps1, ps2);        
        kz_s1s2.push_back(0.5*std::abs(ppuu_posk._pp._re));
        avg_kz_s1s2 += 0.5*std::abs(ppuu_posk._pp._re);
        EWD::KVector kvec_pos = EWD::KVector(+1*k,0);
        EWD::KVector kvec_neg = EWD::KVector(-1*k,0);
        kvecs_2_0.push_back(kvec_pos);
        kvecs_2_0.push_back(kvec_neg);
    }
    avg_kz_s1s2 /= _NC_max;
    
    double kxyz_s1s2_norm;
    bool neutral_mode = false;
    if (pow(avg_kx_s1s2*avg_ky_s1s2*avg_kz_s1s2,2./3.) > 1e-100)
        kxyz_s1s2_norm = 1./pow(avg_kx_s1s2*avg_ky_s1s2*avg_kz_s1s2,2./3.) * EWD::int2eV / _LxLyLz;
    else {
        LOG(logDEBUG,*_log)
            << "    - Symptoms of a neutral system: Use Ark2Expk2 grading." << flush;
        kxyz_s1s2_norm = EWD::int2eV / _LxLyLz;
        neutral_mode = true;
    }
        
    kx_s1s2[0] = pow(avg_ky_s1s2*avg_kz_s1s2,1./6.)*pow(avg_kx_s1s2,2./3.);
    ky_s1s2[0] = pow(avg_kz_s1s2*avg_kx_s1s2,1./6.)*pow(avg_ky_s1s2,2./3.);
    kz_s1s2[0] = pow(avg_kx_s1s2*avg_ky_s1s2,1./6.)*pow(avg_kz_s1s2,2./3.);
    
    // ONE COMPONENT ZERO, TWO NON-ZERO
    LOG(logINFO,*_log)
        << "  o K-planes through origin: Applying K resonances" << flush;
    
    vector< EWD::KVector >::iterator kvit;
    int kx, ky, kz;
    kx = 0;
    for (ky = -_NB_max; ky < _NB_max+1; ++ky) {
        if (ky == 0) continue;
        for (kz = -_NC_max; kz < _NC_max+1; ++kz) {
            if (kz == 0) continue;
            vec k = kx*_A + ky*_B + kz*_C;
            double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] * ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] * kxyz_s1s2_norm;
            if (neutral_mode) grade = _ewdactor.Ark2Expk2(k) * EWD::int2eV / _LxLyLz;
            EWD::KVector kvec = EWD::KVector(k,grade);
            kvecs_1_0.push_back(kvec);
        }
    }
    ky = 0;
    for (kx = -_NA_max; kx < _NA_max+1; ++kx) {
        if (kx == 0) continue;
        for (kz = -_NC_max; kz < _NC_max+1; ++kz) {
            if (kz == 0) continue;
            vec k = kx*_A + ky*_B + kz*_C;
            double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] * ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] * kxyz_s1s2_norm;
            if (neutral_mode) grade = _ewdactor.Ark2Expk2(k) * EWD::int2eV / _LxLyLz;
            EWD::KVector kvec = EWD::KVector(k,grade);
            kvecs_1_0.push_back(kvec);
        }
    }
    kz = 0;
    for (kx = -_NA_max; kx < _NA_max+1; ++kx) {
        if (kx == 0) continue;
        for (ky = -_NB_max; ky < _NB_max+1; ++ky) {
            if (ky == 0) continue;
            vec k = kx*_A + ky*_B + kz*_C;
            double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] * ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] * kxyz_s1s2_norm;
            if (neutral_mode) grade = _ewdactor.Ark2Expk2(k) * EWD::int2eV / _LxLyLz;
            EWD::KVector kvec = EWD::KVector(k,grade);
            kvecs_1_0.push_back(kvec);
        }
    }
    _kvecsort._p = 1e-300;
    std::sort(kvecs_1_0.begin(), kvecs_1_0.end(), _kvecsort);
    //for (kvit = kvecs_1_0.begin(); kvit < kvecs_1_0.end(); ++kvit) {
    //    EWD::KVector kvec = *kvit;
    //    cout << endl << std::scientific << kvec.getX() << " " << kvec.getY() << " " << kvec.getZ() << " grade " << kvec.getGrade() << flush;
    //}
    
    
    // ZERO COMPONENTS ZERO, THREE NON-ZERO
    LOG(logINFO,*_log)
        << "  o K-space (off-axis): Applying K resonances" << flush;
    
    for (kx = -_NA_max; kx < _NA_max+1; ++kx) {
        if (kx == 0) continue;
        for (ky = -_NB_max; ky < _NB_max+1; ++ky) {
            if (ky == 0) continue;
            for (kz = -_NC_max; kz < _NC_max+1; ++kz) {
                if (kz == 0) continue;
                vec k = kx*_A + ky*_B + kz*_C;
                double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] * ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] * kxyz_s1s2_norm;
                if (neutral_mode) grade = _ewdactor.Ark2Expk2(k) * EWD::int2eV / _LxLyLz;
                EWD::KVector kvec = EWD::KVector(k,grade);
                kvecs_0_0.push_back(kvec);
            }
        }    
    }
    
    _kvecsort._p = 1e-300;
    std::sort(kvecs_0_0.begin(), kvecs_0_0.end(), _kvecsort);
    //for (kvit = kvecs_0_0.begin(); kvit < kvecs_0_0.end(); ++kvit) {
    //    EWD::KVector kvec = *kvit;
    //    cout << endl << std::scientific << kvec.getX() << " " << kvec.getY() << " " << kvec.getZ() << " grade " << kvec.getGrade() << flush;
    //}    
    
    _kvecs_2_0.clear();
    _kvecs_1_0.clear();
    _kvecs_0_0.clear();
    
    _kvecs_2_0 = kvecs_2_0;
    _kvecs_1_0 = kvecs_1_0;
    _kvecs_0_0 = kvecs_0_0;
    _kxyz_s1s2_norm = kxyz_s1s2_norm;
    
    _did_generate_kvectors = true;
    return;
}


void PEwald3D3D::ScanCutoff() {
    /*
    double sum = 0.0;
    double sum_ppuu = 0.0;
    
    LOG(logDEBUG,*_log) << flush 
        << "Scan cutoff (long-range check)" << flush;
    
    std::ofstream ofs;
    string scan_file = "scan_cut_"
        + boost::lexical_cast<string>(_polar_qm0[0]->getId())
        + "_" + _jobType + ".tab";
    ofs.open(scan_file.c_str(), ofstream::out);
    
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    vector<PolarNb*>::iterator nit;
    vector< vector<PolarSeg*> > ::iterator vsit;
    
    double dR_shell = 0.5;
    double R_overhead = 1.1;
    double R_add = 3;
    double R_max = _R_co*R_overhead+R_add;
    double R_max_shell = R_max+2*_polar_cutoff+_max_int_dist_qm0;
    this->SetupMidground(R_max);
    */
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    double sum = 0.0;
    double sum_ppuu = 0.0;
    
    LOG(logDEBUG,*_log) << flush 
        << "Scan cutoff (long-range check)" << flush;
    
    std::ofstream ofs;
    string scan_file = "scan_cut_"
        + boost::lexical_cast<string>(_polar_qm0[0]->getId())
        + "_" + _jobType + ".tab";
    ofs.open(scan_file.c_str(), ofstream::out);
    
    
    vector<PolarSeg*>::iterator sit;
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    vector<PolarNb*>::iterator nit;
    vector< vector<PolarSeg*> > ::iterator vsit;
    
    double R_factor = 1.;
    double dR_shell = 0.5;
    double R_overhead = 1.1;
    double R_add = 3;
    double R_max = R_factor*_R_co*R_overhead+R_add;
    double R_max_shell = R_max+2*_polar_cutoff+_max_int_dist_qm0;

	R_max = 64;
	R_max_shell = 70;
    
    vector<TinyNeighbour*> nbs;
    vector<TinyNeighbour*>::iterator tnit;
    nbs.reserve(1000000);
    
    int scan_na_max = ceil((R_max_shell)/maxnorm(_a)-0.5)+1;
    int scan_nb_max = ceil((R_max_shell)/maxnorm(_b)-0.5)+1;
    int scan_nc_max = ceil((R_max_shell)/maxnorm(_c)-0.5)+1;
    
    if (this->_shape == "xyslab") scan_nc_max = 0;
    
    // loop, create images, bin, interact
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        // Periodic images
        for (int na = -scan_na_max; na < scan_na_max+1; ++na) {
        for (int nb = -scan_nb_max; nb < scan_nb_max+1; ++nb) {
        for (int nc = -scan_nc_max; nc < scan_nc_max+1; ++nc) {
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
                    if (dR_L <= R_max_shell) {
                        is_within_range = true;
                    }
                }
                // Add if appropriate, depolarize = false
                if (is_within_range) {
                    TinyNeighbour *newNb = new TinyNeighbour(pseg, L);
                    nbs.push_back(newNb);
                }
            }
            else ;
        }}} // Loop over na, nb, nc
    } // Loop over BGP
    
    for (sit1 = _fg_C.begin(); sit1 != _fg_C.end(); ++sit1) {

        // Bin midground into shells
        vector< vector<TinyNeighbour*> > shelled_nbs;
        int N_shells = int(R_max_shell/dR_shell)+1;
        shelled_nbs.resize(N_shells);

        for (tnit = nbs.begin(); tnit != nbs.end(); ++tnit) {
            PolarSeg* nb = (*tnit)->_nb;
            double R = votca::tools::abs((*sit1)->getPos()-nb->getPos()-(*tnit)->_L);            
            int shell_idx = int(R/dR_shell);
            shelled_nbs[shell_idx].push_back(*tnit);
        }
        
        // Sum over consecutive shells
        for (int sidx = 0; sidx < N_shells; ++sidx) {
            // Shell rms trackers
            double shell_sum = 0.0;
            double shell_term = 0.0;
            double shell_rms = 0.0;
            int shell_count = 0;
            // Interact with shell
            vector<TinyNeighbour*> &nb_shell = shelled_nbs[sidx];            
            double shell_R = (sidx+1)*dR_shell;            
            if (nb_shell.size() < 1) continue;            
            EWD::triple<double> ppuu(0,0,0);
            for (tnit = nb_shell.begin(); tnit < nb_shell.end(); ++tnit) {
                PolarSeg *nb = (*tnit)->_nb;
                vec L = (*tnit)->_L;
                nb->Translate(L);
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = nb->begin(); pit2 < nb->end(); ++pit2) {
    //                    ppuu = _ewdactor.U12_ERFC(*(*pit1), *(*pit2));
   //                     sum_ppuu += ppuu._pp;
     //                   shell_term = ppuu._pp;
      //                  shell_sum += shell_term;
       //                 shell_rms += shell_term*shell_term;
        //                shell_count += 1;
                        _actor.BiasIndu(*(*pit1), *(*pit2));
                        double ef = _actor.E_f(*(*pit1), *(*pit2));
                        sum_ppuu += ef;
                        shell_term = ef;
                        shell_sum += shell_term;
                        shell_rms += shell_term*shell_term;
                        shell_count += 1;
                    }
                }
                nb->Translate(-L);
            }
            shell_rms = sqrt(shell_rms/shell_count)*EWD::int2eV;
            sum += shell_sum;
            LOG(logDEBUG,*_log)
                << (format("  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %3$5d   ER = %2$+1.7f eV   dER2(sum) = %4$+1.3e eV") 
                % shell_R % (sum*EWD::int2eV) % nb_shell.size() % (shell_rms*shell_count) % (*sit1)->getId()).str() << flush;
            ofs
                << (format("MST DBG ...  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %3$5d   ER = %2$+1.7f eV   dER2(sum) = %4$+1.3e eV") 
                % shell_R % (sum*EWD::int2eV) % nb_shell.size() % (shell_rms*shell_count) % (*sit1)->getId()).str() << endl;
        }
    }
    
    ofs.close();
    return;
    
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    
//    for (int N = 0; N < 10; ++N) {
//        for (int na=-N; na<N+1; ++na) {
//        for (int nb=-N; nb<N+1; ++nb) {
//        for (int nc=-N; nc<N+1; ++nc) {
//            if ( std::abs(na) != N && std::abs(nb) != N && std::abs(nc) != N ) continue;
//            vec L = na*_a + nb*_b + nc*_c;
//            
//            for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
//                PolarSeg *pseg = *sit1;
//                
//                // Skip if N==0 and ID match
//                // Shift
//                // Interact
//                // Shift back
//                
//                
//            }
//        }}}
//    }
    
    
    
    
//    double ax = _a.getX();
//    double by = _b.getY();
//    double cz = _c.getZ();
//    vec ax_by_cz = vec(1./ax,1./by,1./cz);
//    double norm = 1./maxnorm(ax_by_cz);
//    LOG(logDEBUG,*_log) << ax << " " << by << " " << cz << " " << norm << flush;
    
    /*
    // FOR EACH FOREGROUND SEGMENT (FGC) ...
    for (sit1 = _fg_C.begin(); sit1 != _fg_C.end(); ++sit1) {

        // Bin midground into shells
        vector< vector<PolarSeg*> > shelled_mg_N;
        int N_shells = int(R_max_shell/dR_shell)+1;
        shelled_mg_N.resize(N_shells);

        for (sit2 = _mg_N.begin(); sit2 != _mg_N.end(); ++sit2) {
            double R = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
//            vec dr = (*sit1)->getPos()-(*sit2)->getPos();
//            dr = norm*vec(dr.getX()/ax, dr.getY()/by, dr.getZ()/cz);
//            double R = votca::tools::maxnorm(dr);
            
            int shell_idx = int(R/dR_shell);
            shelled_mg_N[shell_idx].push_back(*sit2);
        }

        // Sum over consecutive shells
        for (int sidx = 0; sidx < N_shells; ++sidx) {
            // Shell rms trackers
            double shell_sum = 0.0;
            double shell_term = 0.0;
            double shell_rms = 0.0;
            int shell_count = 0;
            // Interact with shell
            vector<PolarSeg*> &shell_mg = shelled_mg_N[sidx];            
            double shell_R = (sidx+1)*dR_shell;            
            if (shell_mg.size() < 1) continue;            
            EWD::triple<double> ppuu(0,0,0);
            for (sit2 = shell_mg.begin(); sit2 < shell_mg.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        _actor.BiasIndu(*(*pit1), *(*pit2));
                        double ef = _actor.E_f(*(*pit1), *(*pit2));
                        sum_ppuu += ef;
                        shell_term = ef;
                        shell_sum += shell_term;
                        shell_rms += shell_term*shell_term;
                        shell_count += 1;
                    }
                }
            }
            shell_rms = sqrt(shell_rms/shell_count)*EWD::int2eV;
            sum += shell_sum;
            LOG(logDEBUG,*_log)
                << (format("  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %3$5d   ER = %2$+1.7f eV   dER2(sum) = %4$+1.3e eV") 
                % shell_R % (sum*EWD::int2eV) % shell_mg.size() % (shell_rms*shell_count) % (*sit1)->getId()).str() << flush;
            ofs
                << (format("MST DBG ...  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %3$5d   ER = %2$+1.7f eV   dER2(sum) = %4$+1.3e eV") 
                % shell_R % (sum*EWD::int2eV) % shell_mg.size() % (shell_rms*shell_count) % (*sit1)->getId()).str() << endl;
        }
    }
    
    ofs.close();
    */
    return;
}


EWD::triple<> PEwald3D3D::ConvergeRealSpaceSum(vector<PolarSeg*> &target) {
    
    double sum = 0.0;
    double sum_pp = 0.0;
    double sum_pu = 0.0;
    double sum_uu = 0.0;
    _converged_R = false;
    
    LOG(logDEBUG,*_log) << flush 
        << "R-space energy via midground" << flush;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    vector<PolarNb*>::iterator nit;
    vector< vector<PolarSeg*> > ::iterator vsit;
    
    // ENERGY - REUSE NEIGHBOURS ?
    if (_did_field_pin_R_shell) {
        EWD::triple<double> ppuu(0,0,0);
        for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
            vector<PolarNb*> &nbs = (*sit1)->PolarNbs();
            for (nit = nbs.begin(); nit != nbs.end(); ++nit) {
                PolarSeg *nb = (*nit)->getNb();
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = nb->begin(); pit2 < nb->end(); ++pit2) {
                        ppuu = _ewdactor.U12_ERFC(*(*pit1), *(*pit2));
                        sum_pp += ppuu._pp;
                        sum_pu += ppuu._pu;
                        sum_uu += ppuu._uu;
//                        _actor.BiasIndu(*(*pit1), *(*pit2));
//                        double ef = _actor.E_f(*(*pit1), *(*pit2));
//                        sum_pp += ef;
                    }
                }
            }
            if (tools::globals::verbose) {
                LOG(logDEBUG,*_log)
                << (format("  o Id = %5$-4d Rc = %1$+02.7f   |MGN| = %2$5d   dF(rms) = %3$+1.3e V/m   [1eA => %4$+1.3e eV]") 
                % -1.0 % nbs.size() % -1.0  % -1.0 % ((*sit1)->getId())).str() << flush;
            }
        }
        _converged_R = true;
//        cout << endl << "XINTERACTOR " << _actor.getEPP() << " " << _actor.getEPU() << " "<< _actor.getEUU() << flush;
    }
    
    // ENERGY - REGENERATE NEIGHBOURS ?
    else {
        double dR_shell = 0.5;
        double R_overhead = 1.1;
        double R_add = 3;
        double R_max = _R_co*R_overhead+R_add;
        double R_max_shell = R_max+2*_polar_cutoff+_max_int_dist_qm0;
        this->SetupMidground(R_max);
        
        // FOR EACH FOREGROUND SEGMENT (FGC) ...
        unsigned energy_converged_count = 0;
        for (sit1 = target.begin(); sit1 != target.end(); ++sit1) {        
            (*sit1)->ClearPolarNbs();

            // Bin midground into shells
            vector< vector<PolarSeg*> > shelled_mg_N;
            int N_shells = int(R_max_shell/dR_shell)+1;
            shelled_mg_N.resize(N_shells);

            for (sit2 = _mg_N.begin(); sit2 != _mg_N.end(); ++sit2) {
                double R = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
                int shell_idx = int(R/dR_shell);
                shelled_mg_N[shell_idx].push_back(*sit2);
            }

            // Sum over consecutive shells
            for (int sidx = 0; sidx < N_shells; ++sidx) {
                // Shell rms trackers
                double shell_sum = 0.0;
                double shell_term = 0.0;
                double shell_rms = 0.0;
                int shell_count = 0;
                // Interact with shell
                vector<PolarSeg*> &shell_mg = shelled_mg_N[sidx];            
                double shell_R = (sidx+1)*dR_shell;            
                if (shell_mg.size() < 1) continue;            
                EWD::triple<double> ppuu(0,0,0);
                for (sit2 = shell_mg.begin(); sit2 < shell_mg.end(); ++sit2) {
                    for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                            ppuu = _ewdactor.U12_ERFC(*(*pit1), *(*pit2));
                            sum_pp += ppuu._pp;
                            sum_pu += ppuu._pu;
                            sum_uu += ppuu._uu;
                            shell_term = ppuu._pp + ppuu._pu + ppuu._uu;
                            shell_sum += shell_term;
                            shell_rms += shell_term*shell_term;
                            shell_count += 1;
//                            _actor.BiasIndu(*(*pit1), *(*pit2));
//                            double ef = _actor.E_f(*(*pit1), *(*pit2));
//                            sum_pp += ef;
//                            shell_term = ef;
//                            shell_sum += shell_term;
//                            shell_rms += shell_term*shell_term;
//                            shell_count += 1;
                        }
                    }
                }
                shell_rms = sqrt(shell_rms/shell_count)*EWD::int2eV;
                sum += shell_sum;
                if (tools::globals::verbose) {
                    LOG(logDEBUG,*_log)
                    << (format("  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %3$5d   ER = %2$+1.7f eV   dER2(sum) = %4$+1.3e eV") 
                    % shell_R % (sum*EWD::int2eV) % shell_mg.size() % (shell_rms*shell_count) % (*sit1)->getId()).str() << flush;
                }

                if (shell_rms*shell_count <= _crit_dE && shell_R >= _R_co) {
                    energy_converged_count += 1;
                    if (tools::globals::verbose){
                        LOG(logDEBUG,*_log)  
                        << (format("  :: ID = %2$-4d : Converged to precision as of Rc = %1$+1.3f nm") 
                        % shell_R % (*sit1)->getId()) << flush;
                    }
                    break;
                }
            }
        }
        
//        cout << endl << "XINTERACTOR " << _actor.getEPP() << " " << _actor.getEPU() << " "<< _actor.getEUU() << flush;

        if (energy_converged_count == target.size()) {
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision (%1$d items)") 
                % energy_converged_count) << flush;
            _converged_R = true;
        }
        else if (energy_converged_count < target.size()) {
            LOG(logERROR,*_log) << "ERROR Energy not converged on " 
                << target.size() - energy_converged_count << " counts." << flush;
            _converged_R = false;
        }
        else {
            assert(false);
        }
    }
    
    //boost::timer::auto_cpu_timer t0(*_log);
    //t0.start();
    //t0.stop();
    //t0.report();
    
    return EWD::triple<>(sum_pp, sum_pu, sum_uu);
}


EWD::triple<> PEwald3D3D::ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target) {
    
    // ATTENTION K-vectors are generated based on an interaction-energy
    //           criterion between FGC and BGP. Hence, the <target> density
    //           should at the very least be located within the space
    //           covered by FGC.
    if (!_did_generate_kvectors)
        this->GenerateKVectors(_fg_C, _bg_P);
    vector< EWD::KVector >::iterator kvit;
    
    double sum_re = 0.0;
    double sum_re_pp = 0.0;
    double sum_re_pu = 0.0;
    double sum_re_uu = 0.0;
    double sum_im = 0.0;
    _converged_K = false;
    
    // TWO COMPONENTS ZERO, ONE NON-ZERO
    LOG(logINFO,*_log) << flush 
        << "K-lines through origin: Checking K resonances" << flush;
    for (kvit = _kvecs_2_0.begin(); kvit < _kvecs_2_0.end(); ++kvit) {
        EWD::KVector kvec = *kvit;
        double Ak = _ewdactor.Ark2Expk2(kvec.getK());
        EWD::triple<EWD::cmplx> ppuu = _ewdactor.S1S2(kvec.getK(), target, _bg_P);        
        sum_re_pp += Ak*ppuu._pp._re;
        sum_re_pu += Ak*ppuu._pu._re;
        sum_re_uu += Ak*ppuu._uu._re;        
        sum_re += Ak*(ppuu._pp._re + ppuu._pu._re + ppuu._uu._re);        
        sum_im += Ak*(ppuu._pp._im + ppuu._pu._im + ppuu._uu._im);
    }
    
    LOG(logINFO,*_log)
        << (format("  :: RE %1$+1.7e IM %2$+1.7e")
            % (sum_re/_LxLyLz*EWD::int2eV)
            % (sum_im/_LxLyLz*EWD::int2eV)).str() << flush;
    
    // ONE COMPONENT ZERO, TWO NON-ZERO
    LOG(logINFO,*_log)
        << "K-planes through origin: Applying K resonances" << flush;
    
    double crit_grade = 1. * _kxyz_s1s2_norm;
    bool converged12 = false;
    kvit = _kvecs_1_0.begin();
    while (!converged12 && kvit < _kvecs_1_0.end()) {
        
        double de_this_shell = 0.0;
        int shell_count = 0;
        bool increment_grade = true;
        
        while (kvit < _kvecs_1_0.end()) {
            EWD::KVector kvec = *kvit;
            if (kvec.getGrade() < crit_grade) break;
            if (shell_count > 1000) {
                increment_grade = false;
                break;
            }
            EWD::triple<EWD::cmplx> ppuu = _ewdactor.AS1S2(kvec.getK(), target, _bg_P);
            
            sum_re_pp += ppuu._pp._re;
            sum_re_pu += ppuu._pu._re;
            sum_re_uu += ppuu._uu._re;            
            sum_re += ppuu._pp._re + ppuu._pu._re + ppuu._uu._re;
            sum_im += ppuu._pp._im + ppuu._pu._im + ppuu._uu._im;
            
            //de_this_shell += sqrt(as1s2._re*as1s2._re + as1s2._im*as1s2._im);
            de_this_shell += ppuu._pp._re + ppuu._pu._re + ppuu._uu._re;
            //cout << endl << std::showpos << std::scientific 
            //   << kvec.getX() << " " << kvec.getY() << " " << kvec.getZ() 
            //   << " grade " << kvec.getGrade() << " re " << (as1s2._re/_LxLyLz*_ewdactor.int2eV) << flush;            
            ++kvit;
            ++shell_count;
        }
        de_this_shell = (de_this_shell < 0.) ? -de_this_shell : de_this_shell;
        
        if (shell_count > 0){
            LOG(logDEBUG,*_log)
             << (format("  o M = %1$04d   G = %2$+1.3e   dE(rms) = %3$+1.3e eV")
             % shell_count
             % crit_grade
             % (de_this_shell/_LxLyLz*EWD::int2eV)).str() << flush;
        }
        
        if (shell_count > 10 && de_this_shell/_LxLyLz*EWD::int2eV < _crit_dE) {
            LOG(logINFO,*_log)
                << (format("  :: RE %1$+1.7e IM %2$+1.7e") 
                % (sum_re/_LxLyLz*EWD::int2eV)
                % (sum_im/_LxLyLz*EWD::int2eV)).str() << flush;
            converged12 = true;
        }
        
        if (increment_grade) crit_grade /= 10.0;
    }
    
    
    // ZERO COMPONENTS ZERO, THREE NON-ZERO
    LOG(logINFO,*_log)
        << "K-space (off-axis): Applying K resonances" << flush;    
    
    crit_grade = 1. * _kxyz_s1s2_norm;
    double converged03 = false;
    kvit = _kvecs_0_0.begin();
    while (!converged03 && kvit < _kvecs_0_0.end()) {
        
        double de_this_shell = 0.0;
        int shell_count = 0;
        bool increment_grade = true;
        
        while (kvit < _kvecs_0_0.end()) {
            EWD::KVector kvec = *kvit;
            if (kvec.getGrade() < crit_grade) break;
            if (shell_count > 1000) {
                increment_grade = false;
                break;
            }
            EWD::triple<EWD::cmplx> ppuu = _ewdactor.AS1S2(kvec.getK(), target, _bg_P);
            
            sum_re_pp += ppuu._pp._re;
            sum_re_pu += ppuu._pu._re;
            sum_re_uu += ppuu._uu._re;            
            sum_re += ppuu._pp._re + ppuu._pu._re + ppuu._uu._re;
            sum_im += ppuu._pp._im + ppuu._pu._im + ppuu._uu._im;

            //de_this_shell += sqrt(as1s2._re*as1s2._re + as1s2._im*as1s2._im);
            de_this_shell += ppuu._pp._re + ppuu._pu._re + ppuu._uu._re;
            //cout << endl << std::showpos << std::scientific 
            //   << kvec.getX() << " " << kvec.getY() << " " << kvec.getZ() 
            //   << " grade " << kvec.getGrade() << " re " << (as1s2._re/_LxLyLz*_ewdactor.int2eV) << flush;            
            ++kvit;
            ++shell_count;
        }
        de_this_shell = (de_this_shell < 0.) ? -de_this_shell : de_this_shell;
        
        if (shell_count > 0) {
            LOG(logDEBUG,*_log)
             << (format("  o M = %1$04d   G = %2$+1.3e   dE(rms) = %3$+1.3e eV")
             % shell_count
             % crit_grade
             % (de_this_shell/_LxLyLz*EWD::int2eV)).str() << flush;
        }
        if (shell_count > 10 && de_this_shell/_LxLyLz*EWD::int2eV < _crit_dE) {
            LOG(logINFO,*_log)
                << (format("  :: RE %1$+1.7e IM %2$+1.7e") 
                % (sum_re/_LxLyLz*EWD::int2eV)
                % (sum_im/_LxLyLz*EWD::int2eV)).str() << flush;
            converged03 = true;
        }
        
        if (increment_grade) crit_grade /= 10.0;
    }
    
    _converged_K = converged12 && converged03;
    
    if (_converged_K)
        LOG(logINFO,*_log)
            << (format(":::: Converged to precision, {0-2}, {1-2}, {0-3}."))
            << flush;
    else ;

    return EWD::triple<>(sum_re_pp/_LxLyLz, sum_re_pu/_LxLyLz, sum_re_uu/_LxLyLz);
}


EWD::triple<> PEwald3D3D::CalculateShapeCorrection(vector<PolarSeg*> &target) {
    
    LOG(logDEBUG,*_log) << flush
        << "Energy correction terms" << flush;
    LOG(logDEBUG,*_log)
        << "  o Shape-correction to energy, using '" << _shape << "'" << flush;
    
    EWD::triple<double> ppuu = _ewdactor.U12_ShapeTerm(target, _bg_P,
        _shape, _LxLyLz, _log);
    double sum_pp = ppuu._pp;
    double sum_pu = ppuu._pu;
    double sum_uu = ppuu._uu;
    
//    vector<PolarSeg*>::iterator sit1; 
//    vector<APolarSite*> ::iterator pit1;
//    vector<PolarSeg*>::iterator sit2; 
//    vector<APolarSite*> ::iterator pit2;
//    
//    double EJ = 0.0;
    
//    if (_shape == "xyslab") {
//        EWD::triple<double> ppuu(0,0,0);
//        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
//           for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
//              for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//                 for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
//                    ppuu = _ewdactor.U12_XYSlab(*(*pit1), *(*pit2));
//                    sum_pp += ppuu._pp;
//                    sum_pu += ppuu._pu;
//                    sum_uu += ppuu._uu;
//                 }
//              }
//           }
//        }
//        sum_pp *= -2*M_PI/_LxLyLz;
//        sum_pu *= -2*M_PI/_LxLyLz;
//        sum_uu *= -2*M_PI/_LxLyLz;
//        EJ = sum_pp + sum_pu + sum_uu;
//        
//        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//        EWD::triple<double> ppuu_fast = _ewdactor.U12_XYSlab_Factorized(_fg_C, _bg_P);
//        double sum_pp_fast = ppuu_fast._pp;
//        double sum_pu_fast = ppuu_fast._pu;
//        double sum_uu_fast = ppuu_fast._uu;
//        sum_pp_fast *= -4*M_PI/_LxLyLz;
//        sum_pu_fast *= -4*M_PI/_LxLyLz;
//        sum_uu_fast *= -4*M_PI/_LxLyLz;
//        cout << endl << "Slow: " << sum_pp << " " << sum_pu << " " << sum_uu << flush;
//        cout << endl << "Fast: " << sum_pp_fast << " " << sum_pu_fast << " " << sum_uu_fast << flush;
//        
//        EWD::triple<double> ppuu_cube = _ewdactor.U12_ShapeTerm_Factorized(_fg_C, _bg_P, _shape, _LxLyLz);
//        double sum_pp_cube = ppuu_cube._pp;
//        double sum_pu_cube = ppuu_cube._pu;
//        double sum_uu_cube = ppuu_cube._uu;
//        cout << endl << "Cube: " << sum_pp_cube << " " << sum_pu_cube << " " << sum_uu_cube << flush;
//        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//    }
//    else {
//        LOG(logERROR,*_log)
//            << (format("Shape %1$s not implemented. Setting EJ = 0.0 ...") 
//            % _shape) << flush;
//        EJ = 0.0;
//    }
    
    return EWD::triple<>(sum_pp, sum_pu, sum_uu);
    //return EJ;
}


EWD::triple<> PEwald3D3D::CalculateForegroundCorrection(vector<PolarSeg*> &target) {
    
    LOG(logDEBUG,*_log)
        << "  o Foreground-correction to energy via FGN" << flush;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    //double EC = 0.0;
    double sum_pp = 0.0;
    double sum_pu = 0.0;
    double sum_uu = 0.0;
    
    EWD::triple<double> ppuu(0,0,0);
    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    ppuu = _ewdactor.U12_ERF(*(*pit1), *(*pit2));
                    sum_pp += ppuu._pp;
                    sum_pu += ppuu._pu;
                    sum_uu += ppuu._uu;
                }
            }
        }
    }
    
    //EC = sum_pp + sum_pu + sum_uu;    
    return EWD::triple<>(sum_pp, sum_pu, sum_uu);
    //return EC;
}


void PEwald3D3D::Field_ConvergeRealSpaceSum() {
    
    //double sum = 0.0;
    _field_converged_R = false;
    
    LOG(logDEBUG,*_log) << flush 
        << "R-space fields via midground" << flush;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    // GENERATE MIDGROUND & ASSEMBLE IT INTO SHELLS
    double dR_shell = 0.5;
    double R_overhead = 1.1;
    double R_add = 3;
    double R_max = _R_co*R_overhead+R_add;
    double R_max_shell = R_max+2*_polar_cutoff+_max_int_dist_qm0;
    this->SetupMidground(R_max);

    
    // FOR EACH FOREGROUND SEGMENT (FGC) ...
    unsigned field_converged_count = 0;
    for (sit1 = _fg_C.begin(); sit1 != _fg_C.end(); ++sit1) {
        (*sit1)->ClearPolarNbs();
        
        // Bin midground into shells
        vector< vector<PolarSeg*> > shelled_mg_N;
        int N_shells = int(R_max_shell/dR_shell)+1;
        shelled_mg_N.resize(N_shells);
        
        for (sit2 = _mg_N.begin(); sit2 != _mg_N.end(); ++sit2) {
            double R = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
            int shell_idx = int(R/dR_shell);
            shelled_mg_N[shell_idx].push_back(*sit2);
        }        
        
        // Sum over consecutive shells
        for (int sidx = 0; sidx < N_shells; ++sidx) {
            // Shell rms trackers
            double shell_rms = 0.0;
            int shell_count = 0;
            // Interact with shell
            vector<PolarSeg*> &shell_mg = shelled_mg_N[sidx];
            double shell_R = (sidx+1)*dR_shell;            
            if (shell_mg.size() < 1) continue;            
            for (sit2 = shell_mg.begin(); sit2 != shell_mg.end(); ++sit2) {
                if (_save_nblist) (*sit1)->AddNewPolarNb(*sit2);
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        shell_rms += _ewdactor.FPU12_ERFC_At_By(*(*pit1), *(*pit2));
                        shell_count += 1;
//                        vec f0 = (*pit1)->getFieldP();
//                        _actor.BiasStat(*(*pit1), *(*pit2));
//                        _actor.FieldPerm(*(*pit1), *(*pit2));
//                        vec f1 = (*pit1)->getFieldP();
//                        vec df = f1-f0;
//                        shell_rms += df*df;
//                        shell_count += 1;
                    }
                }
            }
            
            // Assert convergence: Energy of dipole of size 0.1*e*nm summed over shell
            shell_rms = sqrt(shell_rms/shell_count)*EWD::int2V_m;
            double e_measure = shell_rms*1e-10*shell_count; // 
        
            if (tools::globals::verbose){
                LOG(logDEBUG,*_log)
                << (format("  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %2$5d   dF(rms) = %3$+1.3e V/m   [1eA => %4$+1.3e eV]") 
                % shell_R % shell_mg.size() % shell_rms  % e_measure % ((*sit1)->getId())).str() << flush;
            }
            
            if (e_measure <= _crit_dE && shell_R >= _R_co) {
                field_converged_count += 1;
                if (tools::globals::verbose){
                    LOG(logDEBUG,*_log)
                    << (format("  :: ID = %2$-4d Converged to precision as of Rc = %1$+1.3f nm") 
                    % shell_R % (*sit1)->getId()) << flush;
                }
                break;
            }
        }
    }
    
    if (field_converged_count == _fg_C.size()) {
        LOG(logDEBUG,*_log)  
            << (format(":::: Converged to precision (%1$d items)") 
            % field_converged_count) << flush;
        _field_converged_R = true;
    }
    else if (field_converged_count < _fg_C.size()) {
        LOG(logERROR,*_log) << "ERROR Field not converged on " 
            << _fg_C.size() - field_converged_count << " counts." << flush;
        _field_converged_R = false;
    }
    else {
        assert(false);
    }
   	
	if (_save_nblist) 
        _did_field_pin_R_shell = true;
    else
        _did_field_pin_R_shell = false;
    
//    // THERE IS A MEMORY ISSUE HERE - VERY STRANGE
//    // Change 18000 to 20000 and the leak disappears!?
//    double sum = 0.0;
//    _field_converged_R = false;
//    LOG(logDEBUG,*_log) << flush 
//        << "R-space fields via midground" << flush;
//    
//    vector<PolarSeg*>::iterator sit1;
//    
//    // FOR EACH FOREGROUND SEGMENT (FGC) ...
//    int field_converged_count = 0;
//    for (sit1 = _fg_C.begin(); sit1 != _fg_C.end(); ++sit1) {
//        
//        for (int cnt = 0; cnt < 18000; ++cnt) {
//            PolarNb *new_nb = (*sit1)->AddNewPolarNb(*sit1);
//        }
//        
//        (*sit1)->ClearPolarNbs();
//        assert((*sit1)->PolarNbs().size() == 0);
//        cout << endl << "sit1 = ID " << (*sit1)->getId() << flush;
//        
//        
//    }
//    _did_field_pin_R_shell = true;    
    
    return;
}


void PEwald3D3D::Field_ConvergeReciprocalSpaceSum() {

    this->GenerateKVectors(_fg_C, _bg_P);
    double sum_re = 0.0;
    double sum_im = 0.0;
    _field_converged_K = false;
    double rV = 1./_LxLyLz;
    
    vector< EWD::KVector >::iterator kvit;
    
    // TWO COMPONENTS ZERO, ONE NON-ZERO
    LOG(logINFO,*_log) << flush 
        << "K-lines through origin: Checking K resonances" << flush;
    for (kvit = _kvecs_2_0.begin(); kvit < _kvecs_2_0.end(); ++kvit) {
        EWD::KVector kvec = *kvit;
        EWD::cmplx f_as1s2 = _ewdactor.FPU12_AS1S2_At_By(kvec.getK(), _fg_C, _bg_P, rV);
        sum_re += sqrt(f_as1s2._re);
        sum_im += f_as1s2._im;
    }
    
    LOG(logINFO,*_log)
        << (format("  :: RE %1$+1.7e IM %2$+1.7e")
            % (sum_re*EWD::int2V_m)
            % (sum_im*EWD::int2V_m)).str() << flush;
    
    // ONE COMPONENT ZERO, TWO NON-ZERO
    LOG(logINFO,*_log)
        << "K-planes through origin: Applying K resonances" << flush;    
    
    double crit_grade = 1. * _kxyz_s1s2_norm;
    bool converged12 = false;
    kvit = _kvecs_1_0.begin();
    while (!converged12 && kvit < _kvecs_1_0.end()) {
        
        double shell_rms = 0.0;
        int rms_count = 0;
        
        while (kvit < _kvecs_1_0.end()) {
            EWD::KVector kvec = *kvit;
            if (kvec.getGrade() < crit_grade) break;
            EWD::cmplx f_as1s2 = _ewdactor.FPU12_AS1S2_At_By(kvec.getK(), _fg_C, _bg_P, rV);
            sum_re += f_as1s2._re;
            sum_im += f_as1s2._im;
            shell_rms += f_as1s2._re;
            //cout << endl << std::showpos << std::scientific 
            //   << kvec.getX() << " " << kvec.getY() << " " << kvec.getZ() 
            //   << " grade " << kvec.getGrade() << " re " << (as1s2._re/_LxLyLz*_ewdactor.int2eV) << flush;            
            ++kvit;
            ++rms_count;
        }
        shell_rms = (rms_count > 0) ? sqrt(shell_rms/rms_count)*EWD::int2V_m : 0.0;
        double e_measure = shell_rms*1e-10*rms_count;
        
        if (rms_count > 0) {
            LOG(logDEBUG,*_log)
             << (format("  o M = %1$04d   G = %2$+1.3e   dF(rms) = %3$+1.3e V/m   [1eA => %4$+1.3e eV]")
             % rms_count
             % crit_grade
             % shell_rms
             % e_measure).str() << flush;
        }
        
        if (rms_count > 10 && e_measure <= _crit_dE) {
            LOG(logINFO,*_log)
                << (format("  :: RE %1$+1.7e IM %2$+1.7e") 
                % (sqrt(sum_re)*EWD::int2V_m)
                % (sum_im*EWD::int2V_m)).str() << flush;
            converged12 = true;
        }
        
        crit_grade /= 10.0;
    }
    
    // ZERO COMPONENTS ZERO, THREE NON-ZERO
    LOG(logINFO,*_log)
        << "K-space (off-axis): Applying K resonances" << flush;
    
    crit_grade = 1. * _kxyz_s1s2_norm;
    double converged03 = false;
    kvit = _kvecs_0_0.begin();
    while (!converged03 && kvit < _kvecs_0_0.end()) {
        
        double shell_rms = 0.0;
        int rms_count = 0;
        
        while (kvit < _kvecs_0_0.end()) {
            EWD::KVector kvec = *kvit;
            if (kvec.getGrade() < crit_grade) break;
            EWD::cmplx f_as1s2 = _ewdactor.FPU12_AS1S2_At_By(kvec.getK(), _fg_C, _bg_P, rV);
            sum_re += f_as1s2._re;
            sum_im += f_as1s2._im;
            shell_rms += f_as1s2._re;
            //cout << endl << std::showpos << std::scientific 
            //   << kvec.getX() << " " << kvec.getY() << " " << kvec.getZ() 
            //   << " grade " << kvec.getGrade() << " re " << (as1s2._re/_LxLyLz*_ewdactor.int2eV) << flush;            
            ++kvit;
            ++rms_count;
        }
        shell_rms = (rms_count > 0) ? sqrt(shell_rms/rms_count)*EWD::int2V_m : 0.0;
        double e_measure = shell_rms*1e-10*rms_count;
        
        if (rms_count > 0) {
            LOG(logDEBUG,*_log)
             << (format("  o M = %1$04d   G = %2$+1.3e   dF(rms) = %3$+1.3e V/m   [1eA => %4$+1.3e eV]")
             % rms_count
             % crit_grade
             % shell_rms
             % e_measure).str() << flush;
        }
        
        if (rms_count > 10 && e_measure <= _crit_dE) {
            LOG(logINFO,*_log)
                << (format("  :: RE %1$+1.7e IM %2$+1.7e") 
                % (sqrt(sum_re)*EWD::int2V_m)
                % (sum_im*EWD::int2V_m)).str() << flush;
            converged03 = true;
        }
        
        crit_grade /= 10.0;
    }
    
    _field_converged_K = converged12 && converged03;
    
    if (_field_converged_K) {
        LOG(logINFO,*_log)
            << (format(":::: Converged to precision, {0-2}, {1-2}, {0-3}."))
            << flush;
    }
    else ;

    return;
}


void PEwald3D3D::Field_CalculateForegroundCorrection() {
    
    LOG(logDEBUG,*_log)
        << "  o Foreground-correction to fields via FGN" << flush;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    double rms = 0.0;
    int rms_count = 0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    rms += _ewdactor.FPU12_ERF_At_By(*(*pit1), *(*pit2));
                    rms_count += 1;
                }
            }
        }
    }
    rms = sqrt(rms/rms_count)*EWD::int2V_m;
    
    return;
}


void PEwald3D3D::Field_CalculateShapeCorrection() {

    LOG(logDEBUG,*_log) << flush
        << "Field correction terms" << flush;
    LOG(logDEBUG,*_log)
        << "  o Shape-correction to fields, using '" << _shape << "'" << flush;
    
    _ewdactor.FPU12_ShapeField_At_By(_fg_C, _bg_P, _shape, _LxLyLz);
    
//    vector<PolarSeg*>::iterator sit1; 
//    vector<APolarSite*> ::iterator pit1;
//    vector<PolarSeg*>::iterator sit2; 
//    vector<APolarSite*> ::iterator pit2;
//    
//    double rms = 0.0;
//    int rms_count = 0;    
//    if (_shape == "xyslab") {
//        double TwoPi_V = 2*M_PI/_LxLyLz;
//        
////        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
////           for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
////              for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
////                 for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
////                    rms += _ewdactor.F12_XYSlab_At_By(*(*pit1), *(*pit2), TwoPi_V);
////                    rms_count += 1;
////                 }
////              }
////           }
////        }
//        
//        _ewdactor.FPU12_XYSlab_ShapeField_At_By(_fg_C, _bg_P, TwoPi_V);
//        
//    }
//    else {
//        LOG(logERROR,*_log)
//            << (format("Shape %1$s not implemented. Setting EJ = 0.0 ...") 
//            % _shape) << flush;
//    }
//    rms = sqrt(rms/rms_count)*EWD::int2V_m;
    
    return;
}


void PEwald3D3D::Potential_ConvergeRealSpaceSum(vector<PolarSeg*> &target) {
    
    double sum = 0.0;
    double sum_phi = 0.0;
    _potential_converged_R = false;
    
    LOG(logDEBUG,*_log) << flush 
        << "R-space potentials via midground" << flush;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    vector<PolarNb*>::iterator nit;
    vector< vector<PolarSeg*> > ::iterator vsit;
    
    bool neighbours_stored = false;
    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        int nb_count = (*sit1)->PolarNbs().size();
        if (nb_count > 0) {
            neighbours_stored = true;
            break;
        }
    }
    
    // ENERGY - REUSE NEIGHBOURS ?
    if (neighbours_stored) {
        for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
            vector<PolarNb*> &nbs = (*sit1)->PolarNbs();
            for (nit = nbs.begin(); nit != nbs.end(); ++nit) {
                PolarSeg *nb = (*nit)->getNb();
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = nb->begin(); pit2 < nb->end(); ++pit2) {
                        double phi = _ewdactor.PhiPU12_ERFC_At_By(*(*pit1), *(*pit2));
                        sum_phi += phi;
                    }
                }
            }
            if (tools::globals::verbose) {
                LOG(logDEBUG,*_log)
                << (format("  o Id = %5$-4d Rc = %1$+02.7f   |MGN| = %2$5d   dF(rms) = %3$+1.3e V/m   [1eA => %4$+1.3e eV]") 
                % -1.0 % nbs.size() % -1.0  % -1.0 % ((*sit1)->getId())).str() << flush;
            }
        }
        _potential_converged_R = true;
    }
    
    // ENERGY - REGENERATE NEIGHBOURS ?
    else {
        double dR_shell = 0.5;
        double R_overhead = 1.1;
        double R_add = 3;
        double R_max = _R_co*R_overhead+R_add;
        double R_max_shell = R_max+2*_polar_cutoff+_max_int_dist_qm0;
        this->SetupMidground(R_max);
        
        // FOR EACH FOREGROUND SEGMENT (FGC) ...
        unsigned energy_converged_count = 0;
        for (sit1 = target.begin(); sit1 != target.end(); ++sit1) {        
            (*sit1)->ClearPolarNbs();

            // Bin midground into shells
            vector< vector<PolarSeg*> > shelled_mg_N;
            int N_shells = int(R_max_shell/dR_shell)+1;
            shelled_mg_N.resize(N_shells);

            for (sit2 = _mg_N.begin(); sit2 != _mg_N.end(); ++sit2) {
                double R = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
                int shell_idx = int(R/dR_shell);
                shelled_mg_N[shell_idx].push_back(*sit2);
            }

            // Sum over consecutive shells
            for (int sidx = 0; sidx < N_shells; ++sidx) {
                // Shell rms trackers
                double shell_sum = 0.0;
                double shell_term = 0.0;
                double shell_rms = 0.0;
                int shell_count = 0;
                // Interact with shell
                vector<PolarSeg*> &shell_mg = shelled_mg_N[sidx];            
                double shell_R = (sidx+1)*dR_shell;            
                if (shell_mg.size() < 1) continue;            
                EWD::triple<double> ppuu(0,0,0);
                for (sit2 = shell_mg.begin(); sit2 < shell_mg.end(); ++sit2) {
                    for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                            double phi = _ewdactor.PhiPU12_ERFC_At_By(*(*pit1), *(*pit2));
                            sum_phi += phi;
                            shell_term = phi;
                            shell_sum += shell_term;
                            shell_rms += shell_term*shell_term;
                            shell_count += 1;
                        }
                    }
                }
                shell_rms = sqrt(shell_rms/shell_count)*EWD::int2eV;
                sum += shell_sum;
                if (tools::globals::verbose) {
                    LOG(logDEBUG,*_log)
                    << (format("  o ID = %5$-4d Rc = %1$+02.7f   |MGN| = %3$5d   ER = %2$+1.7f V   dER2(sum) = %4$+1.3e V") 
                    % shell_R % (sum*EWD::int2eV) % shell_mg.size() % (shell_rms*shell_count) % (*sit1)->getId()).str() << flush;
                }

                if (shell_rms*shell_count <= _crit_dE && shell_R >= _R_co) {
                    energy_converged_count += 1;
                    if (tools::globals::verbose) {
                        LOG(logDEBUG,*_log)  
                        << (format("  :: ID = %2$-4d : Converged to precision as of Rc = %1$+1.3f nm") 
                        % shell_R % (*sit1)->getId()) << flush;
                    }
                    break;
                }
            }
        }
        
        if (energy_converged_count == target.size()) {
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision (%1$d items)") 
                % energy_converged_count) << flush;
            _potential_converged_R = true;
        }
        else if (energy_converged_count < target.size()) {
            LOG(logERROR,*_log) << "ERROR Energy not converged on " 
                << target.size() - energy_converged_count << " counts." << flush;
            _potential_converged_R = false;
        }
        else {
            assert(false);
        }
    }
    
    //boost::timer::auto_cpu_timer t0(*_log);
    //t0.start();
    //t0.stop();
    //t0.report();
    
    return;
}


void PEwald3D3D::Potential_ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target) {
    
    // ATTENTION K-vectors are generated based on an interaction-energy
    //           criterion between FGC and BGP. Hence, the <target> density
    //           should at the very least be located within the space
    //           covered by FGC.
    if (!_did_generate_kvectors)
        this->GenerateKVectors(_fg_C, _bg_P);
    vector< EWD::KVector >::iterator kvit;
    
    double sum_re = 0.0;
    double sum_im = 0.0;
    _potential_converged_K = false;
    double rV = 1./_LxLyLz;
    
    // TWO COMPONENTS ZERO, ONE NON-ZERO
    LOG(logINFO,*_log) << flush 
        << "K-lines through origin: Checking K resonances" << flush;
    for (kvit = _kvecs_2_0.begin(); kvit < _kvecs_2_0.end(); ++kvit) {
        EWD::KVector kvec = *kvit;
        EWD::cmplx f_as1s2 = _ewdactor.PhiPU12_AS1S2_At_By(kvec.getK(), target, _bg_P, rV);
        sum_re += sqrt(f_as1s2._re);
        sum_im += f_as1s2._im;
    }
    
    LOG(logINFO,*_log)
        << (format("  :: RE %1$+1.7e IM %2$+1.7e")
            % (sum_re*EWD::int2eV)
            % (sum_im*EWD::int2eV)).str() << flush;
    
    // ONE COMPONENT ZERO, TWO NON-ZERO
    LOG(logINFO,*_log)
        << "K-planes through origin: Applying K resonances" << flush;    
    
    double crit_grade = 1. * _kxyz_s1s2_norm;
    bool converged12 = false;
    kvit = _kvecs_1_0.begin();
    while (!converged12 && kvit < _kvecs_1_0.end()) {
        
        double shell_rms = 0.0;
        int rms_count = 0;
        
        while (kvit < _kvecs_1_0.end()) {
            EWD::KVector kvec = *kvit;
            if (kvec.getGrade() < crit_grade) break;
            EWD::cmplx f_as1s2 = _ewdactor.PhiPU12_AS1S2_At_By(kvec.getK(), target, _bg_P, rV);
            sum_re += f_as1s2._re;
            sum_im += f_as1s2._im;
            shell_rms += f_as1s2._re;
            ++kvit;
            ++rms_count;
        }
        shell_rms = (rms_count > 0) ? sqrt(shell_rms/rms_count)*EWD::int2eV : 0.0;
        double e_measure = shell_rms*rms_count;
        
        if (rms_count > 0) {
            LOG(logDEBUG,*_log)
             << (format("  o M = %1$04d   G = %2$+1.3e   dPhi(rms) = %3$+1.3e V   [1e => %4$+1.3e eV]")
             % rms_count
             % crit_grade
             % shell_rms
             % e_measure).str() << flush;
        }
        
        if (rms_count > 10 && e_measure <= _crit_dE) {
            LOG(logINFO,*_log)
                << (format("  :: RE %1$+1.7e IM %2$+1.7e") 
                % (sqrt(sum_re)*EWD::int2eV)
                % (sum_im*EWD::int2eV)).str() << flush;
            converged12 = true;
        }
        
        crit_grade /= 10.0;
    }
    
    // ZERO COMPONENTS ZERO, THREE NON-ZERO
    LOG(logINFO,*_log)
        << "K-space (off-axis): Applying K resonances" << flush;
    
    crit_grade = 1. * _kxyz_s1s2_norm;
    double converged03 = false;
    kvit = _kvecs_0_0.begin();
    while (!converged03 && kvit < _kvecs_0_0.end()) {
        
        double shell_rms = 0.0;
        int rms_count = 0;
        
        while (kvit < _kvecs_0_0.end()) {
            EWD::KVector kvec = *kvit;
            if (kvec.getGrade() < crit_grade) break;
            EWD::cmplx f_as1s2 = _ewdactor.PhiPU12_AS1S2_At_By(kvec.getK(), target, _bg_P, rV);
            sum_re += f_as1s2._re;
            sum_im += f_as1s2._im;
            shell_rms += f_as1s2._re;
            ++kvit;
            ++rms_count;
        }
        shell_rms = (rms_count > 0) ? sqrt(shell_rms/rms_count)*EWD::int2eV : 0.0;
        double e_measure = shell_rms*1e-10*rms_count;
        
        if (rms_count > 0) {
            LOG(logDEBUG,*_log)
             << (format("  o M = %1$04d   G = %2$+1.3e   dF(rms) = %3$+1.3e V   [1e => %4$+1.3e eV]")
             % rms_count
             % crit_grade
             % shell_rms
             % e_measure).str() << flush;
        }
        
        if (rms_count > 10 && e_measure <= _crit_dE) {
            LOG(logINFO,*_log)
                << (format("  :: RE %1$+1.7e IM %2$+1.7e") 
                % (sqrt(sum_re)*EWD::int2eV)
                % (sum_im*EWD::int2eV)).str() << flush;
            converged03 = true;
        }
        
        crit_grade /= 10.0;
    }
    
    _potential_converged_K = converged12 && converged03;
    
    if (_potential_converged_K)
        LOG(logINFO,*_log)
            << (format(":::: Converged to precision, {0-2}, {1-2}, {0-3}."))
            << flush;
    else ;

    return;
}


void PEwald3D3D::Potential_CalculateForegroundCorrection(vector<PolarSeg*> &target) {
    LOG(logDEBUG,*_log)
        << "  o Foreground-correction to potentials via FGN" << flush;
    
    vector<PolarSeg*>::iterator sit1;
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2;
    vector<APolarSite*> ::iterator pit2;
    
    double rms = 0.0;
    int rms_count = 0;
    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    rms += _ewdactor.PhiPU12_ERF_At_By(*(*pit1), *(*pit2));
                    rms_count += 1;
                }
            }
        }
    }
    rms = sqrt(rms/rms_count)*EWD::int2eV;
    
    return;
}


void PEwald3D3D::Potential_CalculateShapeCorrection(vector<PolarSeg*> &target) {
    LOG(logDEBUG,*_log) << flush
        << "Potential correction terms" << flush;
    LOG(logDEBUG,*_log)
        << "  o Shape-correction to potentials, using '" << _shape << "'" << flush;
    
    _ewdactor.PhiPU12_ShapeField_At_By(target, _bg_P, _shape, _LxLyLz);
    return;
}


}}
