#include <votca/xtp/xinductor.h>
#include <boost/format.hpp>
#include <vector>
#include <boost/timer/timer.hpp>

using boost::format;

namespace votca { namespace xtp {

    
XInductor::~XInductor() {
    vector<InduWorker*>::iterator wit;
    for (wit = _indus.begin(); wit != _indus.end(); ++wit) {
        delete *wit;
    }
    _indus.clear();
}
    
    
void XInductor::Evaluate(XJob *job) {    
    
    this->Configure(job);
    
//    _qm0.clear();
//    _mm1.clear();
//    _mm2.clear();
//    _qmm.clear();
//    
//    _job = job;
//    _isConverged = !this->_induce;
//    
//    // QMM = [ QM0 --- MM1 ] --- MM2 = [ MM2 ]    
//    _qm0 = job->getPolarTop()->QM0();
//    _mm1 = job->getPolarTop()->MM1();
//    _mm2 = job->getPolarTop()->MM2();    
//
//    _qmm.reserve(_qm0.size()+_mm1.size());
//    _qmm.insert(_qmm.end(),_qm0.begin(),_qm0.end());
//    _qmm.insert(_qmm.end(),_mm1.begin(),_mm1.end());

    // ++++++++++++++++++++++++++ //
    // (De-)polarize, charge to N //
    // ++++++++++++++++++++++++++ //

    if (job->StartFromCPT()) {
        // Permanent fields already computed, do not zero these out ...
        LOG(logDEBUG,*_log) << "Carry out partial depolarization from CPT." << flush;
        vector< PolarSeg* >   ::iterator sit;
        vector< APolarSite* > ::iterator pit;

        // Partially depolarize inner sphere
        for (sit = _qmm.begin(); sit < _qmm.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->ResetFieldU();
            (*pit)->ResetU1Hist();
        }}

        // Partially depolarize outer shell
        for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->ResetFieldU();
            (*pit)->ResetU1Hist();
        }}
    }
    else {        
        LOG(logDEBUG,*_log) << "Carry out full depolarization." << flush;        
        vector< PolarSeg* >   ::iterator sit;
        vector< APolarSite* > ::iterator pit;

        // Depolarize inner sphere
        for (sit = _qmm.begin(); sit < _qmm.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }}

        // Depolarize outer shell
        for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->Depolarize();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }}
    }
    
//    // +++++++++++++++++ //
//    // Induction workers //
//    // +++++++++++++++++ //




//
//    for (int id = 0; id < this->_subthreads; ++id) {
//        InduWorker *newIndu = new InduWorker(id, _top, this);
//        _indus.push_back(newIndu);
//        newIndu->InitSpheres(&_qmm, &_mm2);
//        newIndu->SetSwitch(1);
//    }
//    
//    this->InitChunks();
    
    // ++++++++++++++++++++++++++ //
    // Compute state energy       //
    // ++++++++++++++++++++++++++ //

    //double  E_state  = 0.0;
    int     iter     = 0;
    
    boost::timer::cpu_timer cpu_t;
    cpu_t.start();
    boost::timer::cpu_times t0 = cpu_t.elapsed();    
    if (this->_induce) iter      = this->Induce(job);
    boost::timer::cpu_times t1 = cpu_t.elapsed();
    //if (this->_induce) E_state   = this->Energy(job);
    //else               E_state   = this->EnergyStatic(job);
    if (this->_induce) this->Energy(job);
    else                this->EnergyStatic(job);
    boost::timer::cpu_times t2 = cpu_t.elapsed();
    
    double t_indu = (t1.wall - t0.wall)/1e9/60.;
    double t_ener = (t2.wall - t1.wall)/1e9/60.;
    LOG(logINFO,*_log) << (format("  o Total:     %1$1.2f min")
        % (t_ener+t_indu)) << flush;
    LOG(logINFO,*_log) << (format("  o Induction: %1$1.2f min")
        % (t_indu)) << flush;
    LOG(logINFO,*_log) << (format("  o Energy:    %1$1.2f min")
        % (t_ener)) << flush;

    job->setInduIter(iter);
    return;
}


void XInductor::Configure(XJob *job) {
    // SETUP MULTIPOLAR DENSITIES
    _qm0.clear();
    _mm1.clear();
    _mm2.clear();
    _qmm.clear();
    
    _job = job;
    _isConverged = !this->_induce;
    
    // QMM = [ QM0 --- MM1 ] --- MM2 = [ MM2 ]    
    _qm0 = job->getPolarTop()->QM0();
    _mm1 = job->getPolarTop()->MM1();
    _mm2 = job->getPolarTop()->MM2();    

    
    //_qm0[0]->WriteMPS("crap.mps");
    
    
    _qmm.reserve(_qm0.size()+_mm1.size());
    _qmm.insert(_qmm.end(),_qm0.begin(),_qm0.end());
    _qmm.insert(_qmm.end(),_mm1.begin(),_mm1.end());
    
    //_qmm[0]->WriteMPS("crap.mps");
    
    
    // INDUCTION & ENERGY WORKERS
    // Delete previous ...
    vector<InduWorker*>::iterator wit;
    for (wit = _indus.begin(); wit != _indus.end(); ++wit) {
        delete *wit;
    }
    _indus.clear();
    // ... and allocate new workers.
    for (int id = 0; id < this->_subthreads; ++id) {
        InduWorker *newIndu = new InduWorker(id, _top, this);
        _indus.push_back(newIndu);
        newIndu->InitSpheres(&_qmm, &_mm2);
        newIndu->SetSwitch(1);



    }
    
    this->InitChunks();
    
    return;
}


int XInductor::Induce(XJob *job) {
    
    for (int id = 0; id < this->_subthreads; ++id) {
        _indus[id]->SetSwitch(1);
    }

    vector< PolarSeg* >   ::iterator sit1;
    vector< PolarSeg* >   ::iterator sit2;
    vector< APolarSite* > ::iterator pit1;
    vector< APolarSite* > ::iterator pit2;
    
    // CONVERGENCE PARAMETERS
    double wSOR = this->_wSOR_N;
    double eTOL = this->_epsTol;
    int    maxI = this->_maxIter;
    
    // Adapt wSOR if any segment in QM region is charged
    for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
        double Q = (*sit1)->CalcTotQ();
        if (Q*Q > 0.5) { 
            wSOR = this->_wSOR_C;
            break;
        }
    }
    
    LOG(logINFO,*_log) << "Inductor: Using WSOR = " << wSOR 
        << ", ASHARP = " << _aDamp << flush;

    // Intra-pair induction ...
    bool   induce_intra_pair = this->_induce_intra_pair;
    // ... change this for jobs of type "site":
    if (_qmm.size() == 1) { induce_intra_pair = true; }
    
    
    // Compute CoM positions, generate coarsegrained sites
    vector<PolarFrag*>::iterator fit;
    for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
        (*sit1)->CalcPos();
        (*sit1)->GeneratePermInduCgSite(false);
        for (fit = (*sit1)->PolarFrags().begin(); 
            fit < (*sit1)->PolarFrags().end(); ++fit) {
            (*fit)->CalcPosCenterOfGeom();
            (*fit)->GeneratePermInduCgSite(false);
        }
    }
    
    double r_switch_cg_frag = 100.; // 2.5;
    double r_switch_cg_seg  = 100.; // 4.5;
    
    // ++++++++++++++++++++++++++++++++++++++++++++++ //
    // Inter-site fields (arising from perm. m'poles) //
    // ++++++++++++++++++++++++++++++++++++++++++++++ //
    
    boost::timer::cpu_timer cpu_t_perm;
    cpu_t_perm.start();
    boost::timer::cpu_times t_perm_0 = cpu_t_perm.elapsed();
    
    for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
    for (sit2 = sit1 + 1; sit2 < _qmm.end(); ++sit2) {

        // Intra-pair permanent induction field?
         if ( !induce_intra_pair ) {
             if ( job->isInCenter((*sit1)->getId())
               && job->isInCenter((*sit2)->getId()) ) {
                 continue;
             }
         }
//         for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//         for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
//             _actor.BiasStat(*(*pit1), *(*pit2));
//             _actor.FieldPerm(*(*pit1), *(*pit2));
//         }}
         
         double dr = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
         if (dr > r_switch_cg_seg) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                _actor.BiasStat(*(*pit1), *((*sit2)->getPermCgSite()));
                _actor.FieldPermAsPerm_At_By(*(*pit1), *((*sit2)->getPermCgSite()));
            }
            for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                _actor.BiasStat(*(*pit2), *((*sit1)->getPermCgSite()));
                _actor.FieldPermAsPerm_At_By(*(*pit2), *((*sit1)->getPermCgSite()));
            }
        }
        // Interaction atom <> fragment (cg)
        else if (dr > r_switch_cg_frag) {
            for (fit = (*sit2)->PolarFrags().begin(); fit < (*sit2)->PolarFrags().end(); ++fit) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    _actor.BiasStat(*(*pit1), *((*fit)->getPermCgSite()));
                    _actor.FieldPermAsPerm_At_By(*(*pit1), *((*fit)->getPermCgSite()));
                }
            }
            for (fit = (*sit1)->PolarFrags().begin(); fit < (*sit1)->PolarFrags().end(); ++fit) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit2), *((*fit)->getPermCgSite()));
                    _actor.FieldPermAsPerm_At_By(*(*pit2), *((*fit)->getPermCgSite()));
                }
            }
        }
        // Interaction atom <> atom
        else {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                _actor.BiasStat(*(*pit1),*(*pit2));
                _actor.FieldPerm(*(*pit1), *(*pit2));
            }}
        }
    }}
    
    // Permanent fields generated by outer shell    
    // (Outer shell itself is treated as non-polarizable)
    for (sit2 = _mm2.begin();
         sit2 < _mm2.end();
         ++sit2) {
    for (sit1 = _qmm.begin(); 
         sit1 < _qmm.end();
         ++sit1) {
         for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
         for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
             _actor.BiasStat(*(*pit1), *(*pit2));
             _actor.FieldPerm(*(*pit1), *(*pit2));
         }}
    }}
    
    boost::timer::cpu_times t_perm_1 = cpu_t_perm.elapsed();
    double t_perm = (t_perm_1.wall-t_perm_0.wall)/1e9;
    LOG(logINFO,*_log) << (format("  PERM      |  T=%1$1.2fs ") 
        % (t_perm));
    LOG(logINFO,*_log) << flush;

    // +++++++++++++++++++ //
    // 1st-order induction //
    // +++++++++++++++++++ //

    // Direct induction. Could also be restored from file (in the case of
    // iterative QM/MM being performed)
    if (true || !job->StartFromCPT()) {
        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                (*pit1)->InduceDirect();
            }
        }
    }
    else {
        LOG(logINFO,*_log) << "Job started from archive, reusing existing "
            << "state as starting configuration" << flush;
    }
    
    
    
    // Transition ranges atm<>atm, atm<>frag, atm<>seg
    vector<double> rcs_frag;
    vector<double> rcs_seg;
    vector<double> ets;
    //rcs_frag.push_back(1.5);
    //rcs_frag.push_back(2.0);
    rcs_frag.push_back(r_switch_cg_frag);
    //rcs_seg.push_back(3.0);
    //rcs_seg.push_back(3.75);
    rcs_seg.push_back(r_switch_cg_seg);
    //ets.push_back(100*eTOL);
    //ets.push_back(50*eTOL);
    ets.push_back(eTOL);
    
    int iter_cg = 0;
    
    for (int ridx = 0; ridx < 1; ++ridx) {
        double rc_frag = rcs_frag[ridx];
        double rc_seg = rcs_seg[ridx];
        double rc_eTOL = ets[ridx];
        iter_cg = 0;
        for ( ; iter_cg < maxI; ++iter_cg) {

            boost::timer::cpu_timer cpu_t;
            cpu_t.start();
            boost::timer::cpu_times t0 = cpu_t.elapsed();

            // Reset fields FUx, FUy, FUz
            for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    (*pit1)->ResetFieldU();
                }
            }

            boost::timer::cpu_times t1 = cpu_t.elapsed();
            // Coarsegrain polar segments & fragments
            for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
                (*sit1)->GeneratePermInduCgSite(false);
                for (fit = (*sit1)->PolarFrags().begin(); 
                    fit < (*sit1)->PolarFrags().end(); ++fit) {
                    (*fit)->GeneratePermInduCgSite(false);
                }
            }


            boost::timer::cpu_times t2 = cpu_t.elapsed();
            // Intra-site contribution to induction field
            for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = pit1 + 1;        pit2 < (*sit1)->end(); ++pit2) {
                    _actor.BiasIndu(*(*pit1),*(*pit2));
                    _actor.FieldIndu(*(*pit1),*(*pit2));
                }}
            }

            boost::timer::cpu_times t3 = cpu_t.elapsed();
            int count_atm_atm   = 0;
            int count_atm_frag  = 0;
            int count_atm_seg   = 0;
            // Inter-site contribution to induction field
            for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
            for (sit2 = sit1 + 1; sit2 < _qmm.end(); ++sit2) {
                double dr = votca::tools::abs((*sit1)->getPos()-(*sit2)->getPos());
                // Interaction atom <> segment (cg)
                if (dr > rc_seg) {
                    count_atm_seg += 1;
                    for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                        _actor.BiasIndu(*(*pit1), *((*sit2)->getInduCgSite()));
                        _actor.FieldPermAsIndu_At_By(*(*pit1), *((*sit2)->getInduCgSite()));
                    }
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        _actor.BiasIndu(*(*pit2), *((*sit1)->getInduCgSite()));
                        _actor.FieldPermAsIndu_At_By(*(*pit2), *((*sit1)->getInduCgSite()));
                    }
                }
                // Interaction atom <> fragment (cg)
                else if (dr > rc_frag) {
                    count_atm_frag += 1;
                    for (fit = (*sit2)->PolarFrags().begin(); fit < (*sit2)->PolarFrags().end(); ++fit) {
                        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                            _actor.BiasIndu(*(*pit1), *((*fit)->getInduCgSite()));
                            _actor.FieldPermAsIndu_At_By(*(*pit1), *((*fit)->getInduCgSite()));
                        }
                    }
                    for (fit = (*sit1)->PolarFrags().begin(); fit < (*sit1)->PolarFrags().end(); ++fit) {
                        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                            _actor.BiasIndu(*(*pit2), *((*fit)->getInduCgSite()));
                            _actor.FieldPermAsIndu_At_By(*(*pit2), *((*fit)->getInduCgSite()));
                        }
                    }
                }
                 // Interaction atom <> atom
                else {
                    count_atm_atm += 1;
                    for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        _actor.BiasIndu(*(*pit1),*(*pit2));
                        _actor.FieldIndu(*(*pit1), *(*pit2));
                    }}
                }
            }}

            boost::timer::cpu_times t4 = cpu_t.elapsed();

            // Induce again
            for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
                 for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                     (*pit1)->Induce(wSOR);                                         
                 }
            }

            // Check for convergence
            bool    converged       = true;
            double  maxdU_U         = -1;
            double  avgdU_U         = 0.0;
            double  rmsdU           = 0.0;
            int     baseN           = 0;
            for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
                 for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                     double dU_U = (*pit1)->HistdU();
                     avgdU_U += dU_U;
                     double dU2 = (*pit1)->HistdU2();
                     rmsdU += dU2;
                     ++baseN;
                     if ( dU_U > maxdU_U ) { maxdU_U = dU_U; }
                     if ( dU_U > rc_eTOL) { converged = false; }
                 }
            }
            avgdU_U /= baseN;
            rmsdU /= baseN;
            rmsdU = sqrt(rmsdU);
            if (avgdU_U < rc_eTOL/10.) { converged = true; }


            boost::timer::cpu_times t5 = cpu_t.elapsed();


            double t_total = (t5.wall-t0.wall)/1e9;
            double t_cg    = (t2.wall-t1.wall)/1e9;
            double t_intra = (t3.wall-t2.wall)/1e9;
            double t_inter = (t4.wall-t3.wall)/1e9;

            LOG(logINFO,*_log) << (boost::format(
                "  ITER %1$3d  |  AVG(dU/U) %3$1.7e  EPS %11$1.1e  |  NN(0<>%7$1.1f<>%6$1.1f) %8$04d|%9$04d|%10$04d") 
                % iter_cg % maxdU_U % avgdU_U % rmsdU % baseN % rc_seg 
                % rc_frag % count_atm_atm % count_atm_frag % count_atm_seg 
                % eTOL).str();
            LOG(logINFO,*_log) << (format("  |  T=%1$1.2fs CG/I/O %2$2.2f%% %3$2.2f%% %4$2.2f%%") 
                % (t_total) % (100*t_cg/t_total) % (100*t_intra/t_total) 
                % (100*t_inter/t_total));
            LOG(logINFO,*_log) << flush;


            // Break if converged
            if      (converged) { 
                _isConverged = true;
                break; 
            }
            else if (iter_cg == maxI - 1) {
                _isConverged = false;
                this->setError((boost::format("Did not converge to precision "
                   "(%1$d steps, AVG(dU:U) = %2$1.3e)") % maxI % avgdU_U).str());
                break;
            }
        }
    }

    // Reset U1 history
    for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            (*pit1)->ResetU1Hist();
        }
    }

    return iter_cg;
//    assert(false);
//    
//    
//    // ++++++++++++++++++++++ //
//    // Higher-order induction //
//    // ++++++++++++++++++++++ //
//
//    
//    int iter = 0;
//    for ( ; iter < maxI; ++iter) {
//        // Reset fields FUx, FUy, FUz
//        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
//            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//                (*pit1)->ResetFieldU();
//            }
//        }
//
//        // Intra-site contribution to induction field
//        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
//            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//            for (pit2 = pit1 + 1;        pit2 < (*sit1)->end(); ++pit2) {
//                _actor.BiasIndu(*(*pit1),*(*pit2));
//                _actor.FieldIndu(*(*pit1),*(*pit2));
//            }}
//        }
//
//        // Inter-site contribution to induction field
//        //for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
//        //for (sit2 = sit1 + 1; sit2 < _qmm.end(); ++sit2) {
//        //    for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//        //    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
//        //        _actor.FieldIndu(*(*pit1), *(*pit2));
//        //    }}
//        //}}
//
//        for (int id = 0; id < this->_subthreads; ++id) { 
//            _indus[id]->Start();
//        }
//
//        for (int id = 0; id < this->_subthreads; ++id) {
//            _indus[id]->WaitDone();
//        }
//
//        this->ClearTodoTable();
//
//
//
//        // Induce again
//        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
//             for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//                 (*pit1)->Induce(wSOR);                                         
//             }
//        }
//
//        // Check for convergence
//        bool    converged       = true;
//        double  maxdU_U         = -1;
//        double  avgdU_U         = 0.0;
//        double  rmsdU           = 0.0;
//        int     baseN           = 0;
//        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
//             for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
//                 double dU_U = (*pit1)->HistdU();
//                 avgdU_U += dU_U;
//                 double dU2 = (*pit1)->HistdU2();
//                 rmsdU += dU2;
//                 ++baseN;
//                 if ( dU_U > maxdU_U ) { maxdU_U = dU_U; }
//                 if ( dU_U > eTOL) { converged = false; }
//             }
//        }
//        avgdU_U /= baseN;
//        rmsdU /= baseN;
//        rmsdU = sqrt(rmsdU);
//        if (avgdU_U < eTOL/10.) { converged = true; }
//        
//        LOG(logINFO,*_log) << (boost::format(
//            "  Iter %1$3d | max(dU/U) %2$1.7e  avg(dU/U) %3$1.7e  rms(dU) %4$1.7e  N %5$d") 
//            % iter % maxdU_U % avgdU_U % rmsdU % baseN).str() << flush;
//
//        // Break if converged
//        if      (converged) { 
//            _isConverged = true;
//            break; 
//        }
//        else if (iter == maxI - 1) {
//            _isConverged = false;
//            //this->setError("Did not converge to precision (" 
//            //    + boost::lexical_cast<string>(maxI) + " steps, AVG dU:U " 
//            //    + boost::lexical_cast<string>(avgdU) + ")");
//            this->setError((boost::format("Did not converge to precision "
//               "(%1$d steps, AVG(dU:U) = %2$1.3e)") % maxI % avgdU_U).str());
//            break;
//        }
//    }    
//    
//    return iter;

}


double XInductor::Energy(XJob *job) {

    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;    

    _actor.ResetEnergy();
    
    // Energy splittings =======================================================
    // PAIR/SITE        <->        SPH1         <->          SPH2 = OUT       //
    double E_Tot = 0.0;
    // ... 0th kind    
    double E_Pair_Pair = 0.0;    
    double E_Pair_Sph1 = 0.0;
    double E_Sph1_Sph1 = 0.0;    
    double E_Pair_Sph2 = 0.0;
    double E_Sph1_Sph2 = 0.0;
    // ... 1st kind
    double eu_inter = 0.0;
    double eu_intra = 0.0;
    double e_perm   = 0.0;
    // ... 2nd kind
    double epp      = 0.0;
    double epu      = 0.0;
    double euu      = 0.0;
    // ... 3rd kind
    double e_f_c_c          = 0.0;
    double e_f_c_non_c      = 0.0;
    double e_f_c_out        = 0.0;
    double e_f_non_c_non_c  = 0.0;   
    double e_f_non_c_out    = 0.0;
    double e_m_c            = 0.0;
    double e_m_c_out        = 0.0;
    double e_m_non_c        = 0.0;
    double e_m_non_c_out    = 0.0;
    double e_m_out          = 0.0;
    // =========================================================================

    vector< PolarSeg* >     ::iterator      sit1;
    vector< PolarSeg* >     ::iterator      sit2;
    vector< APolarSite* >   ::iterator      pit1;
    vector< APolarSite* >   ::iterator      pit2;

    
    // =============================================================== //
    // System Energy | QM | MM1 | MM2 |                                //
    // =============================================================== //
    

    for (int id = 0; id < this->_subthreads; ++id) {
        _indus[id]->SetSwitch(0);
    }

    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Inter-site energy comprising central + first polarization shell //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //


    for (int id = 0; id < this->_subthreads; ++id) {
        _indus[id]->Start();
    }

    for (int id = 0; id < this->_subthreads; ++id) {
        _indus[id]->WaitDone();
    }

    for (int id = 0; id < this->_subthreads; ++id) {
        E_Pair_Pair += _indus[id]->GetEPairPair();
        E_Pair_Sph1 += _indus[id]->GetEPairSph1();
        E_Sph1_Sph1 += _indus[id]->GetESph1Sph1();

        eu_inter += _indus[id]->GetActor().getEU_INTER();
        eu_intra += _indus[id]->GetActor().getEU_INTRA();
        e_perm   += _indus[id]->GetActor().getEP();

        epp += _indus[id]->GetActor().getEPP();
        epu += _indus[id]->GetActor().getEPU();
        euu += _indus[id]->GetActor().getEUU();

        e_f_c_c             += _indus[id]->GetE_f_C_C();
        e_f_c_non_c         += _indus[id]->GetE_f_C_non_C();
        e_f_non_c_non_c     += _indus[id]->GetE_f_non_C_non_C();            
        e_m_c               += _indus[id]->GetE_m_C();
        e_m_non_c           += _indus[id]->GetE_m_non_C();
    }

    this->ClearTodoTable(); 


    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Inter-site energy resulting from interaction with static shell  //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    // Interaction between central and static shell
    for (sit1 = _qm0.begin(); sit1 < _qm0.end(); ++sit1) {
    for (sit2 = _mm2.begin(); sit2 < _mm2.end(); ++sit2) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_c_out += _actor.E_f(*(*pit1), *(*pit2));
            e_m_c_out += _actor.E_m(*(*pit1), *(*pit2));            
        }}
    }}
    
    // Interaction between polarizable and static shell
    for (sit1 = _mm1.begin(); sit1 < _mm1.end(); ++sit1) {
    for (sit2 = _mm2.begin(); sit2 < _mm2.end(); ++sit2) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_non_c_out += _actor.E_f(*(*pit1), *(*pit2));
            e_m_non_c_out += _actor.E_m(*(*pit1), *(*pit2));            
        }}
    }}

    // Increment energies
    // ... 0th kind        
    E_Pair_Sph2 += e_f_c_out + e_m_c_out;
    E_Sph1_Sph2 += e_f_non_c_out + e_m_non_c_out;
    // ... 1st kind
    e_perm      += _actor.getEP();
    eu_inter    += _actor.getEU_INTER();
    // ... 2nd kind
    epp += _actor.getEPP();
    epu += _actor.getEPU();
    euu += _actor.getEUU();
    // ... 3rd kind
    // ... ... -> done in loop above, but need to summarize e_m_*
    e_m_c      += e_m_c_out;
    e_m_non_c  += e_m_non_c_out;


    E_Tot = E_Pair_Pair + E_Pair_Sph1 + E_Sph1_Sph1 + E_Pair_Sph2;
    
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Intramolecular field interaction                                //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    double e_f_intra_0 = 0.0;
    double e_m_intra_0 = 0.0;
    for (sit1 = _qm0.begin(); sit1 < _qm0.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = pit1+1; pit2 < (*sit1)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_intra_0 += _actor.E_f_intra(*(*pit1), *(*pit2));
            e_m_intra_0 += _actor.E_m_intra(*(*pit1), *(*pit2));
            _actor.RevBias();
            e_m_intra_0 += _actor.E_m_intra(*(*pit2), *(*pit1));
        }}
    }
    double e_f_intra_1 = 0.0;
    double e_m_intra_1 = 0.0;
    for (sit1 = _mm1.begin(); sit1 < _mm1.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = pit1+1; pit2 < (*sit1)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_intra_1 += _actor.E_f_intra(*(*pit1), *(*pit2));
            e_m_intra_1 += _actor.E_m_intra(*(*pit1), *(*pit2));
            _actor.RevBias();
            e_m_intra_1 += _actor.E_m_intra(*(*pit2), *(*pit1));
        }}
    }
    double e_f_intra_2 = 0.0;
    double e_m_intra_2 = 0.0;
    for (sit1 = _mm2.begin(); sit1 < _mm2.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = pit1+1; pit2 < (*sit1)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_intra_2 += _actor.E_f_intra(*(*pit1), *(*pit2));
            e_m_intra_2 += _actor.E_m_intra(*(*pit1), *(*pit2));
            _actor.RevBias();
            e_m_intra_2 += _actor.E_m_intra(*(*pit2), *(*pit1));
        }}
    }
    
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Total ind. work (to be used with pre-generated perm. fields)    //
    // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    double e_m_0 = 0.0;
    double e_m_1 = 0.0;
    double e_m_2 = 0.0;
    for (sit1 = _qm0.begin(); sit1 < _qm0.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            e_m_0 += (*pit1)->InductionWork();
        }
    }
    for (sit1 = _mm1.begin(); sit1 < _mm1.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            e_m_1 += (*pit1)->InductionWork();
        }
    }
    for (sit1 = _mm2.begin(); sit1 < _mm2.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            e_m_2 += (*pit1)->InductionWork();
        }
    }
    
//    cout << endl << "E_f_intra_0 " << e_f_intra_0*int2eV << flush;
//    cout << endl << "E_m_intra_0 " << e_m_intra_0*int2eV << flush;
//    cout << endl << "E_f_intra_1 " << e_f_intra_1*int2eV << flush;
//    cout << endl << "E_m_intra_1 " << e_m_intra_1*int2eV << flush;
//    cout << endl << "E_f_intra_2 " << e_f_intra_2*int2eV << flush;
//    cout << endl << "E_m_intra_2 " << e_m_intra_2*int2eV << flush;
//    cout << endl << "E_m_0 " << e_m_0*int2eV << flush;
//    cout << endl << "E_m_1 " << e_m_1*int2eV << flush;
//    cout << endl << "E_m_2 " << e_m_2*int2eV << flush;
    
    // This is important if permanent/induction fields have been applied
    // that do not originate in QM0, MM1, MM2
    
    e_m_c     = e_m_0 + e_f_intra_0;
    e_m_non_c = e_m_1 + e_f_intra_1;
    e_m_out   = e_m_2 + e_f_intra_2;
    
    //e_f_c_c += e_f_intra_0;
    //e_f_non_c_non_c += e_f_intra_1;
    
    
    
    
    
    
    
    

    // =============================================================== //
    // Energy Output                                                   //
    // =============================================================== //

    // ... 0th kind
    E_Tot = E_Pair_Pair
          + E_Pair_Sph1
          + E_Sph1_Sph1
          + E_Pair_Sph2
          + E_Sph1_Sph2;

//    LOG(logINFO,*_log) 
//             << "... E(X) = " << E_Tot * int2eV << " eV "
//             << flush << "...      = (Site, Site) " << E_Pair_Pair * int2eV
//             << flush << "...      + (Site, Sph1) " << E_Pair_Sph1 * int2eV
//             << flush << "...      + (Sph1, Sph1) " << E_Sph1_Sph1 * int2eV
//             << flush << "...      + (Site, Sph2) " << E_Pair_Sph2 * int2eV
//             << flush << "...      + (Sph1, Sph2) " << E_Sph1_Sph2 * int2eV
//             << flush;

    // ... 1st kind
    double E_PPUU = epp 
                  + epu 
                  + euu;

//    LOG(logINFO,*_log)
//             << "... E(X) = " << E_PPUU * int2eV << " eV " 
//             << flush << "...      = (PP) "    << epp  * int2eV
//             << flush << "...      + (PU) "    << epu  * int2eV
//             << flush << "...      + (UU) "    << euu  * int2eV
//             << flush;

    // ... 2nd kind
    double E_f_m = e_f_c_c 
                 + e_f_c_non_c
                 + e_f_c_out 
                 + e_f_non_c_non_c 
                 + e_f_non_c_out
                 + e_m_c 
                 + e_m_non_c
                 + e_m_out;

    LOG(logINFO,*_log)
        << (format("PP-PU-UU Splitting")).str()
        << flush << (format("  + U [Q  -  Q]       = %1$+1.7f eV") % (epp      * int2eV)).str()
        << flush << (format("  + U [Q  - dQ]       = %1$+1.7f eV") % (epu      * int2eV)).str()       
        << flush << (format("  + U [dQ - dQ]       = %1$+1.7f eV") % (euu      * int2eV)).str()
        << flush << (format("  = ------------------------------")).str()
        << flush << (format("  + SUM(E)            = %1$+1.7f eV") % (E_PPUU   * int2eV)).str()
        << flush;
    LOG(logINFO,*_log)
        << (format("QM0-MM1-MM2 Splitting")).str()
        << flush << (format("  + Field-term [0-0]  = %1$+1.7f eV (%2$+1.7f)") % (e_f_c_c          * int2eV) % (e_f_intra_0*int2eV)).str()
        << flush << (format("  + Field-term [0-1]  = %1$+1.7f eV") % (e_f_c_non_c      * int2eV)).str()       
        << flush << (format("  + Field-term [0-2]  = %1$+1.7f eV") % (e_f_c_out        * int2eV)).str()
        << flush << (format("  + Field-term [1-1]  = %1$+1.7f eV (%2$+1.7f)") % (e_f_non_c_non_c  * int2eV) % (e_f_intra_1*int2eV)).str()
        << flush << (format("  + Field-term [1-2]  = %1$+1.7f eV") % (e_f_non_c_out    * int2eV)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + Work-term  [-0-]  = %1$+1.7f eV") % (e_m_c            * int2eV)).str()
        << flush << (format("  + Work-term  [-1-]  = %1$+1.7f eV") % (e_m_non_c        * int2eV)).str()
        << flush << (format("  + Work-term  [-2-]  = %1$+1.7f eV") % (e_m_out          * int2eV)).str()
        << flush << (format("  = ------------------------------")).str()
        << flush << (format("    SUM(E)            = %1$+1.7f eV") % (E_f_m               *int2eV)).str()
        << flush;
    
    
//    LOG(logINFO,*_log)
//             << "... E(X) = " << E_f_m * int2eV << " eV " 
//             << flush << "...      = (f,0-0) " << e_f_c_c          * int2eV
//             << flush << "...      + (f,0-1) " << e_f_c_non_c      * int2eV
//             << flush << "...      + (f,0-2) " << e_f_c_out        * int2eV
//             << flush << "...      + (f,1-1) " << e_f_non_c_non_c  * int2eV
//             << flush << "...      + (f,1-2) " << e_f_non_c_out    * int2eV
//             << flush << "...      + (m,-0-) " << e_m_c            * int2eV
//             << flush << "...      + (m,-1-) " << e_m_non_c        * int2eV
//             << flush << "...      + (m,-2-) " << e_m_out          * int2eV
//             << flush;

    // Forward results to job
    job->setEnergy(E_Tot            *int2eV,           
                   E_Pair_Pair      *int2eV,
                   E_Pair_Sph1      *int2eV,
                   E_Pair_Sph2      *int2eV, 
                   E_Sph1_Sph1      *int2eV,
                   E_Sph1_Sph2      *int2eV,                       
                   e_perm           *int2eV,
                   eu_inter         *int2eV);

    job->setEnergy_PPUU(epp         *int2eV,
                        epu         *int2eV,
                        euu         *int2eV);

    job->setEnergy_f_m(e_f_c_c         *int2eV,
                       e_f_c_non_c     *int2eV,
                       e_f_c_out       *int2eV,
                       e_f_non_c_non_c *int2eV, 
                       e_f_non_c_out   *int2eV,
                       e_m_c           *int2eV, 
                       e_m_non_c       *int2eV,
                       e_m_out         *int2eV);

    return E_Tot;
}


double XInductor::EnergyStatic(XJob *job) {
    
    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;

    _actor.ResetEnergy();
    
    // Energy splittings =======================================================
    // PAIR/SITE        <->        SPH1         <->          SPH2 = OUT       //
    double E_Tot = 0.0;
    // ... 0th kind    
    double E_Pair_Pair = 0.0;    
    double E_Pair_Sph1 = 0.0;
    double E_Sph1_Sph1 = 0.0;    
    double E_Pair_Sph2 = 0.0;
    double E_Sph1_Sph2 = 0.0;
    // ... 1st kind
    double eu_inter = 0.0;
    //double eu_intra = 0.0;
    double e_perm   = 0.0;
    // ... 2nd kind
    double epp      = 0.0;
    double epu      = 0.0;
    double euu      = 0.0;
    // ... 3rd kind
    double e_f_c_c          = 0.0;
    double e_f_c_non_c      = 0.0;
    double e_f_c_out        = 0.0;
    double e_f_non_c_non_c  = 0.0;   
    double e_f_non_c_out    = 0.0;
    double e_m_c            = 0.0;
    //double e_m_c_out        = 0.0;
    double e_m_non_c        = 0.0;
    //double e_m_non_c_out    = 0.0;
    double e_m_out          = 0.0;
    // =========================================================================

    
    vector< PolarSeg* >    ::iterator      sit1;
    vector< PolarSeg* >    ::iterator      sit2;
    vector< APolarSite* >  ::iterator      pit1;
    vector< APolarSite* >  ::iterator      pit2;

        
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Interaction pair <-> inner cut-off, without intra-pair interaction //
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    for (sit1 = _qm0.begin(); sit1 < _qm0.end(); ++sit1) {
    for (sit2 = _mm1.begin(); sit2 < _mm1.end(); ++sit2) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_c_non_c += _actor.E_f(*(*pit1), *(*pit2));
        }}
    }}


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Interaction pair <-> outer cut-off                                 //
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

    for (sit1 = _qm0.begin(); sit1 < _qm0.end(); ++sit1) {
    for (sit2 = _mm2.begin(); sit2 < _mm2.end(); ++sit2) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_c_out += _actor.E_f(*(*pit1), *(*pit2));            
        }}
    }}


    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    // Intra-pair interaction                                             //
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
    
    for (sit1 = _qm0.begin(); sit1 < _qm0.end(); ++sit1) {
    for (sit2 = sit1 + 1; sit2 < _qm0.end(); ++sit2) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
            _actor.BiasIndu(*(*pit1), *(*pit2));
            e_f_c_c += _actor.E_f(*(*pit1), *(*pit2));            
        }}
    }}


    // Increment energies
    // ... 0th kind
    E_Pair_Pair     += e_f_c_c;
    E_Pair_Sph1     += e_f_c_non_c;
    E_Pair_Sph2     += e_f_c_out;
    // ... 1st kind
    e_perm          += _actor.getEP();
    eu_inter        += _actor.getEU_INTER();
    // ... 2nd kind
    epp             += _actor.getEPP();
    epu             += _actor.getEPU();
    euu             += _actor.getEUU();
    // ... 3rd kind
    // ... ... -> done in loops above        


    // =============================================================== //
    // Energy Output                                                   //
    // =============================================================== //

    // ... 0th kind
    E_Tot = E_Pair_Pair
          + E_Pair_Sph1
          + E_Sph1_Sph1
          + E_Pair_Sph2
          + E_Sph1_Sph2;
    
//    LOG(logINFO,*_log) 
//             << "... E(X) = " << E_Tot * int2eV << " eV "
//             << flush << "...      = (Site, Site) " << E_Pair_Pair * int2eV
//             << flush << "...      + (Site, Sph1) " << E_Pair_Sph1 * int2eV
//             << flush << "...      + (Sph1, Sph1) " << E_Sph1_Sph1 * int2eV
//             << flush << "...      + (Site, Sph2) " << E_Pair_Sph2 * int2eV
//             << flush << "...      + (Sph1, Sph2) " << E_Sph1_Sph2 * int2eV
//             << flush;

    // ... 1st kind
    double E_PPUU = epp 
                  + epu 
                  + euu;

//    LOG(logINFO,*_log)
//             << "... E(X) = " << E_PPUU * int2eV << " eV " 
//             << flush << "...      = (PP) "    << epp  * int2eV
//             << flush << "...      + (PU) "    << epu  * int2eV
//             << flush << "...      + (UU) "    << euu  * int2eV
//             << flush;

    // ... 2nd kind
    double E_f_m = e_f_c_c 
                 + e_f_c_non_c
                 + e_f_c_out 
                 + e_f_non_c_non_c 
                 + e_f_non_c_out
                 + e_m_c 
                 + e_m_non_c
                 + e_m_out;

//    LOG(logINFO,*_log)
//             << "... E(X) = " << E_f_m * int2eV << " eV " 
//             << flush << "...      = (f,0-0) " << e_f_c_c          * int2eV
//             << flush << "...      + (f,0-1) " << e_f_c_non_c      * int2eV
//             << flush << "...      + (f,0-2) " << e_f_c_out        * int2eV
//             << flush << "...      + (f,1-1) " << e_f_non_c_non_c  * int2eV
//             << flush << "...      + (f,1-2) " << e_f_non_c_out    * int2eV
//             << flush << "...      + (m,-0-) " << e_m_c            * int2eV
//             << flush << "...      + (m,-1-) " << e_m_non_c        * int2eV
//             << flush << "...      + (m,-2-) " << e_m_out          * int2eV
//             << flush;

    
    LOG(logINFO,*_log)
        << (format("PP-PU-UU Splitting")).str()
        << flush << (format("  + U [Q  -  Q]       = %1$+1.7f eV") % (epp      * int2eV)).str()
        << flush << (format("  + U [Q  - dQ]       = %1$+1.7f eV") % (epu      * int2eV)).str()       
        << flush << (format("  + U [dQ - dQ]       = %1$+1.7f eV") % (euu      * int2eV)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("    SUM(E)            = %1$+1.7f eV") % (E_PPUU   * int2eV)).str()
        << flush;
    LOG(logINFO,*_log)
        << (format("QM0-MM1-MM2 Splitting")).str()
        << flush << (format("  + Field-term [0-0]  = %1$+1.7f eV") % (e_f_c_c          * int2eV)).str()
        << flush << (format("  + Field-term [0-1]  = %1$+1.7f eV") % (e_f_c_non_c      * int2eV)).str()       
        << flush << (format("  + Field-term [0-2]  = %1$+1.7f eV") % (e_f_c_out        * int2eV)).str()
        << flush << (format("  + Field-term [1-1]  = %1$+1.7f eV") % (e_f_non_c_non_c  * int2eV)).str()
        << flush << (format("  + Field-term [1-2]  = %1$+1.7f eV") % (e_f_non_c_out    * int2eV)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + Work-term  [-0-]  = %1$+1.7f eV") % (e_m_c            * int2eV)).str()
        << flush << (format("  + Work-term  [-1-]  = %1$+1.7f eV") % (e_m_non_c        * int2eV)).str()
        << flush << (format("  + Work-term  [-2-]  = %1$+1.7f eV") % (e_m_out          * int2eV)).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("    SUM(E)            = %1$+1.7f eV") % (E_f_m               *int2eV)).str()
        << flush;
    
    
    // Forward results to job
    job->setEnergy(E_Tot            *int2eV,           
                   E_Pair_Pair      *int2eV,
                   E_Pair_Sph1      *int2eV,
                   E_Pair_Sph2      *int2eV, 
                   E_Sph1_Sph1      *int2eV,
                   E_Sph1_Sph2      *int2eV,                       
                   e_perm           *int2eV,
                   eu_inter         *int2eV);

    job->setEnergy_PPUU(epp         *int2eV,
                        epu         *int2eV,
                        euu         *int2eV);

    job->setEnergy_f_m(e_f_c_c         *int2eV,
                       e_f_c_non_c     *int2eV,
                       e_f_c_out       *int2eV,
                       e_f_non_c_non_c *int2eV, 
                       e_f_non_c_out   *int2eV,
                       e_m_c           *int2eV, 
                       e_m_non_c       *int2eV,
                       e_m_out         *int2eV);    
    
    return E_Tot;
}


XInductor::XInductor(Topology *top, Property *opt, 
                     string sfx, int nst, bool mav)
                  : _subthreads(nst), _maverick(mav) {
    
    string key = sfx + ".tholemodel";

        if ( opt->exists(key+".induce") ) {
            int induce = opt->get(key+".induce").as< int >();
            _induce = (induce == 0) ? false : true;
        }
        else { _induce = true; }

        if ( opt->exists(key+".induce_intra_pair") ) {
            int induce = opt->get(key+".induce_intra_pair").as< int >();
            _induce_intra_pair = (induce == 0) ? false : true;
        }
        else { _induce_intra_pair = true; }

        if ( opt->exists(key+".exp_damp") ) {
            _aDamp = opt->get(key+".exp_damp").as< double >();
        }
        else {
            cout << endl << "... ... WARNING: No sharpness parameter supplied";
            cout << endl << "... ... ... Using default a = 0.39";
            _aDamp = 0.39;
        }

    key = sfx + ".convergence";

        if ( opt->exists(key+".wSOR_N") ) {
            _wSOR_N = opt->get(key+".wSOR_N").as< float >();
        }
        else { _wSOR_N = 0.75; }
        if ( opt->exists(key+".wSOR_C") ) {
            _wSOR_C = opt->get(key+".wSOR_C").as< float >();
        }
        else { _wSOR_C = 0.75; }

        if ( opt->exists(key+".max_iter") ) {
            _maxIter = opt->get(key+".max_iter").as< int >();
        }
        else { _maxIter = 512; }

        if ( opt->exists(key+".tolerance") ) {
            _epsTol = opt->get(key+".tolerance").as< double >();
        }
        else { _epsTol = 0.001; }
    
    _actor = XInteractor(NULL, _aDamp);
    
}


}}