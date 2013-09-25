#include <votca/ctp/xinductor.h>
#include <boost/format.hpp>
#include <vector>

using boost::format;

namespace votca { namespace ctp {


void XInductor::Evaluate(XJob *job) {    
    
    _qm0.clear();
    _mm1.clear();
    _mm2.clear();
    _qmm.clear();
    
    _job = job;
    
    // QMM = [ QM0 --- MM1 ] --- MM2 = [ MM2 ]    
    _qm0 = job->getPolarTop()->QM0();
    _mm1 = job->getPolarTop()->MM1();
    _mm2 = job->getPolarTop()->MM2();    

    _qmm.reserve(_qm0.size()+_mm1.size());
    _qmm.insert(_qmm.end(),_qm0.begin(),_qm0.end());
    _qmm.insert(_qmm.end(),_mm1.begin(),_mm1.end());

    // ++++++++++++++++++++++++++ //
    // (De-)polarize, charge to N //
    // ++++++++++++++++++++++++++ //

    if (job->StartFromCPT()) {
        // Permanent fields already computed, do not zero these out ...
        LOG(logDEBUG,*_log) << "Carry out partial depolarization." << flush;
        vector< PolarSeg* >   ::iterator sit;
        vector< APolarSite* > ::iterator pit;

        // Partially depolarize inner sphere
        for (sit = _qmm.begin(); sit < _qmm.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->ResetFieldU();
            (*pit)->ResetU1Hist();
            (*pit)->ResetU1();
            (*pit)->Charge(0); // <- Not necessarily neutral state
        }}

        // Partially depolarize outer shell
        for (sit = _mm2.begin(); sit < _mm2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->ResetFieldU();
            (*pit)->ResetU1Hist();
            (*pit)->ResetU1();
            (*pit)->Charge(0); // <- Not necessarily neutral state
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
    
    // +++++++++++++++++ //
    // Induction workers //
    // +++++++++++++++++ //

    for (int id = 0; id < this->_subthreads; ++id) {
        InduWorker *newIndu = new InduWorker(id, _top, this);
        _indus.push_back(newIndu);
        newIndu->InitSpheres(&_qmm, &_mm2);
        newIndu->SetSwitch(1);
    }
    
    this->InitChunks();
    
    // ++++++++++++++++++++++++++ //
    // Compute state energy       //
    // ++++++++++++++++++++++++++ //

    double  E_state  = 0.0;
    int     iter     = 0;

    if (this->_induce) iter      = this->Induce(job);
    if (this->_induce) E_state   = this->Energy(job);
    else               E_state   = this->EnergyStatic(job);

    job->setInduIter(iter);
    
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
    
    
    // ++++++++++++++++++++++++++++++++++++++++++++++ //
    // Inter-site fields (arising from perm. m'poles) //
    // ++++++++++++++++++++++++++++++++++++++++++++++ //
    
    for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
    for (sit2 = sit1 + 1; sit2 < _qmm.end(); ++sit2) {

        // Intra-pair permanent induction field?
         if ( !induce_intra_pair ) {
             if ( job->isInCenter((*sit1)->getId())
               && job->isInCenter((*sit2)->getId()) ) {
                 continue;
             }
         }
         for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
         for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
             _actor.BiasStat(*(*pit1), *(*pit2));
             _actor.FieldPerm(*(*pit1), *(*pit2));
         }}
    }}
    
    // Permanent fields generated by outer shell    
    // (Outer shell itself is treated as non-polarizable)
    for (sit1 = _qmm.begin(); 
         sit1 < _qmm.end();
         ++sit1) {
    for (sit2 = _mm2.begin();
         sit2 < _mm2.end();
         ++sit2) {
         for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
         for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
             _actor.BiasStat(*(*pit1), *(*pit2));
             _actor.FieldPerm(*(*pit1), *(*pit2));
         }}
    }}
    

    // +++++++++++++++++++ //
    // 1st-order induction //
    // +++++++++++++++++++ //

    // Direct induction. Could also be restored from file (in the case of
    // iterative QM/MM being performed)
    for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
         for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
             (*pit1)->InduceDirect();
         }
    }


    // ++++++++++++++++++++++ //
    // Higher-order induction //
    // ++++++++++++++++++++++ //

    
    int iter = 0;
    for ( ; iter < maxI; ++iter) {

        // Reset fields FUx, FUy, FUz
        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                (*pit1)->ResetFieldU();
            }
        }

        // Intra-site contribution to induction field
        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
            for (pit2 = pit1 + 1;        pit2 < (*sit1)->end(); ++pit2) {
                _actor.BiasIndu(*(*pit1),*(*pit2));
                _actor.FieldIndu(*(*pit1),*(*pit2));
            }}
        }

        // Inter-site contribution to induction field
        //for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
        //for (sit2 = sit1 + 1; sit2 < _qmm.end(); ++sit2) {
        //    for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        //    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
        //        _actor.FieldIndu(*(*pit1), *(*pit2));
        //    }}
        //}}

        for (int id = 0; id < this->_subthreads; ++id) { 
            _indus[id]->Start();
        }

        for (int id = 0; id < this->_subthreads; ++id) {
            _indus[id]->WaitDone();
        }

        this->ClearTodoTable();



        // Induce again
        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
             for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                 (*pit1)->Induce(wSOR);                                         
             }
        }

        // Check for convergence
        bool    converged       = true;
        double  maxdU           = -1;
        double  avgdU           = 0.0;
        int     baseN           = 0;
        for (sit1 = _qmm.begin(); sit1 < _qmm.end(); ++sit1) {
             for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                 double dU = (*pit1)->HistdU();
                 avgdU += dU;
                 ++baseN;
                 if ( dU > maxdU ) { maxdU = dU; }
                 if ( dU > eTOL ) { converged = false; }
             }
        }
        avgdU /= baseN;
        if (avgdU < eTOL/10.) { converged = true; }

//        cout << " | MAX dU " << maxdU
//             << " | AVG dU " << avgdU
//             << " | SOR " << wSOR << flush;

        // Break if converged
        if      (converged) { 
            _isConverged = true;
            break; 
        }
        else if (iter == maxI - 1) {
            _isConverged = false;
            //this->setError("Did not converge to precision (" 
            //    + boost::lexical_cast<string>(maxI) + " steps, AVG dU:U " 
            //    + boost::lexical_cast<string>(avgdU) + ")");
            this->setError((boost::format("Did not converge to precision "
               "(%1$d steps, AVG(dU:U) = %2$1.3e)") % maxI % avgdU).str());
            break;
        }
    }

    return iter;

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