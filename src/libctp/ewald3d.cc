#include <votca/ctp/ewald3d.h>
#include <boost/format.hpp>
#include <algorithm>


using boost::format;


namespace votca { namespace ctp {


Ewald3D::~Ewald3D() {
    vector< PolarSeg* >::iterator sit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();
}
    
    
Ewald3D::Ewald3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
    : _top(top), _ptop(ptop), _log(log) {
    
    // EVALUATE OPTIONS
    _R_co = opt->get("options.ewald.coulombmethod.cutoff").as<double>();
    _shape = opt->get("options.ewald.coulombmethod.shape").as<string>();
    _crit_dE = opt->get("options.ewald.convergence.energy").as<double>();    
    
    // EWALD INTERACTION PARAMETERS (GUESS ONLY)
    _K_co = 100/_R_co;
    _alpha = 3.5/_R_co;
    
    // SET-UP REAL & RECIPROCAL SPACE
    _a = _top->getBox().getCol(0);
    _b = _top->getBox().getCol(1);
    _c = _top->getBox().getCol(2);
    _LxLyLz = _a*(_b^_c);
    _LxLy = abs(_a ^ _b);
    
    _A = 2*M_PI/_LxLyLz * _b^_c;
    _B = 2*M_PI/_LxLyLz * _c^_a;
    _C = 2*M_PI/_LxLyLz * _a^_b;

    _na_max = ceil(_R_co/maxnorm(_a)-0.5)+1;
    _nb_max = ceil(_R_co/maxnorm(_b)-0.5)+1;
    _nc_max = ceil(_R_co/maxnorm(_c)-0.5)+1;

    _NA_max = ceil(_K_co/maxnorm(_A));
    _NB_max = ceil(_K_co/maxnorm(_B));
    _NC_max = ceil(_K_co/maxnorm(_C));
    
    // SET-UP POLAR GROUNDS (FORE-, MID-, BACK-)
    _fg_C.clear();
    _fg_N.clear();
    _mg_N.clear();
    _bg_N.clear();
    _bg_P.clear();

    _fg_C = ptop->FGC();
    _fg_N = ptop->FGN();
    _bg_N = ptop->BGN();        
    _bg_P.insert(_bg_P.end(), _fg_N.begin(), _fg_N.end());
    _bg_P.insert(_bg_P.end(), _bg_N.begin(), _bg_N.end());    
    _inForeground.resize(_bg_P.size()+1,false);
    
    // SET-UP MIDGROUND (INCLUDING PERIODIC IMAGES IF REQUIRED)
    LOG(logINFO,*_log) << flush;
    LOG(logINFO,*_log) << "Generate periodic images. ";
    this->SetupMidground(_R_co);
    
    // CALCULATE COG POSITIONS, NET CHARGE; SET-UP BOOLEAN FG TABLE
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    double Q_fg_C = 0.0;
    double Q_fg_N = 0.0;
    double Q_mg_N = 0.0;
    double Q_bg_N = 0.0;
    double Q_bg_P = 0.0;  
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_C += (*sit)->CalcTotQ();
        _inForeground[(*sit)->getId()] = true;
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_fg_N += (*sit)->CalcTotQ();
        assert(_inForeground[(*sit)->getId()] == true);
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_mg_N += (*sit)->CalcTotQ();
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_N += (*sit)->CalcTotQ();
        assert(_inForeground[(*sit)->getId()] == false);
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        (*sit)->CalcPos();
        Q_bg_P += (*sit)->CalcTotQ();
    }
    
    LOG(logINFO,*_log)
        << (format("Net ground charge and size")).str()
        << flush << (format("  o Q(FGC) = %1$+1.3fe |FGC| = %2$+5d") % Q_fg_C % _fg_C.size()).str()
        << flush << (format("  o Q(FGN) = %1$+1.3fe |FGN| = %2$+5d") % Q_fg_N % _fg_N.size()).str()
        << flush << (format("  o Q(MGN) = %1$+1.3fe |MGN| = %2$+5d") % Q_mg_N % _mg_N.size()).str()
        << flush << (format("  o Q(BGN) = %1$+1.3fe |BGN| = %2$+5d") % Q_bg_N % _bg_N.size()).str()
        << flush << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P % _bg_P.size()).str()
        << flush;
    
    if (std::abs(Q_bg_P) > 1e-2) {
        cout << endl;
        cout << endl << format("***************************** ERROR ******************************");
        cout << endl << format("       Background charge |Q(BGP)| is larger than 0.01e.");
        cout << endl << format("       Ignore: e.g. rounding error?");
        cout << endl << format("       Or think again: e.g. erroneous parametrization?");
        cout << endl << format("******************************************************************");
        cout << endl;
    }
    
    // CHARGE APPROPRIATELY & DEPOLARIZE
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _fg_N.begin(); sit < _fg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _bg_N.begin(); sit < _bg_N.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            (*pit)->Depolarize();
        }
    }
    
    // CALCULATE NET DIPOLE OF BGP & FGC
    vec netdpl_bgP = vec(0,0,0);
    vec netdpl_fgC = vec(0,0,0);
    double qzz_bgP = 0.0;
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_bgP += (*pit)->getPos() * (*pit)->getQ00();
            qzz_bgP += (*pit)->getQ00() * ((*pit)->getPos().getZ() * (*pit)->getPos().getZ());
        }
    }
    for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {        
        PolarSeg* pseg = *sit;        
        for (pit = pseg->begin(); pit < pseg->end(); ++pit) {
            netdpl_fgC += (*pit)->getPos() * (*pit)->getQ00();
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
    
    return;
}


void Ewald3D::SetupMidground(double R_co) {
    // SET-UP MIDGROUND
    // TODO Extend this to several molecules in the foreground
    assert(_fg_C.size() == 1);
    _center = _fg_C[0]->getPos();
    // NOTE No periodic-boundary correction here: We require that all
    //      polar-segment CoM coords. be folded with respect to the central
    //      segment
    // NOTE Excludes interaction of polar segment with neutral self in real-
    //      space sum
    // NOTE Includes periodic images if within cut-off
    
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;
    for (sit = _mg_N.begin(); sit < _mg_N.end(); ++sit)
        delete (*sit);
    _mg_N.clear();
    
    // IMAGE BOXES TO CONSIDER
    vector< int > nas;
    vector< int > nbs;
    vector< int > ncs;
    for (int na = -_na_max; na < _na_max+1; ++na)
        nas.push_back(na);
    for (int nb = -_nb_max; nb < _nb_max+1; ++nb)
        nbs.push_back(nb);
    for (int nc = -_nc_max; nc < _nc_max+1; ++nc)
        ncs.push_back(nc);
    vector< int >::iterator ait;
    vector< int >::iterator bit;
    vector< int >::iterator cit;
    
    // SAMPLE MIDGROUND FROM BGP EXCLUDING CENTRAL SEG.
    assert(_fg_N.size() == _fg_C.size());
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        PolarSeg *pseg = *sit;
        Segment *other = _top->getSegment(pseg->getId());
        // Periodic images
        for (ait = nas.begin(); ait < nas.end(); ++ait) {
            int na = *ait;
            for (bit = nbs.begin(); bit < nbs.end(); ++bit) {                
                int nb = *bit;
                for (cit = ncs.begin(); cit < ncs.end(); ++cit) {
                    int nc = *cit;
                    if (na == 0 && nb == 0 && nc == 0 && _inForeground[other->getId()])
                        continue;
                    vec L = na*_a + nb*_b + nc*_c;
                    vec pos_L = pseg->getPos() + L;
                    double dR_L = abs(_center-pos_L);
                    if (dR_L <= R_co) {
                        PolarSeg *newSeg = new PolarSeg(pseg);
                        newSeg->Translate(L);
                        _mg_N.push_back(newSeg);
                    }
                } // Loop over nc
            } // Loop over nb
        } // Loop over na
    } // Loop over BGP
    return;
}


void Ewald3D::WriteDensitiesPDB(string pdbfile) {
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


double Ewald3D::ConvergeRealSpaceSum() {
    
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
    for (int i = 0; i < 100; ++i) {
        double Rc = _R_co + (i-1)*dR;
        // Set-up midground
        this->SetupMidground(Rc);
        // Calculate interaction energy
        this_ER = 0.;
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        this_ER += _actor.E_QQ_ERFC(*(*pit1), *(*pit2), _alpha);
                    }
                }
            }
        }
        double dER_rms = sqrt((this_ER-prev_ER)*(this_ER-prev_ER))*_actor.int2eV;
        LOG(logDEBUG,*_log)
            << (format("Rc = %1$+1.7f   |MGN| = %3$5d nm   ER = %2$+1.7f eV   dER(rms) = %4$+1.7f") 
            % Rc % (this_ER*_actor.int2eV) % _mg_N.size() % dER_rms).str() << flush;
        if (i > 0 && dER_rms < _crit_dE) {
            _converged_R = true;
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision as of Rc = %1$+1.3f nm") 
                % Rc ) << flush;
            break;
        }
        prev_ER = this_ER;
    }
    return this_ER;
}


double Ewald3D::ConvergeReciprocalSpaceSum() {
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;    
    
    // CELLS OF THE RECIPROCAL LATTICE: SUM OVER CUBOIDAL SHELLS
    // ... Shell-size increment is magnitude of largest reciprocal cell vector
    // ... Tschebyschow norm is used to group k-vectors into k-shells 
    
    vector< vector<vec> > shell_ks;
    vector< vector<vec> >::iterator shellit;
    double shell_dk = (maxnorm(_A) > maxnorm(_B)) ?
        ((maxnorm(_A) > maxnorm(_C)) ? maxnorm(_A) : maxnorm(_C)) 
      : ((maxnorm(_B) > maxnorm(_C)) ? maxnorm(_B) : maxnorm(_C));
    // Determine all k-vectors within k-space cut-off
    vector< vec >::iterator kit;
    vector< vec > ks;
    for (int kx = -_NA_max; kx < _NA_max+1; ++kx) {
        for (int ky = -_NB_max; ky < _NB_max+1; ++ky) {
            for (int kz = -_NC_max; kz < _NC_max+1; ++kz) {
                if (kx == 0 && ky == 0 && kz == 0) continue;
                vec k = kx*_A + ky*_B + kz*_C;
                if (maxnorm(k) > _K_co) continue;
                ks.push_back(k);
            }
        }
    }
    // Sort according to magnitude
    std::sort(ks.begin(), ks.end(), _maxsort);
    // Group into shells
    int shell_idx = 0;
    kit = ks.begin();
    shell_ks.resize(int(_K_co/shell_dk+0.5)+1);
    double shell_k = shell_dk;
    while (shell_k <= _K_co) {
        for ( ; kit < ks.end(); ++kit) {
            vec k = *kit;
            if (maxnorm(k) <= shell_k) {
                // Add to current shell
                shell_ks[shell_idx].push_back(k);
            }
            else {
                // Open new shell
                shell_k += shell_dk;
                shell_idx += 1;
                shell_ks[shell_idx].push_back(k);
                break;
            }
        }
        // All k-vectors consumed?
        if (kit == ks.end()) break;
    }
//    for (int i = 0; i < shell_ks.size(); ++i) {
//        ofstream ofs;
//        string outfile = (format("shell_%1$d.out") % (i+1)).str();
//        ofs.open(outfile.c_str(), ofstream::out);
//        for (kit = shell_ks[i].begin(); kit < shell_ks[i].end(); ++kit) {
//            ofs << (*kit).getX() << " " << (*kit).getY() << " " << (*kit).getZ() << endl;
//        }
//        ofs.close();
//    }
     
        
    // K=K TERM
    LOG(logDEBUG,*_log) << flush;
    double EKK_fgC_bgP = 0.0;    
    int N_EKK_memory = int(0.5*(_NA_max+_NB_max)+0.5);
    int N_K_proc = 0;
    int N_shells_proc = 0;
    vector< double > dEKKs;
    _converged_K = false;
    
    double re_E = 0.0;
    double im_E = 0.0;
    
    for (shellit = shell_ks.begin(); shellit < shell_ks.end(); ++shellit, ++N_shells_proc) {
        
        for (kit = (*shellit).begin(); kit < (*shellit).end(); ++kit, ++N_K_proc) {
            vec k = *kit;

            // K-DEPENDENT FACTOR
            double K = abs(k);
            double expkk_k = 4*M_PI*exp(-K*K/(4*_alpha*_alpha)) / (K*K);        

            LOG(logDEBUG,*_log)
                << (format("k[%5$d] = %1$+1.3f %2$+1.3f %3$+1.3f   |K| = %4$+1.3f 1/nm") 
                % (k.getX()) % (k.getY()) % (k.getZ()) % K % (N_shells_proc+1));

            // STRUCTURE FACTORS
            // Calculate structure factor S(k) for FGC
            double qcos_fgC = 0.0;
            double qsin_fgC = 0.0;
            for (sit = _fg_C.begin(); sit < _fg_C.end(); ++sit) {
                for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
                    qcos_fgC += (*pit)->Q00 * cos(k * (*pit)->getPos());
                    qsin_fgC += (*pit)->Q00 * sin(k * (*pit)->getPos());
                }
            }        
            // Calculate structure factor S(-k) for BGP
            double qcos_bgP = 0.0;
            double qsin_bgP = 0.0;
            for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
                for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
                    qcos_bgP += (*pit)->Q00 * cos(-k * (*pit)->getPos());
                    qsin_bgP += (*pit)->Q00 * sin(-k * (*pit)->getPos());
                }
            }
            // Structure-factor product
            double re_s1s2 = qcos_fgC*qcos_bgP - qsin_fgC*qsin_bgP;
            double im_s1s2 = qsin_fgC*qcos_bgP + qcos_fgC*qsin_bgP;

            // REAL & IMAGINARY ENERGY
            double re_dE = expkk_k * re_s1s2;
            double im_dE = expkk_k * im_s1s2;        
            re_E += re_dE;
            im_E += im_dE;

            LOG(logDEBUG,*_log)
                << (format("    Re(dE) = %1$+1.7f") 
                % (re_dE/_LxLyLz*_actor.int2eV));

            LOG(logDEBUG,*_log)
                << (format("    Re(E) = %1$+1.7f Im(E) = %2$+1.7f")
                % (re_E/_LxLyLz*_actor.int2eV) % (im_E/_LxLyLz*_actor.int2eV));        

            // CONVERGED?
            double dEKK = sqrt(re_dE*re_dE + im_dE*im_dE);
            double dEKK_rms = 0.0;
            if (dEKKs.size() < N_EKK_memory) {
                dEKKs.resize(N_EKK_memory,dEKK);
            }
            else {
                dEKKs[N_K_proc % N_EKK_memory] = dEKK;
            }
            for (int i = 0; i < dEKKs.size(); ++i) {
                dEKK_rms += dEKKs[i]*dEKKs[i];
            }
            dEKK_rms /= dEKKs.size();
            dEKK_rms = sqrt(dEKK_rms);

            LOG(logDEBUG,*_log)
                << (format("   RMS(%2$d) = %1$+1.7f") 
                % (dEKK_rms/_LxLyLz*_actor.int2eV) % N_EKK_memory) << flush;

            if (dEKK_rms/_LxLyLz*_actor.int2eV <= _crit_dE && N_K_proc > 2 && N_shells_proc > 0) {
                _converged_K = true;
                LOG(logDEBUG,*_log)
                    << (format(":::: Converged to precision as of |K| = %1$+1.3f 1/nm") 
                    % K ) << flush;
                break;
            }
        } // Sum over k's in k-shell
        if (_converged_K) break;
    } // Sum over k-shells
    
    EKK_fgC_bgP = re_E/_LxLyLz;
    return EKK_fgC_bgP;
}


double Ewald3D::CalculateShapeCorrection(string shape) {
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    double EJ = 0.0;
    
    if (shape == "xyslab") {
        // DIRECT CALCULATION VIA DOUBLE LOOP
        // TODO The double-loop can be avoided, but direct summation as below 
        //      appears to be more stable from a numerical point of view
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
           for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
              for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                 for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    double za = (*pit1)->getPos().getZ();
                    double zb = (*pit2)->getPos().getZ();
                    EJ += (za-zb)*(za-zb) * (*pit1)->getQ00()*(*pit2)->getQ00();
                 }
              }
           }
        }

        //    // DECOMPOSITION INTO CHARGE-DENSITY MULTIPOLE MOMENTS
        //    vec DA = vec(0,0,0);
        //    vec DA = vec(0,0,0);
        //    double QA = 0.0;
        //    double DZZBG = 0.0;
        //    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        //        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        //            DA += (*pit1)->getQ00() * (*pit1)->getPos();
        //            QA += (*pit1)->getQ00();
        //        }
        //    }
        //    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        //        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        //            double zb = (*pit1)->getPos().getZ();
        //            DB += (*pit1)->getQ00() * (*pit1)->getPos();
        //            DZZB += (*pit1)->getQ00() * zb*zb;
        //        }
        //    }
        //    EJ = QA*DZZB - 2*(DA.getZ())*(DB.getZ());

        EJ *= - 2*M_PI/_LxLyLz;
    }
    else {
        LOG(logERROR,*_log)
            << (format("Shape %1$s not implemented. Setting EJ = 0.0 ...") 
            % shape) << flush;
        EJ = 0.0;
    }
        
    return EJ;
}


double Ewald3D::CalculateSq2(vec &k) {
    vector<PolarSeg*>::iterator sit; 
    vector<APolarSite*> ::iterator pit;    
    double cs = 0.0;
    double ss = 0.0;    
    for (sit = _bg_P.begin(); sit < _bg_P.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            cs += (*pit)->Q00 * cos(k * (*pit)->getPos());
            ss += (*pit)->Q00 * sin(k * (*pit)->getPos());
        }
    }    
    return cs*cs + ss*ss;
}


void Ewald3D::Evaluate() {
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logDEBUG,*_log) << "System & Ewald parameters (3D x 3D)" << flush;
    LOG(logDEBUG,*_log) << "  o Real-space unit cell:      " << _a << " x " << _b << " x " << _c << flush;
    LOG(logDEBUG,*_log) << "  o Real-space c/o (guess):    " << _R_co << " nm" << flush;
    LOG(logDEBUG,*_log) << "  o na(max), nb(max), nc(max): " << _na_max << ", " << _nb_max << ", " << _nc_max << flush;
    LOG(logDEBUG,*_log) << "  o 1st Brillouin zone:        " << _A << " x " << _B << " x " << _C << flush;
    LOG(logDEBUG,*_log) << "  o Reciprocal-space c/o:      " << _K_co << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o R-K switching param.       " << _alpha << " 1/nm" << flush;
    LOG(logDEBUG,*_log) << "  o Unit-cell volume:          " << _LxLyLz << " nm**3" << flush;
    LOG(logDEBUG,*_log) << "  o LxLy (3D, hence omitted):  " << _LxLy << " nm**2" << flush;
    LOG(logDEBUG,*_log) << "  o kx(max), ky(max), kz(max): " << _NA_max << ", " << _NB_max << ", " << _NC_max << flush;
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;    
    
    // REAL-SPACE CONTRIBUTION
    double EPP_fgC_mgN = ConvergeRealSpaceSum();    
    
    // RECIPROCAL-SPACE CONTRIBUTION
    double EKK_fgC_bgP = ConvergeReciprocalSpaceSum();
    
    // SHAPE-CORRECTION (FROM K=0 TERM)
    double EJ_fgC_bgP = CalculateShapeCorrection("xyslab");
    
    // FOREGROUND CORRECTION
    double EPP_fgC_fgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_fgN += _actor.E_QQ_ERF(*(*pit1), *(*pit2), _alpha);
                }
            }
        }
    }
    
    // REAL-SPACE HIGHER-RANK CORRECTION
    double EDQ_fgC_mgN = 0.0;
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    _actor.BiasStat(*(*pit1), *(*pit2));
                    EDQ_fgC_mgN += _actor.E_Q0_DQ(*(*pit1), *(*pit2));
                }
            }
        }
    }
    
    _ER = EPP_fgC_mgN*_actor.int2eV;
    _EC = EPP_fgC_fgN*_actor.int2eV;
    _EK = EKK_fgC_bgP*_actor.int2eV;
    _EJ = EJ_fgC_bgP*_actor.int2eV;
    _E0 = 0.0; // Contained in shape-correction _EJ
    _ET = _ER + _EK + _E0 - _EC;
    _EDQ = EDQ_fgC_mgN*_actor.int2eV;
    
    LOG(logDEBUG,*_log) << flush;
    LOG(logINFO,*_log)
        << (format("Interaction FGC -> ***")).str()
        << flush << (format("  + EPP(FGC->MGN)  = %1$+1.7f eV") % _ER).str()
        << flush << (format("  + EK0(FGC->BGP)  = %1$+1.7f eV") % _E0).str()
        << flush << (format("  + EKK(FGC->BGP)  = %1$+1.7f eV") % _EK).str()
        << flush << (format("  - EPP(FGC->FGN)  = %1$+1.7f eV") % _EC).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("    SUM(E)         = %1$+1.7f eV") % _ET).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + EDQ(FGC->MGN)  = %1$+1.7f eV") % _EDQ).str()
        << flush << (format("  + EJ (%2$s)    = %1$+1.7f eV") % _EJ % _shape).str()
        << flush << (format("    ------------------------------")).str()
        << flush << (format("  + SUM(E) (+DQ,J) = %1$+1.7f eV") % (_ET+_EJ+_EDQ)).str()
        << flush;
    LOG(logDEBUG,*_log) << flush;    
    return;
}


string Ewald3D::GenerateErrorString() {
    string rstr;
    rstr += (format("Converged R-sum = %1$s, converged K-sum = %2$s")
        % ((_converged_R) ? "true" : "false")
        % ((_converged_K) ? "true" : "false")).str();
    return rstr;
}


string Ewald3D::GenerateOutputString() {
    string rstr;
    rstr += (format("XYZ %1$+1.7f %2$+1.7f %3$+1.7f ") 
        % _center.getX() % _center.getY() % _center.getZ()).str();
    rstr += (format("ET %1$+1.7f ER %2$+1.7f EK %3$+1.7f E0 %4$+1.7f EC %5$+1.7f EDQ %6$+1.7f ") 
        % _ET % _ER % _EK % _E0 % _EC %_EDQ).str();
    rstr += (format("FGC %1$1d FGN %2$1d MGN %3$3d BGN %4$4d BGP %5$4d") 
        % _fg_C.size() % _fg_N.size() % _mg_N.size() % _bg_N.size() 
        % _bg_P.size()).str();
    return rstr;
}
    
    
    
    
    
    
}}