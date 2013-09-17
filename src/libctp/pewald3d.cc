#include <votca/ctp/pewald3d.h>
#include <boost/format.hpp>
#include <algorithm>
#include <boost/date_time/posix_time/posix_time.hpp>


using boost::format;


namespace votca { namespace ctp {


PEwald3D3D::~PEwald3D3D() { ; }
    
    
PEwald3D3D::PEwald3D3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
  : Ewald3DnD(top, ptop, opt, log) {
    _shape = opt->get("options.ewald.coulombmethod.shape").as<string>();
}


double PEwald3D3D::ConvergeRealSpaceSum() {
    
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
        
        
        
        //boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
        //string pdbfile = (format("%1$s") % now.time_of_day()).str();        
        //this->WriteDensitiesPDB(pdbfile);
        
        // Calculate interaction energy
        this_ER = 0.;        

        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
            for (sit2 = _mg_N.begin(); sit2 < _mg_N.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                        this_ER += _ewdactor.U12_ERFC(*(*pit1), *(*pit2));
                    }
                }
                //cout << endl << "I " << (*sit1)->getId() << " --- " << (*sit2)->getId() << " " << showpos << scientific << this_ER*_ewdactor.int2eV << flush;
                //break;
            }
        }
        
        double dER_rms = sqrt((this_ER-prev_ER)*(this_ER-prev_ER))*_ewdactor.int2eV;
        LOG(logDEBUG,*_log)
            << (format("Rc = %1$+1.7f   |MGN| = %3$5d nm   ER = %2$+1.7f eV   dER(rms) = %4$+1.7f") 
            % Rc % (this_ER*_ewdactor.int2eV) % _mg_N.size() % dER_rms).str() << flush;
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


double PEwald3D3D::ConvergeReciprocalSpaceSum() {
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;    
    
    // CELLS OF THE RECIPROCAL LATTICE: SUM OVER ELLIPSOIDAL SHELLS
    // ... Shell-size increment is magnitude of largest reciprocal cell vector
    // ... Tschebyschow/Euclidean norm is used to group k-vectors into k-shells 
    
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
    std::sort(ks.begin(), ks.end(), _eucsort);
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
        
        double this_shell_dE = 0.0;
        double this_shell_maxK = 0.0;
        
        for (kit = (*shellit).begin(); kit < (*shellit).end(); ++kit, ++N_K_proc) {
            vec k = *kit;
            double K = abs(k);
            double Kxy = sqrt(k.getX()*k.getX()+k.getY()*k.getY());
            double Kz = sqrt(k.getZ()*k.getZ());
            
            LOG(logDEBUG,*_log)
                << (format("k[%5$d] = %1$+1.3f %2$+1.3f %3$+1.3f   |Kxy| = %4$+1.3f 1/nm") 
                % (k.getX()) % (k.getY()) % (k.getZ()) % (Kxy/K) % (N_shells_proc+1));
            
            EwdInteractor::cmplx as1s2 = _ewdactor.AS1S2(k, _fg_C, _bg_P);
            

            
            

            // REAL & IMAGINARY ENERGY
            double re_dE = as1s2._re;
            double im_dE = as1s2._im;
            re_E += re_dE;
            im_E += im_dE;            
            this_shell_dE += re_dE;
            this_shell_maxK = K;

            LOG(logDEBUG,*_log)
                << (format("    Re(dE) = %1$+1.7f")
                % (re_dE/_LxLyLz*_ewdactor.int2eV));

            LOG(logDEBUG,*_log)
                << (format("    Re(E) = %1$+1.7f Im(E) = %2$+1.7f")
                % (re_E/_LxLyLz*_ewdactor.int2eV) % (im_E/_LxLyLz*_ewdactor.int2eV));        

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
                % (dEKK_rms/_LxLyLz*_ewdactor.int2eV) % N_EKK_memory) << flush;

            if (false && dEKK_rms/_LxLyLz*_ewdactor.int2eV <= _crit_dE && N_K_proc > 2 && N_shells_proc > 0) {
                _converged_K = true;
                LOG(logDEBUG,*_log)
                    << (format(":::: Converged to precision as of |K| = %1$+1.3f 1/nm") 
                    % K ) << flush;
                break;
            }
        } // Sum over k's in k-shell
        if (_converged_K) break;
        
        double abs_dE = sqrt(this_shell_dE*this_shell_dE);
        LOG(logDEBUG,*_log)
            << (format("===> SHELL    maxK = %1$+1.7f   |dE| = %2$+1.7f") 
                % this_shell_maxK % (this_shell_dE/_LxLyLz*_ewdactor.int2eV)) << flush;
        if (abs_dE/_LxLyLz*_ewdactor.int2eV <= _crit_dE) {
            _converged_K = true;
            LOG(logDEBUG,*_log)
                << (format(":::: SHELL Converged to precision as of |K| = %1$+1.3f 1/nm") 
                % this_shell_maxK ) << flush;
            break;
        }
        
    } // Sum over k-shells
    
    EKK_fgC_bgP = re_E/_LxLyLz;
    return EKK_fgC_bgP;
}


double PEwald3D3D::CalculateShapeCorrection() {
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    double EJ = 0.0;
    
    if (_shape == "xyslab") {
        for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
           for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
              for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                 for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EJ += _ewdactor.U12_XYSlab(*(*pit1), *(*pit2));
                 }
              }
           }
        }
        EJ *= - 2*M_PI/_LxLyLz;
    }
    else {
        LOG(logERROR,*_log)
            << (format("Shape %1$s not implemented. Setting EJ = 0.0 ...") 
            % _shape) << flush;
        EJ = 0.0;
    }
        
    return EJ;
}


double PEwald3D3D::CalculateForegroundCorrection() {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EPP_fgC_fgN = 0.0;
    
    for (sit1 = _fg_C.begin(); sit1 < _fg_C.end(); ++sit1) {
        for (sit2 = _fg_N.begin(); sit2 < _fg_N.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EPP_fgC_fgN += _ewdactor.U12_ERF(*(*pit1), *(*pit2));
                }
            }
        }
    }
    
    return EPP_fgC_fgN;
}
    
    
}}