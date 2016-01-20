#include <votca/xtp/ewald3d.h>
#include <boost/format.hpp>
#include <algorithm>


using boost::format;


namespace votca { namespace xtp {


Ewald3D3D::~Ewald3D3D() { ; }
    
    
Ewald3D3D::Ewald3D3D(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
  : Ewald3DnD(top, ptop, opt, log) {}


EWD::triple<> Ewald3D3D::ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target) {

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
    unsigned int N_EKK_memory = (unsigned int)(0.5*(_NA_max+_NB_max)+0.5);
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
            for (sit = target.begin(); sit < target.end(); ++sit) {
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
                % (re_dE/_LxLyLz*int2eV));

            LOG(logDEBUG,*_log)
                << (format("    Re(E) = %1$+1.7f Im(E) = %2$+1.7f")
                % (re_E/_LxLyLz*int2eV) % (im_E/_LxLyLz*int2eV));        

            // CONVERGED?
            double dEKK = sqrt(re_dE*re_dE + im_dE*im_dE);
            double dEKK_rms = 0.0;
            if (dEKKs.size() < N_EKK_memory) {
                dEKKs.resize(N_EKK_memory,dEKK);
            }
            else {
                dEKKs[N_K_proc % N_EKK_memory] = dEKK;
            }
            for (unsigned int i = 0; i < dEKKs.size(); ++i) {
                dEKK_rms += dEKKs[i]*dEKKs[i];
            }
            dEKK_rms /= dEKKs.size();
            dEKK_rms = sqrt(dEKK_rms);

            LOG(logDEBUG,*_log)
                << (format("   RMS(%2$d) = %1$+1.7f") 
                % (dEKK_rms/_LxLyLz*int2eV) % N_EKK_memory) << flush;

            if (dEKK_rms/_LxLyLz*int2eV <= _crit_dE && N_K_proc > 2 && N_shells_proc > 0) {
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
    return EWD::triple<>(EKK_fgC_bgP,0,0);
}


EWD::triple<> Ewald3D3D::CalculateShapeCorrection(vector<PolarSeg*> &target) {
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    
    double EJ = 0.0;
    
    if (_shape == "xyslab") {
        // DIRECT CALCULATION VIA DOUBLE LOOP
        // TODO The double-loop can be avoided, but direct summation as below 
        //      appears to be more stable from a numerical point of view
        for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
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
        //    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
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
            % _shape) << flush;
        EJ = 0.0;
    }
        
    return EWD::triple<>(EJ,0,0);
}


double Ewald3D3D::CalculateSq2(vec &k) {
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
    
    
}}
