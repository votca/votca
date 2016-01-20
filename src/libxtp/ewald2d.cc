#include <votca/xtp/ewald2d.h>
#include <boost/format.hpp>
#include <algorithm>


using boost::format;


namespace votca { namespace xtp {


Ewald3D2D::~Ewald3D2D() { ; }
    
    
Ewald3D2D::Ewald3D2D(Topology *top, PolarTop *ptop, Property *opt, Logger *log) 
  : Ewald3DnD(top, ptop, opt, log) {
    _nc_max = 0;
}


EWD::triple<> Ewald3D2D::ConvergeReciprocalSpaceSum(vector<PolarSeg*> &target) {
    
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;    
    
    // CELLS OF THE RECIPROCAL LATTICE OVER WHICH TO SUM
    //    ---------------------------------
    //      | x | x | x | x |   |   |   |
    //    ---------------------------------            
    //      | x | x | x | x |   |   |   |
    //    ---------------------------------
    //      | x | x | x | O |   |   |   |
    //    ---------------------------------
    //      | x | x | x |   |   |   |   |  
    //    ---------------------------------
    //      | x | x | x |   |   |   |   |  
    //    ---------------------------------
    
    vector< vec > Ks;
    vector< vec >::iterator kit;
    //int kz = 0; // This is 2D Ewald
    int kx = 0;
    int ky = 0;
    for (ky = 1; ky < _NB_max+1; ++ky) {        
        vec K = kx*_A + ky*_B;
        if (abs(K)>_K_co) { continue; }
        Ks.push_back(K);
    }
    for (kx = 1; kx < _NA_max+1; ++kx) {
        for (ky = -_NB_max; ky < _NB_max+1; ++ky) {
            vec K = kx*_A + ky*_B;            
            if (abs(K)>_K_co) { continue; }            
            Ks.push_back(K);
        }
    }
    std::sort(Ks.begin(), Ks.end(), _eucsort);  
        
    // K=K TERM
    LOG(logDEBUG,*_log) << flush;
    double EKK_fgC_bgP = 0.0;    
    unsigned N_EKK_memory = unsigned(0.5*(_NA_max+_NB_max)+0.5);
    int N_K_proc = 0;
    vector< double > dEKKs;
    _converged_K = false;
    
    for (kit = Ks.begin(); kit < Ks.end(); ++kit, ++N_K_proc) {
        vec k = *kit;
        double K = abs(k);
        double dEKK = 0.0;
        LOG(logDEBUG,*_log)  
            << (format("k = %1$+1.3f %2$+1.3f %3$+1.3f   |K| = %4$+1.3f 1/nm") 
            % (k.getX()) % (k.getY()) % (k.getZ()) % K);
        for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
            for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
                for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                    for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {                                                
                        dEKK += _actor.E_QQ_KK(*(*pit1), *(*pit2), _alpha, k);
                    }
                }
            }
        }
        
        EKK_fgC_bgP += dEKK;
        
        // CONVERGED?
        double dEKK_rms = 0.0;
        if (dEKKs.size() < N_EKK_memory) {
            dEKKs.resize(N_EKK_memory,dEKK);
        }
        else {
            dEKKs[N_K_proc % N_EKK_memory] = dEKK;
        }
        for (unsigned i = 0; i < dEKKs.size(); ++i) {
            dEKK_rms += dEKKs[i]*dEKKs[i];
        }
        dEKK_rms /= dEKKs.size();
        dEKK_rms = sqrt(dEKK_rms);
        
        LOG(logDEBUG,*_log)  
            << (format("   dE = %1$+1.7f   EKK = %2$+1.7f   RMS(%4$d) = %3$+1.7f") 
            % (dEKK*2*M_PI/_LxLy*int2eV) % (EKK_fgC_bgP*2*M_PI/_LxLy*int2eV)
            % (dEKK_rms*2*M_PI/_LxLy*int2eV) % N_EKK_memory) << flush;
        
//        cout << endl;
//        for (int i = 0; i < dEKKs.size(); ++i) {
//            cout << "| " << dEKKs[i];
//        }
//        cout << flush;
        
        if (dEKK_rms*2*M_PI/_LxLy*int2eV <= _crit_dE && N_K_proc > 10) {
        //if (dEKK_rms*2*M_PI/LxLy*int2eV <= dEKK_rms_crit) {
            _converged_K = true;
            LOG(logDEBUG,*_log)  
                << (format(":::: Converged to precision as of |K| = %1$+1.3f 1/nm") 
                % K ) << flush;
            break;
        }
        
        
    }
    EKK_fgC_bgP *= 2*M_PI/_LxLy;
    return EWD::triple<>(EKK_fgC_bgP,0,0);    
}


EWD::triple<> Ewald3D2D::CalculateK0Correction(vector<PolarSeg*> &target) {
    vector<PolarSeg*>::iterator sit1; 
    vector<APolarSite*> ::iterator pit1;
    vector<PolarSeg*>::iterator sit2; 
    vector<APolarSite*> ::iterator pit2;
    double EK0_fgC_bgP = 0.0;
    for (sit1 = target.begin(); sit1 < target.end(); ++sit1) {
        for (sit2 = _bg_P.begin(); sit2 < _bg_P.end(); ++sit2) {
            for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
                for (pit2 = (*sit2)->begin(); pit2 < (*sit2)->end(); ++pit2) {
                    EK0_fgC_bgP += _actor.E_QQ_K0(*(*pit1), *(*pit2), _alpha);
                }
            }
        }
    }
    EK0_fgC_bgP *= 2*M_PI/_LxLy;
    return EWD::triple<>(EK0_fgC_bgP,0,0);
}
    
    
}}