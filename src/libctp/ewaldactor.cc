#include <votca/ctp/ewaldactor.h>


namespace votca {
namespace ctp {

    
void EwdInteractor::FPU12_XYSlab_ShapeField_At_By(vector<PolarSeg*> &at, 
    vector<PolarSeg*> &by, double &TwoPi_V) {
    // This function requires neutrality of &s but gets around
    // the double sum (see ::F12_XYSlab_At_By)
    // ATTENTION Only increments FPz, not FUz
    double fz = 0.0;
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    for (sit = by.begin(); sit < by.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            fz += 2*( (*pit)->Q00*(*pit)->getPos().getX() + (*pit)->U1z);
            if ((*pit)->_rank > 0) 
                fz += 2*(*pit)->Q1z;
        }
    }
    fz *= TwoPi_V;
    // Increment fields
    for (sit = at.begin(); sit < at.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->FPz += fz;
        }
    }
    return;
}


void EwdInteractor::FP12_XYSlab_ShapeField_At_By(vector<PolarSeg*> &at, 
    vector<PolarSeg*> &by, double &TwoPi_V) {
    // This function requires neutrality of &s but gets around
    // the double sum (see ::F12_XYSlab_At_By)
    double fz = 0.0;
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    for (sit = by.begin(); sit < by.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            fz += 2*( (*pit)->Q00*(*pit)->getPos().getX() );
            if ((*pit)->_rank > 0) 
                fz += 2*(*pit)->Q1z;
        }
    }
    fz *= TwoPi_V;
    // Increment fields
    for (sit = at.begin(); sit < at.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->FPz += fz;
        }
    }
    return;
}


void EwdInteractor::FU12_XYSlab_ShapeField_At_By(vector<PolarSeg*> &at, 
    vector<PolarSeg*> &by, double &TwoPi_V) {
    // This function requires neutrality of &s but gets around
    // the double sum (see ::F12_XYSlab_At_By)
    double fz = 0.0;
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    for (sit = by.begin(); sit < by.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            fz += 2*(*pit)->U1z;
        }
    }
    fz *= TwoPi_V;
    // Increment fields
    for (sit = at.begin(); sit < at.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->FUz += fz;
        }
    }
    return;
}



// ============================ RECIPROCAL SPACE ============================ //
//                                 S-FACTORS                                  //

EWD::cmplx EwdInteractor::PUStructureAmplitude(vector<PolarSeg*> &s) {
    double re_S = 0.0;
    double im_S = 0.0;
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    for (sit = s.begin(); sit < s.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            PUApplyBiasK(*(*pit));
            re_S += re_s;
            re_S += u_re_s;
            im_S += im_s;
            im_S += u_im_s;
        }
    }
    return EWD::cmplx(re_S, im_S);
}


EWD::cmplx EwdInteractor::PStructureAmplitude(vector<PolarSeg*> &s) {
    double re_S = 0.0;
    double im_S = 0.0;
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    for (sit = s.begin(); sit < s.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            PApplyBiasK(*(*pit));
            re_S += re_s;
            im_S += im_s;
        }
    }
    return EWD::cmplx(re_S, im_S);
}


EWD::cmplx EwdInteractor::UStructureAmplitude(vector<PolarSeg*> &s) {
    double re_S = 0.0;
    double im_S = 0.0;
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    for (sit = s.begin(); sit < s.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            UApplyBiasK(*(*pit));
            re_S += u_re_s;
            im_S += u_im_s;
        }
    }
    return EWD::cmplx(re_S, im_S);
}


// ============================ RECIPROCAL SPACE ============================ //
//                                   FIELD                                    //


EWD::cmplx EwdInteractor::FPU12_AS1S2_At_By(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increment PERMANENT fields of s1
    // ATTENTION Structure factors include PERMANENT & INDUCED moments of s2    
    
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    double sum_re_f_rms = 0.0;
    double sum_im_f_xyz = 0.0;
    int rms_count = 0;
        
    // Structure amplitude S2* from s2 = B*
    EWD::cmplx cmplx_S2 = PUStructureAmplitude(s2).Conjugate();
    double re_S2 = cmplx_S2._re;
    double im_S2 = cmplx_S2._im;
    
    // Compute k-component of fields acting on s1 = A(c)
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            kr = kx * (*pit)->getPos().getX()
               + ky * (*pit)->getPos().getY()
               + kz * (*pit)->getPos().getZ();
            coskr = cos(kr);
            sinkr = sin(kr);
            
            // Real component
            double fx = -rV*AK * kx * (sinkr*re_S2 + coskr*im_S2);
            double fy = -rV*AK * ky * (sinkr*re_S2 + coskr*im_S2);
            double fz = -rV*AK * kz * (sinkr*re_S2 + coskr*im_S2);
            
            (*pit)->FPx += fx;
            (*pit)->FPy += fy;
            (*pit)->FPz += fz;
            
            // Imaginary component (error check)
            double ifx = -rV*AK * kx * (sinkr*im_S2 - coskr*re_S2);
            double ify = -rV*AK * ky * (sinkr*im_S2 - coskr*re_S2);
            double ifz = -rV*AK * kz * (sinkr*im_S2 - coskr*re_S2);
            
            rms_count += 1;
            sum_re_f_rms += fx*fx + fy*fy + fz*fz;
            sum_im_f_xyz += ifx + ify + ifz;            
        }
    }
    
    sum_re_f_rms /= rms_count;
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    return EWD::cmplx(sum_re_f_rms, sum_im_f_xyz);    
}


EWD::cmplx EwdInteractor::FP12_AS1S2_At_By(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increment PERMANENT fields of s1
    // ATTENTION Structure factors include PERMANENT moments of s2    
    
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    double sum_re_f_rms = 0.0;
    double sum_im_f_xyz = 0.0;
    int rms_count = 0;
        
    // Structure amplitude S2* from s2 = B*
    EWD::cmplx cmplx_S2 = PStructureAmplitude(s2).Conjugate();
    double re_S2 = cmplx_S2._re;
    double im_S2 = cmplx_S2._im;
    
    // Compute k-component of fields acting on s1 = A(c)
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            kr = kx * (*pit)->getPos().getX()
               + ky * (*pit)->getPos().getY()
               + kz * (*pit)->getPos().getZ();
            coskr = cos(kr);
            sinkr = sin(kr);
            
            // Real component
            double fx = -rV*AK * kx * (sinkr*re_S2 + coskr*im_S2);
            double fy = -rV*AK * ky * (sinkr*re_S2 + coskr*im_S2);
            double fz = -rV*AK * kz * (sinkr*re_S2 + coskr*im_S2);
            
            (*pit)->FPx += fx;
            (*pit)->FPy += fy;
            (*pit)->FPz += fz;
            
            // Imaginary component (error check)
            double ifx = -rV*AK * kx * (sinkr*im_S2 - coskr*re_S2);
            double ify = -rV*AK * ky * (sinkr*im_S2 - coskr*re_S2);
            double ifz = -rV*AK * kz * (sinkr*im_S2 - coskr*re_S2);
            
            rms_count += 1;
            sum_re_f_rms += fx*fx + fy*fy + fz*fz;
            sum_im_f_xyz += ifx + ify + ifz;            
        }
    }
    
    sum_re_f_rms /= rms_count;
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    return EWD::cmplx(sum_re_f_rms, sum_im_f_xyz);    
}


EWD::cmplx EwdInteractor::FU12_AS1S2_At_By(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increment INDUCDED fields of s1
    // ATTENTION Structure factors include INDUCED moments of s2    
    
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    double sum_re_f_rms = 0.0;
    double sum_im_f_xyz = 0.0;
    int rms_count = 0;
        
    // Structure amplitude S2* from s2 = B*
    EWD::cmplx cmplx_S2 = UStructureAmplitude(s2).Conjugate();
    double re_S2 = cmplx_S2._re;
    double im_S2 = cmplx_S2._im;
    
    // Compute k-component of fields acting on s1 = A(c)
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            kr = kx * (*pit)->getPos().getX()
               + ky * (*pit)->getPos().getY()
               + kz * (*pit)->getPos().getZ();
            coskr = cos(kr);
            sinkr = sin(kr);
            
            // Real component
            double fx = -rV*AK * kx * (sinkr*re_S2 + coskr*im_S2);
            double fy = -rV*AK * ky * (sinkr*re_S2 + coskr*im_S2);
            double fz = -rV*AK * kz * (sinkr*re_S2 + coskr*im_S2);
            
            (*pit)->FUx += fx;
            (*pit)->FUy += fy;
            (*pit)->FUz += fz;
            
            // Imaginary component (error check)
            double ifx = -rV*AK * kx * (sinkr*im_S2 - coskr*re_S2);
            double ify = -rV*AK * ky * (sinkr*im_S2 - coskr*re_S2);
            double ifz = -rV*AK * kz * (sinkr*im_S2 - coskr*re_S2);
            
            rms_count += 1;
            sum_re_f_rms += fx*fx + fy*fy + fz*fz;
            sum_im_f_xyz += ifx + ify + ifz;            
        }
    }
    
    sum_re_f_rms /= rms_count;
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    return EWD::cmplx(sum_re_f_rms, sum_im_f_xyz);    
}


// ============================ RECIPROCAL SPACE ============================ //
//                                  ENERGIES                                  //

EWD::triple<EWD::cmplx> EwdInteractor::AS1S2(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2) {
    // NOTE : w/o 1/V
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // Structure amplitude S1
    double re_S1 = 0.0;
    double im_S1 = 0.0;
    double u_re_S1 = 0.0;
    double u_im_S1 = 0.0;
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            PUApplyBiasK(*(*pit));            
            re_S1 += re_s;
            im_S1 += im_s;     // NOTE THE (+)
            u_re_S1 += u_re_s;
            u_im_S1 += u_im_s; // NOTE THE (+)
        }
    }    
    
    // Structure amplitude S2
    double re_S2 = 0.0;
    double im_S2 = 0.0;
    double u_re_S2 = 0.0;
    double u_im_S2 = 0.0;
    for (sit = s2.begin(); sit < s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            PUApplyBiasK(*(*pit));
            re_S2 += re_s;
            im_S2 -= im_s;     // NOTE THE (-)
            u_re_S2 += u_re_s;
            u_im_S2 -= u_im_s; // NOTE THE (-)
        }
    }
    
    double pp_re_AS1S2 = AK * (re_S1*re_S2 - im_S1*im_S2);
    double pp_im_AS1S2 = AK * (re_S1*im_S2 + im_S1*re_S2);
    
    double uu_re_AS1S2 = AK * (u_re_S1*u_re_S2 - u_im_S1*u_im_S2);
    double uu_im_AS1S2 = AK * (u_re_S1*u_im_S2 + u_im_S1*u_re_S2);
    
    double pu_re_AS1S2 = AK * (u_re_S1*re_S2 + re_S1*u_re_S2 - u_im_S1*im_S2 - im_S1*u_im_S2);
    double pu_im_AS1S2 = AK * (u_re_S1*im_S2 + re_S1*u_im_S2 + u_im_S1*re_S2 + im_S1*u_re_S2);
    
    
    return EWD::triple<EWD::cmplx>(EWD::cmplx(pp_re_AS1S2, pp_im_AS1S2),
                                   EWD::cmplx(pu_re_AS1S2, pu_im_AS1S2),
                                   EWD::cmplx(uu_re_AS1S2, uu_im_AS1S2));
}


EWD::triple<EWD::cmplx> EwdInteractor::S1S2(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2) {
    // NOTE : w/o 1/V
    ApplyBiasK(k);    
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // Structure amplitude S1
    double re_S1 = 0.0;
    double im_S1 = 0.0;   
    double u_re_S1 = 0.0;
    double u_im_S1 = 0.0;
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            PUApplyBiasK(*(*pit));            
            re_S1 += re_s;
            im_S1 += im_s;     // NOTE THE (+)
            u_re_S1 += u_re_s;
            u_im_S1 += u_im_s; // NOTE THE (+)
        }
    }    
    
    // Structure amplitude S2
    double re_S2 = 0.0;
    double im_S2 = 0.0;
    double u_re_S2 = 0.0;
    double u_im_S2 = 0.0;
    for (sit = s2.begin(); sit < s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {            
            PUApplyBiasK(*(*pit));
            re_S2 += re_s;
            im_S2 -= im_s;     // NOTE THE (-)
            u_re_S2 += u_re_s;
            u_im_S2 -= u_im_s; // NOTE THE (-)
        }
    }
    
    double pp_re_S1S2 = (re_S1*re_S2 - im_S1*im_S2);
    double pp_im_S1S2 = (re_S1*im_S2 + im_S1*re_S2);
    
    double uu_re_S1S2 = (u_re_S1*u_re_S2 - u_im_S1*u_im_S2);
    double uu_im_S1S2 = (u_re_S1*u_im_S2 + u_im_S1*u_re_S2);
    
    double pu_re_S1S2 = (u_re_S1*re_S2 + re_S1*u_re_S2 - u_im_S1*im_S2 - im_S1*u_im_S2);
    double pu_im_S1S2 = (u_re_S1*im_S2 + re_S1*u_im_S2 + u_im_S1*re_S2 + im_S1*u_re_S2);
    
    return EWD::triple<EWD::cmplx>(EWD::cmplx(pp_re_S1S2, pp_im_S1S2),
                                   EWD::cmplx(pu_re_S1S2, pu_im_S1S2),
                                   EWD::cmplx(uu_re_S1S2, uu_im_S1S2));
}
    




}}