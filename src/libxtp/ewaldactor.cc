#include <votca/xtp/ewaldactor.h>
#include <boost/format.hpp>


namespace votca {
namespace xtp {

    
void EwdInteractor::FPU12_ShapeField_At_By(vector<PolarSeg*> &at,
    vector<PolarSeg*> &by, string shape, double V) {
    // This function requires neutrality of &by but gets around
    // the double sum (see ::F12_XYSlab_At_By)
    // ATTENTION Only increments FP, not FU
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Compute net dipole moment of field-generating density (&by)
    vec Q1 = vec(0,0,0);
    vec U1 = vec(0,0,0);
    for (sit = by.begin(); sit < by.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            vec r = (*pit)->getPos();
            double q0 = (*pit)->Q00;
            vec q1 = vec((*pit)->Q1x, (*pit)->Q1y, (*pit)->Q1z);
            vec u1 = vec((*pit)->U1x, (*pit)->U1y, (*pit)->U1z);
            // Increment moments
            Q1 += q0*r;
            if ((*pit)->getRank() > 0)
                Q1 += q1;
            U1 += u1;
        }
    }
    
    // Shape-dependent fields
    if (shape == "xyslab") {
        fx = 0.0;
        fy = 0.0;
        fz = Q1.getZ() + U1.getZ();
        fx *= 4*M_PI/V;
        fy *= 4*M_PI/V;
        fz *= 4*M_PI/V;
    }
    else if (shape == "cube" || shape == "sphere") {
        fx = Q1.getX() + U1.getX();
        fy = Q1.getY() + U1.getY();
        fz = Q1.getZ() + U1.getZ();
        fx *= 4*M_PI/(3*V);
        fy *= 4*M_PI/(3*V);
        fz *= 4*M_PI/(3*V);
    }
    else if (shape == "none") {
    	fx = 0.0;
    	fy = 0.0;
    	fz = 0.0;
    }
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    // Increment fields of target density (&at)
    for (sit = at.begin(); sit < at.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->FPx += fx;
            (*pit)->FPy += fy;
            (*pit)->FPz += fz;
        }
    }
    return;
}


void EwdInteractor::FP12_ShapeField_At_By(vector<PolarSeg*> &at,
    vector<PolarSeg*> &by, string shape, double V) {
    // This function requires neutrality of &by but gets around
    // the double sum (see ::F12_XYSlab_At_By)
    // ATTENTION Only increments FP, not FU
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Compute net dipole moment of field-generating density (&by)
    vec Q1 = vec(0,0,0);
    for (sit = by.begin(); sit < by.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            vec r = (*pit)->getPos();
            double q0 = (*pit)->Q00;
            vec q1 = vec((*pit)->Q1x, (*pit)->Q1y, (*pit)->Q1z);
            // Increment moments
            Q1 += q0*r;
            if ((*pit)->getRank() > 0)
                Q1 += q1;
        }
    }
    
    // Shape-dependent fields
    if (shape == "xyslab") {
        fx = 0.0;
        fy = 0.0;
        fz = Q1.getZ();
        fx *= 4*M_PI/V;
        fy *= 4*M_PI/V;
        fz *= 4*M_PI/V;
    }
    else if (shape == "cube" || shape == "sphere") {
        fx = Q1.getX();
        fy = Q1.getY();
        fz = Q1.getZ();
        fx *= 4*M_PI/(3*V);
        fy *= 4*M_PI/(3*V);
        fz *= 4*M_PI/(3*V);
    }
    else if (shape == "none") {
    	fx = 0.0;
    	fy = 0.0;
    	fz = 0.0;
    }
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    // Increment fields of target density (&at)
    for (sit = at.begin(); sit < at.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->FPx += fx;
            (*pit)->FPy += fy;
            (*pit)->FPz += fz;
        }
    }
    return;
}


void EwdInteractor::FU12_ShapeField_At_By(vector<PolarSeg*> &at,
    vector<PolarSeg*> &by, string shape, double V) {
    // This function requires neutrality of &by but gets around
    // the double sum (see ::F12_XYSlab_At_By)
    // ATTENTION Only increments FP, not FU
    
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Compute net dipole moment of field-generating density (&by)
    vec U1 = vec(0,0,0);
    for (sit = by.begin(); sit < by.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            vec u1 = vec((*pit)->U1x, (*pit)->U1y, (*pit)->U1z);
            // Increment moments
            U1 += u1;
        }
    }
    
    // Shape-dependent fields
    if (shape == "xyslab") {
        fx = 0.0;
        fy = 0.0;
        fz = U1.getZ();
        fx *= 4*M_PI/V;
        fy *= 4*M_PI/V;
        fz *= 4*M_PI/V;
    }
    else if (shape == "cube" || shape == "sphere") {
        fx = U1.getX();
        fy = U1.getY();
        fz = U1.getZ();
        fx *= 4*M_PI/(3*V);
        fy *= 4*M_PI/(3*V);
        fz *= 4*M_PI/(3*V);
    }
    else if (shape == "none") {
		fx = 0.0;
		fy = 0.0;
		fz = 0.0;
	}
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    // Increment fields of target density (&at)
    for (sit = at.begin(); sit < at.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            (*pit)->FUx += fx;
            (*pit)->FUy += fy;
            (*pit)->FUz += fz;
        }
    }
    return;
}
    

void EwdInteractor::PhiPU12_ShapeField_At_By(vector<PolarSeg*> &s1, 
    vector<PolarSeg*> &s2, string shape, double V) {
    // NOTE : WITH PREFACTOR = -4*PI/V (xyslab) v -4*PI/3V (cube, sphere)    
    // NOTE : Similar operations in PhiPU12_..., PhiP12_..., PhiU12_...
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Charge, dipole, quadrupole for <s2>:
    // Q... <> permanent, U... <> induced
    double Q0_S2 = 0.0;
    votca::tools::vec Q1_S2 = vec(0,0,0);
    votca::tools::vec U1_S2 = vec(0,0,0);
    votca::tools::matrix Q2_S2;
    Q2_S2.ZeroMatrix();
    votca::tools::matrix U2_S2;
    U2_S2.ZeroMatrix();
    
    for (sit = s2.begin(); sit != s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
            vec r = (*pit)->getPos();
            double q0 = (*pit)->Q00;
            vec    q1 = vec((*pit)->Q1x, (*pit)->Q1y, (*pit)->Q1z);
            vec    u1 = vec((*pit)->U1x, (*pit)->U1y, (*pit)->U1z);
            matrix q2 = matrix(
                vec((*pit)->Qxx, (*pit)->Qxy, (*pit)->Qxz),
                vec((*pit)->Qxy, (*pit)->Qyy, (*pit)->Qyz),
                vec((*pit)->Qxz, (*pit)->Qyz, (*pit)->Qzz));
            // Charge
            Q0_S2 += q0;
            // Dipole
            Q1_S2 += q0*r;
            if ((*pit)->getRank() > 0)
                Q1_S2 += q1;
            U1_S2 += u1;
            // Quadrupole
            Q2_S2 += 0.5*q0*(r|r);
            if ((*pit)->getRank() > 0) {
                Q2_S2 += (q1|r);
                if ((*pit)->getRank() > 1) {
                    Q2_S2 += q2;
                }
            }
            U2_S2 += (u1|r);
        }
    }
    
    // Traces
    double TrQ2_S2 = Q2_S2.get(0,0)+Q2_S2.get(1,1)+Q2_S2.get(2,2);
    double TrU2_S2 = U2_S2.get(0,0)+U2_S2.get(1,1)+U2_S2.get(2,2);
    
    
    // Increment potentials
    if (shape == "xyslab") {
        double prefac = -4*M_PI/V;
        for (sit = s1.begin(); sit != s1.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {           
                (*pit)->PhiP += prefac*(Q2_S2.get(2,2) - (*pit)->getPos().getZ()*Q1_S2.getZ());
                (*pit)->PhiU += prefac*(U2_S2.get(2,2) - (*pit)->getPos().getZ()*U1_S2.getZ());  
            }
        }
    }
    else if (shape == "cube" || shape == "sphere") {
        double prefac = -4*M_PI/(3*V);
        for (sit = s1.begin(); sit != s1.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->PhiP += prefac*(TrQ2_S2 - (*pit)->getPos()*Q1_S2);
                (*pit)->PhiU += prefac*(TrU2_S2 - (*pit)->getPos()*U1_S2);
            }
        }
    }
    else if (shape == "none") {
    	; // Nothing to do here
    }
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    return;
}


void EwdInteractor::PhiP12_ShapeField_At_By(vector<PolarSeg*> &s1, 
    vector<PolarSeg*> &s2, string shape, double V) {
    // NOTE : WITH PREFACTOR = -4*PI/V (xyslab) v -4*PI/3V (cube, sphere)
    // NOTE : Similar operations in PhiPU12_..., PhiP12_..., PhiU12_...
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Charge, dipole, quadrupole for <s2>:
    // Q... <> permanent, U... <> induced
    double Q0_S2 = 0.0;
    votca::tools::vec Q1_S2 = vec(0,0,0);
    votca::tools::matrix Q2_S2;
    Q2_S2.ZeroMatrix();
    
    for (sit = s2.begin(); sit != s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
            vec r = (*pit)->getPos();
            double q0 = (*pit)->Q00;
            vec    q1 = vec((*pit)->Q1x, (*pit)->Q1y, (*pit)->Q1z);
            matrix q2 = matrix(
                vec((*pit)->Qxx, (*pit)->Qxy, (*pit)->Qxz),
                vec((*pit)->Qxy, (*pit)->Qyy, (*pit)->Qyz),
                vec((*pit)->Qxz, (*pit)->Qyz, (*pit)->Qzz));
            // Charge
            Q0_S2 += q0;
            // Dipole
            Q1_S2 += q0*r;
            if ((*pit)->getRank() > 0)
                Q1_S2 += q1;
            // Quadrupole
            Q2_S2 += 0.5*q0*(r|r);
            if ((*pit)->getRank() > 0) {
                Q2_S2 += (q1|r);
                if ((*pit)->getRank() > 1) {
                    Q2_S2 += q2;
                }
            }
        }
    }
    
    // Traces
    double TrQ2_S2 = Q2_S2.get(0,0)+Q2_S2.get(1,1)+Q2_S2.get(2,2);    
    
    // Increment potentials
    if (shape == "xyslab") {
        double prefac = -4*M_PI/V;
        for (sit = s1.begin(); sit != s1.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {           
                (*pit)->PhiP += prefac*(Q2_S2.get(2,2) - (*pit)->getPos().getZ()*Q1_S2.getZ());
            }
        }
    }
    else if (shape == "cube" || shape == "sphere") {
        double prefac = -4*M_PI/(3*V);
        for (sit = s1.begin(); sit != s1.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->PhiP += prefac*(TrQ2_S2 - (*pit)->getPos()*Q1_S2);
            }
        }
    }
    else if (shape == "none") {
    	; // Nothing to do here
    }
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    return;
}


void EwdInteractor::PhiU12_ShapeField_At_By(vector<PolarSeg*> &s1, 
    vector<PolarSeg*> &s2, string shape, double V) {
    // NOTE : WITH PREFACTOR = -4*PI/V (xyslab) v -4*PI/3V (cube, sphere)
    // NOTE : Similar operations in PhiPU12_..., PhiP12_..., PhiU12_...
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Charge, dipole, quadrupole for <s2>:
    // Q... <> permanent, U... <> induced
    votca::tools::vec U1_S2 = vec(0,0,0);
    votca::tools::matrix U2_S2;
    U2_S2.ZeroMatrix();
    
    for (sit = s2.begin(); sit != s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
            vec r = (*pit)->getPos();
            vec    u1 = vec((*pit)->U1x, (*pit)->U1y, (*pit)->U1z);
            // Dipole
            U1_S2 += u1;
            // Quadrupole
            U2_S2 += (u1|r);
        }
    }
    
    // Traces
    double TrU2_S2 = U2_S2.get(0,0)+U2_S2.get(1,1)+U2_S2.get(2,2);    
    
    // Increment potentials
    if (shape == "xyslab") {
        double prefac = -4*M_PI/V;
        for (sit = s1.begin(); sit != s1.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {           
                (*pit)->PhiU += prefac*(U2_S2.get(2,2) - (*pit)->getPos().getZ()*U1_S2.getZ());  
            }
        }
    }
    else if (shape == "cube" || shape == "sphere") {
        double prefac = -4*M_PI/(3*V);
        for (sit = s1.begin(); sit != s1.end(); ++sit) {
            for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
                (*pit)->PhiU += prefac*(TrU2_S2 - (*pit)->getPos()*U1_S2);
            }
        }
    }
    else if (shape == "none") {
    	; // Nothing to do here
    }
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    return;
}


//void EwdInteractor::FPU12_XYSlab_ShapeField_At_By(vector<PolarSeg*> &at, 
//    vector<PolarSeg*> &by, double &TwoPi_V) {
//    // This function requires neutrality of &by but gets around
//    // the double sum (see ::F12_XYSlab_At_By)
//    // ATTENTION Only increments FPz, not FUz
//    double fz = 0.0;
//    vector<PolarSeg*>::iterator sit;
//    vector<APolarSite*> ::iterator pit;
//    for (sit = by.begin(); sit < by.end(); ++sit) {
//        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
//            fz += 2*( (*pit)->Q00*(*pit)->getPos().getZ() + (*pit)->U1z);
//            if ((*pit)->_rank > 0) 
//                fz += 2*(*pit)->Q1z;
//        }
//    }
//    fz *= TwoPi_V;
//    // Increment fields
//    for (sit = at.begin(); sit < at.end(); ++sit) {
//        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
//            (*pit)->FPz += fz;
//        }
//    }
//    return;
//}
//
//
//void EwdInteractor::FP12_XYSlab_ShapeField_At_By(vector<PolarSeg*> &at, 
//    vector<PolarSeg*> &by, double &TwoPi_V) {
//    // This function requires neutrality of &by but gets around
//    // the double sum (see ::F12_XYSlab_At_By)
//    double fz = 0.0;
//    vector<PolarSeg*>::iterator sit;
//    vector<APolarSite*> ::iterator pit;
//    for (sit = by.begin(); sit < by.end(); ++sit) {
//        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
//            fz += 2*( (*pit)->Q00*(*pit)->getPos().getZ() );
//            if ((*pit)->_rank > 0) 
//                fz += 2*(*pit)->Q1z;
//        }
//    }
//    fz *= TwoPi_V;
//    // Increment fields
//    for (sit = at.begin(); sit < at.end(); ++sit) {
//        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
//            (*pit)->FPz += fz;
//        }
//    }
//    return;
//}
//
//
//void EwdInteractor::FU12_XYSlab_ShapeField_At_By(vector<PolarSeg*> &at, 
//    vector<PolarSeg*> &by, double &TwoPi_V) {
//    // This function requires neutrality of &by but gets around
//    // the double sum (see ::F12_XYSlab_At_By)
//    double fz = 0.0;
//    vector<PolarSeg*>::iterator sit;
//    vector<APolarSite*> ::iterator pit;
//    for (sit = by.begin(); sit < by.end(); ++sit) {
//        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
//            fz += 2*(*pit)->U1z;
//        }
//    }
//    fz *= TwoPi_V;
//    // Increment fields
//    for (sit = at.begin(); sit < at.end(); ++sit) {
//        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
//            (*pit)->FUz += fz;
//        }
//    }
//    return;
//}



// ============================ RECIPROCAL SPACE ============================ //
//                                 S-FACTORS                                  //

EWD::cmplx EwdInteractor::PUStructureAmplitude(vector<PolarSeg*> &s) {
    // Requires ApplyBias(k) to be called previously
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
    // Requires ApplyBias(k) to be called previously
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
    // Requires ApplyBias(k) to be called previously
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


EWD::cmplx EwdInteractor::PUStructureAmplitude(vector<PolarSeg*> &s, const vec &k) {
    ApplyBiasK(k);
    return PUStructureAmplitude(s);
}


EWD::cmplx EwdInteractor::PStructureAmplitude(vector<PolarSeg*> &s, const vec &k) {
    ApplyBiasK(k);
    return PStructureAmplitude(s);
}


EWD::cmplx EwdInteractor::UStructureAmplitude(vector<PolarSeg*> &s, const vec &k) {
    ApplyBiasK(k);
    return UStructureAmplitude(s);
}


// ============================ RECIPROCAL SPACE ============================ //
//                                   FIELD                                    //


EWD::cmplx EwdInteractor::FP12_At_ByS2(const vec &k, vector<PolarSeg*> &s1, 
    const EWD::cmplx &S2, double &rV) {
    // ATTENTION Increments PERMANENT fields of s1
    // ATTENTION Structure factor S2 from PERM & INDU moments of s2, k    
    double re_S2 = S2._re;
    double im_S2 = - S2._im; // !! NOTE THE (-) !!
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    double sum_re_f_rms = 0.0;
    double sum_im_f_xyz = 0.0;
    int rms_count = 0;
    
    // Compute k-component of fields acting on s1 = A(c)
    ApplyBiasK(k);
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;    
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


EWD::cmplx EwdInteractor::FU12_At_ByS2(const vec &k, vector<PolarSeg*> &s1, 
    const EWD::cmplx &S2, double &rV) {
    // ATTENTION Increments PERMANENT fields of s1
    // ATTENTION Structure factor S2 from PERM & INDU moments of s2, k    
    double re_S2 = S2._re;
    double im_S2 = - S2._im; // !! NOTE THE (-) !!
    
    // NOTE sum_re_f_rms => convergence check (to be performed by caller)
    // NOTE sum_im_f_xyz => sanity check      (to be performed by caller)
    double sum_re_f_rms = 0.0;
    double sum_im_f_xyz = 0.0;
    int rms_count = 0;
    
    // Compute k-component of fields acting on s1 = A(c)
    ApplyBiasK(k);
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;    
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


EWD::cmplx EwdInteractor::FPU12_AS1S2_At_By(const vec &k,
    vector<PolarSeg*> &s1, vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increments PERMANENT fields of s1
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
    // ATTENTION Increments PERMANENT fields of s1
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
    // ATTENTION Increments INDUCED fields of s1
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
//                                 POTENTIALS                                 //

EWD::cmplx EwdInteractor::PhiPU12_AS1S2_At_By(const vec &k, vector<PolarSeg*> &s1, 
    vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increments *permanent* potential PhiP of s1 only
    // ATTENTION Structure factors include PERMANENT & INDUCED moments of s2
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    
    ApplyBiasK(k);
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // NOTE sum_re_phi_rms => convergence check (to be performed by caller)
    // NOTE sum_im_phi_xyz => sanity check      (to be performed by caller)
    double sum_re_phi_ms = 0.0;
    double sum_im_phi = 0.0;
    int rms_count = 0;
        
    // Structure amplitude S2* from s2 = B*
    EWD::cmplx cmplx_S2 = PUStructureAmplitude(s2).Conjugate();
    double re_S2 = cmplx_S2._re;
    double im_S2 = cmplx_S2._im;
    
    // Compute k-component of potential acting on s1 = A(c)
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            kr = kx * (*pit)->getPos().getX()
               + ky * (*pit)->getPos().getY()
               + kz * (*pit)->getPos().getZ();
            coskr = cos(kr);
            sinkr = sin(kr);
            
            // Real component
            double phi  = rV*AK * (coskr*re_S2 - sinkr*im_S2);
            (*pit)->PhiP += phi;
            
            // Imaginary component (error check)
            double iphi = rV*AK * (coskr*im_S2 + sinkr*re_S2);
            
            rms_count += 1;
            sum_re_phi_ms += phi*phi;
            sum_im_phi += iphi;
        }
    }
    
    sum_re_phi_ms /= rms_count;
    
    // NOTE sum_re_phi_rms => convergence check (to be performed by caller)
    // NOTE sum_im_phi     => sanity check      (to be performed by caller)
    return EWD::cmplx(sum_re_phi_ms, sum_im_phi);
}


EWD::cmplx EwdInteractor::PhiP12_AS1S2_At_By(const vec &k, vector<PolarSeg*> &s1, 
    vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increments *permanent* potential PhiP of s1 only
    // ATTENTION Structure factors include PERMANENT moments of s2 only
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    
    ApplyBiasK(k);
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // NOTE sum_re_phi_rms => convergence check (to be performed by caller)
    // NOTE sum_im_phi_xyz => sanity check      (to be performed by caller)
    double sum_re_phi_ms = 0.0;
    double sum_im_phi = 0.0;
    int rms_count = 0;
        
    // Structure amplitude S2* from s2 = B*
    EWD::cmplx cmplx_S2 = PStructureAmplitude(s2).Conjugate();
    double re_S2 = cmplx_S2._re;
    double im_S2 = cmplx_S2._im;
    
    // Compute k-component of potential acting on s1 = A(c)
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            kr = kx * (*pit)->getPos().getX()
               + ky * (*pit)->getPos().getY()
               + kz * (*pit)->getPos().getZ();
            coskr = cos(kr);
            sinkr = sin(kr);
            
            // Real component
            double phi  = rV*AK * (coskr*re_S2 - sinkr*im_S2);
            (*pit)->PhiP += phi;
            
            // Imaginary component (error check)
            double iphi = rV*AK * (coskr*im_S2 + sinkr*re_S2);
            
            rms_count += 1;
            sum_re_phi_ms += phi*phi;
            sum_im_phi += iphi;
        }
    }
    
    sum_re_phi_ms /= rms_count;
    
    // NOTE sum_re_phi_rms => convergence check (to be performed by caller)
    // NOTE sum_im_phi     => sanity check      (to be performed by caller)
    return EWD::cmplx(sum_re_phi_ms, sum_im_phi);
}


EWD::cmplx EwdInteractor::PhiU12_AS1S2_At_By(const vec &k, vector<PolarSeg*> &s1, 
    vector<PolarSeg*> &s2, double &rV) {
    // ATTENTION Increments *permanent* potential PhiP of s1 only
    // ATTENTION Structure factors include PERMANENT moments of s2 only
    // NOTE Analogous operations in PhiPU12_..., PhiP12_..., PhiU12_...
    
    ApplyBiasK(k);
    
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*> ::iterator pit;
    
    // NOTE sum_re_phi_rms => convergence check (to be performed by caller)
    // NOTE sum_im_phi_xyz => sanity check      (to be performed by caller)
    double sum_re_phi_ms = 0.0;
    double sum_im_phi = 0.0;
    int rms_count = 0;
        
    // Structure amplitude S2* from s2 = B*
    EWD::cmplx cmplx_S2 = UStructureAmplitude(s2).Conjugate();
    double re_S2 = cmplx_S2._re;
    double im_S2 = cmplx_S2._im;
    
    // Compute k-component of potential acting on s1 = A(c)
    for (sit = s1.begin(); sit < s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit < (*sit)->end(); ++pit) {
            kr = kx * (*pit)->getPos().getX()
               + ky * (*pit)->getPos().getY()
               + kz * (*pit)->getPos().getZ();
            coskr = cos(kr);
            sinkr = sin(kr);
            
            // Real component
            double phi  = rV*AK * (coskr*re_S2 - sinkr*im_S2);
            (*pit)->PhiU += phi;
            
            // Imaginary component (error check)
            double iphi = rV*AK * (coskr*im_S2 + sinkr*re_S2);
            
            rms_count += 1;
            sum_re_phi_ms += phi*phi;
            sum_im_phi += iphi;
        }
    }
    
    sum_re_phi_ms /= rms_count;
    
    // NOTE sum_re_phi_rms => convergence check (to be performed by caller)
    // NOTE sum_im_phi     => sanity check      (to be performed by caller)
    return EWD::cmplx(sum_re_phi_ms, sum_im_phi);
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


EWD::triple<double> EwdInteractor::U12_ShapeTerm(vector<PolarSeg*> &s1,
    vector<PolarSeg*> &s2, string shape, double V, Logger *log) {
    
    // NOTE : WITH PREFACTOR = -4*PI/V (xylab) v -4*PI/3V (cube, sphere)
    vector<PolarSeg*>::iterator sit;
    vector<APolarSite*>::iterator pit;
    
    // Charge, dipole, quadrupole for <s1>: Q... <> permanent, U... <> induced
    double Q0_S1 = 0.0;
    votca::tools::vec Q1_S1 = vec(0,0,0);
    votca::tools::vec U1_S1 = vec(0,0,0);
    votca::tools::matrix Q2_S1;
    Q2_S1.ZeroMatrix();
    votca::tools::matrix U2_S1;
    U2_S1.ZeroMatrix();
    
    for (sit = s1.begin(); sit != s1.end(); ++sit) {
        for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
            vec    r  = (*pit)->getPos();
            double q0 = (*pit)->Q00;
            vec    q1 = vec((*pit)->Q1x, (*pit)->Q1y, (*pit)->Q1z);
            vec    u1 = vec((*pit)->U1x, (*pit)->U1y, (*pit)->U1z);
            matrix q2 = matrix(
                vec((*pit)->Qxx, (*pit)->Qxy, (*pit)->Qxz),
                vec((*pit)->Qxy, (*pit)->Qyy, (*pit)->Qyz),
                vec((*pit)->Qxz, (*pit)->Qyz, (*pit)->Qzz));
            // Charge
            Q0_S1 += q0;
            // Dipole
            Q1_S1 += q0*r;
            if ((*pit)->getRank() > 0)
                Q1_S1 += q1;
            U1_S1 += u1;
            // Quadrupole
            Q2_S1 += 0.5*q0*(r|r);
            if ((*pit)->getRank() > 0) {
                Q2_S1 += (q1|r);
                if ((*pit)->getRank() > 1) {
                    Q2_S1 += q2;
                }
            }
            U2_S1 += (u1|r);
        }
    }    
    
    // Charge, dipole, quadrupole for <s2>: Q... <> permanent, U... <> induced
    double Q0_S2 = 0.0;
    votca::tools::vec Q1_S2 = vec(0,0,0);
    votca::tools::vec U1_S2 = vec(0,0,0);
    votca::tools::matrix Q2_S2;
    Q2_S2.ZeroMatrix();
    votca::tools::matrix U2_S2;
    U2_S2.ZeroMatrix();
    
    for (sit = s2.begin(); sit != s2.end(); ++sit) {
        for (pit = (*sit)->begin(); pit != (*sit)->end(); ++pit) {
            vec r = (*pit)->getPos();
            double q0 = (*pit)->Q00;
            vec    q1 = vec((*pit)->Q1x, (*pit)->Q1y, (*pit)->Q1z);
            vec    u1 = vec((*pit)->U1x, (*pit)->U1y, (*pit)->U1z);
            matrix q2 = matrix(
                vec((*pit)->Qxx, (*pit)->Qxy, (*pit)->Qxz),
                vec((*pit)->Qxy, (*pit)->Qyy, (*pit)->Qyz),
                vec((*pit)->Qxz, (*pit)->Qyz, (*pit)->Qzz));
            // Charge
            Q0_S2 += q0;
            // Dipole
            Q1_S2 += q0*r;
            if ((*pit)->getRank() > 0)
                Q1_S2 += q1;
            U1_S2 += u1;
            // Quadrupole
            Q2_S2 += 0.5*q0*(r|r);
            if ((*pit)->getRank() > 0) {
                Q2_S2 += (q1|r);
                if ((*pit)->getRank() > 1) {
                    Q2_S2 += q2;
                }
            }
            U2_S2 += (u1|r);
        }
    }
    
    // Traces
    double TrQ2_S1 = Q2_S1.get(0,0)+Q2_S1.get(1,1)+Q2_S1.get(2,2);
    double TrU2_S1 = U2_S1.get(0,0)+U2_S1.get(1,1)+U2_S1.get(2,2);    
    double TrQ2_S2 = Q2_S2.get(0,0)+Q2_S2.get(1,1)+Q2_S2.get(2,2);
    double TrU2_S2 = U2_S2.get(0,0)+U2_S2.get(1,1)+U2_S2.get(2,2);    
    
    if (log != NULL) {
        LOG(logDEBUG, *log) << "S1 moments: " << flush;
        LOG(logDEBUG, *log) << "  Q0   = " 
            << (boost::format("%1$+1.7e") % Q0_S1) << flush;
        LOG(logDEBUG, *log) << "  Q1   = " 
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % Q1_S1.getX() % Q1_S1.getY() % Q1_S1.getZ()) << flush;
        LOG(logDEBUG, *log) << "  U1   = " 
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U1_S1.getX() % U1_S1.getY() % U1_S1.getZ()) << flush;
        
        LOG(logDEBUG, *log) << "  Q2   = "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % Q2_S1.get(0,0) % Q2_S1.get(0,1) % Q2_S1.get(0,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % Q2_S1.get(1,0) % Q2_S1.get(1,1) % Q2_S1.get(1,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e")
                % Q2_S1.get(2,0) % Q2_S1.get(2,1) % Q2_S1.get(2,2)) << flush;
        
        LOG(logDEBUG, *log) << "  U2   = "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U2_S1.get(0,0) % U2_S1.get(0,1) % U2_S1.get(0,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U2_S1.get(1,0) % U2_S1.get(1,1) % U2_S1.get(1,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U2_S1.get(2,0) % U2_S1.get(2,1) % U2_S1.get(2,2)) << flush;
        
        LOG(logDEBUG, *log) << "  TrQ2 = " << TrQ2_S1 << flush;
        LOG(logDEBUG, *log) << "  TrU2 = " << TrU2_S1 << flush;
        
        
        LOG(logDEBUG, *log) << "S2 moments: " << flush;
        LOG(logDEBUG, *log) << "  Q0   = " 
            << (boost::format("%1$+1.7e") % Q0_S2) << flush;
        LOG(logDEBUG, *log) << "  Q1   = " 
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % Q1_S2.getX() % Q1_S2.getY() % Q1_S2.getZ()) << flush;
        LOG(logDEBUG, *log) << "  U1   = " 
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U1_S2.getX() % U1_S2.getY() % U1_S2.getZ()) << flush;
        
        LOG(logDEBUG, *log) << "  Q2   = "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % Q2_S2.get(0,0) % Q2_S2.get(0,1) % Q2_S2.get(0,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % Q2_S2.get(1,0) % Q2_S2.get(1,1) % Q2_S2.get(1,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e")
                % Q2_S2.get(2,0) % Q2_S2.get(2,1) % Q2_S2.get(2,2)) << flush;
        
        LOG(logDEBUG, *log) << "  U2   = "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U2_S2.get(0,0) % U2_S2.get(0,1) % U2_S2.get(0,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U2_S2.get(1,0) % U2_S2.get(1,1) % U2_S2.get(1,2)) << flush;
        LOG(logDEBUG, *log) << "         "
            << (boost::format("%1$+1.7e %2$+1.7e %3$+1.7e") 
                % U2_S2.get(2,0) % U2_S2.get(2,1) % U2_S2.get(2,2)) << flush;
        LOG(logDEBUG, *log) << "  TrQ2 = " << TrQ2_S2 << flush;
        LOG(logDEBUG, *log) << "  TrU2 = " << TrU2_S2 << flush;
    }
    
    // Shape-dependent energies
    double pp = 0.0;
    double pu = 0.0;
    double uu = 0.0;
    
    if (shape == "xyslab") {
        pp = Q0_S1*Q2_S2.get(2,2)
           + Q0_S2*Q2_S1.get(2,2)
           - Q1_S1.getZ()*Q1_S2.getZ();
        pu = Q0_S1*U2_S2.get(2,2)
           + Q0_S2*U2_S1.get(2,2)
           - Q1_S1.getZ()*U1_S2.getZ()
           - Q1_S2.getZ()*U1_S1.getZ();
        uu = - U1_S1.getZ()*U1_S2.getZ();
        pp *= -4*M_PI/V;
        pu *= -4*M_PI/V;
        uu *= -4*M_PI/V;
        
        double dd = - Q1_S1.getZ()*Q1_S2.getZ() 
                    - Q1_S1.getZ()*U1_S2.getZ() 
                    - Q1_S2.getZ()*U1_S1.getZ();
        double qQ =   Q0_S1*Q2_S2.get(2,2) 
                    + Q0_S2*Q2_S1.get(2,2)
                    + Q0_S1*U2_S2.get(2,2)
                    + Q0_S2*U2_S1.get(2,2);
        
        dd *= -4*M_PI/V * EWD::int2eV;
        qQ *= -4*M_PI/V * EWD::int2eV;
        
        LOG(logDEBUG, *log) << (boost::format("  DD %1$+1.7e eV qQ %2$+1.7e eV")
                % dd % qQ) << flush;
                
    }
    else if (shape == "cube" || shape == "sphere") {
        pp = Q0_S1*TrQ2_S2 + Q0_S2*TrQ2_S1 - Q1_S1*Q1_S2;
        pu = Q0_S1*TrU2_S2 + Q0_S2*TrU2_S1 - Q1_S1*U1_S2 - Q1_S2*U1_S1;
        uu = - U1_S1*U1_S2;
        pp *= -4*M_PI/(3*V);
        pu *= -4*M_PI/(3*V);
        uu *= -4*M_PI/(3*V);
        
        if (log != NULL) {
            double pp_qQ = Q0_S1*TrQ2_S2 + Q0_S2*TrQ2_S1;
            double pp_dd =  - Q1_S1*Q1_S2;
            double pu_qQ = Q0_S1*TrU2_S2 + Q0_S2*TrU2_S1;
            double pu_dd =  - Q1_S1*U1_S2 - Q1_S2*U1_S1;
            double uu_dd = - U1_S1*U1_S2;
            pp_qQ *= -4*M_PI/(3*V)*EWD::int2eV;
            pp_dd *= -4*M_PI/(3*V)*EWD::int2eV;
            pu_qQ *= -4*M_PI/(3*V)*EWD::int2eV;
            pu_dd *= -4*M_PI/(3*V)*EWD::int2eV;
            uu_dd *= -4*M_PI/(3*V)*EWD::int2eV;
            LOG(logDEBUG, *log) << (boost::format("Energy-moment decomposition")) << flush;
            LOG(logDEBUG, *log) << (boost::format("  qQ (pp) = %1$+1.7e") % pp_qQ) << flush;
            LOG(logDEBUG, *log) << (boost::format("  dd (pp) = %1$+1.7e") % pp_dd) << flush;
            LOG(logDEBUG, *log) << (boost::format("  qQ (pu) = %1$+1.7e") % pu_qQ) << flush;
            LOG(logDEBUG, *log) << (boost::format("  dd (pu) = %1$+1.7e") % pu_dd) << flush;
            LOG(logDEBUG, *log) << (boost::format("  dd (uu) = %1$+1.7e") % uu_dd) << flush;
        }
    }
    else if (shape == "none") {
    	LOG(logDEBUG, *log) << (boost::format("Assuming isotropic limit")) << flush;
    	pp = 0.0;
    	pu = 0.0;
    	uu = 0.0;
    }
    else {
        cout << endl;
        throw std::runtime_error("Shape '" + shape + "' not implemented");
    }
    
    //cout << endl << "S1 " << Q0_S1 << " " << Q1_S1 << " " << TrQ2_S1 << flush;
    //cout << endl << "S2 " << Q0_S2 << " " << Q1_S2 << " " << TrQ2_S2 << flush;
    
    return EWD::triple<double>(pp, pu, uu);
}

}}
