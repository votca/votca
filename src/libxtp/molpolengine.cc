#include <votca/xtp/molpolengine.h>

namespace votca {
namespace xtp {

    
matrix MolPolEngine::CalculateMolPol(vector<APolarSite*> &poles, bool verbose) {

    vector<APolarSite*> ::iterator pit;

    double axx, axy, axz;             // |
    double ayx, ayy, ayz;             // |-> Polarizability tensor
    double azx, azy, azz;             // |

    double siteU1x, siteU1y, siteU1z; // |-> Molecular ind. dipole

    double extFx, extFy, extFz;       // |-> External applied field

    // +++++++++++++++++++++++++++++ //
    // External field in x-direction //
    // +++++++++++++++++++++++++++++ //

    siteU1x = siteU1y = siteU1z = 0.0;

    extFx = 0.1;
    extFy = 0.0;
    extFz = 0.0;

    // Set permanent field
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        (*pit)->setFieldP(extFx, extFy, extFz);
    }

    // Calculate induction field
    int iter_x = this->SCF_Induce(poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        vec u1 = (*pit)->getU1();
        siteU1x += u1.getX();
        siteU1y += u1.getY();
        siteU1z += u1.getZ();
    }

    // Calculate associated column of polarizability tensor
    axx = - siteU1x / extFx;
    ayx = - siteU1y / extFx;
    azx = - siteU1z / extFx;

    // +++++++++++++++++++++++++++++ //
    // External field in y-direction //
    // +++++++++++++++++++++++++++++ //

    siteU1x = siteU1y = siteU1z = 0.0;

    extFx = 0.0;
    extFy = 0.1;
    extFz = 0.0;

    // Set permanent field
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        (*pit)->setFieldP(extFx, extFy, extFz);
    }

    // Calculate induction field
    int iter_y = this->SCF_Induce(poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        vec u1 = (*pit)->getU1();
        siteU1x += u1.getX();
        siteU1y += u1.getY();
        siteU1z += u1.getZ();
    }

    // Calculate associated column of polarizability tensor
    axy = - siteU1x / extFy;
    ayy = - siteU1y / extFy;
    azy = - siteU1z / extFy;

    // +++++++++++++++++++++++++++++ //
    // External field in z-direction //
    // +++++++++++++++++++++++++++++ //

    siteU1x = siteU1y = siteU1z = 0.0;

    extFx = 0.0;
    extFy = 0.0;
    extFz = 0.1;

    // Set permanent field
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        (*pit)->setFieldP(extFx, extFy, extFz);
    }

    // Calculate induction field
    int iter_z = this->SCF_Induce(poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        vec u1 = (*pit)->getU1();
        siteU1x += u1.getX();
        siteU1y += u1.getY();
        siteU1z += u1.getZ();
    }

    // Calculate associated column of polarizability tensor
    axz = - siteU1x / extFz;
    ayz = - siteU1y / extFz;
    azz = - siteU1z / extFz;

    // +++++++++++++++++++++++++++ //
    // Sum, Trace, Diagonalization //
    // +++++++++++++++++++++++++++ //

    // Sum over atomic polarizabilities
    double SUM_alpha_x = 0.0;
    double SUM_alpha_y = 0.0;
    double SUM_alpha_z = 0.0;
    double SUM_alpha_iso = 0.0;

    vec xdir = vec(1,0,0);
    vec ydir = vec(0,1,0);
    vec zdir = vec(0,0,1);
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        SUM_alpha_x += (*pit)->getProjP(xdir);
        SUM_alpha_y += (*pit)->getProjP(ydir);
        SUM_alpha_z += (*pit)->getProjP(zdir);
        SUM_alpha_iso += (*pit)->getIsoP();
    }

    double NM3_2_A3 = 1000.;
    SUM_alpha_x *= NM3_2_A3;
    SUM_alpha_y *= NM3_2_A3;
    SUM_alpha_z *= NM3_2_A3;
    SUM_alpha_iso *= NM3_2_A3;


    double ISO_alpha = (axx + ayy + azz) / 3.;
    ISO_alpha *= NM3_2_A3;

    // Eigenvalues of polarizability tensor
    matrix alpha = matrix( vec(axx,ayx,azx),
                           vec(axy,ayy,azy),
                           vec(axz,ayz,azz) );
    matrix::eigensystem_t EIGEN;
    alpha.SolveEigensystem(EIGEN);
    
    // Print info if verbose
    if (verbose) {
        
        if (tools::globals::verbose) {
            printf("\n\n");
            printf("CONVERGED IN X %4d  Y %4d  Z %4d", iter_x, iter_y, iter_z);
            printf("\n\n");
        }
        
        printf("\n\nSummary\n");
        printf("  o atomic polarizabilities, xx yy zz (sum, input frame)    "
                "%+1.3e %+1.3e %+1.3e A**3\n",
                SUM_alpha_x,
                SUM_alpha_y,
                SUM_alpha_z);
        printf("  o molecular polarizability xx yy zz (scf, eigen frame)    "
                "%+1.3e %+1.3e %+1.3e A**3\n",
                EIGEN.eigenvalues[0]*NM3_2_A3,
                EIGEN.eigenvalues[1]*NM3_2_A3,
                EIGEN.eigenvalues[2]*NM3_2_A3);

        double eigx = EIGEN.eigenvalues[0];
        double eigy = EIGEN.eigenvalues[1];
        double eigz = EIGEN.eigenvalues[2];
        double min_eig =  (eigx < eigy) ?
                         ((eigx < eigz) ? eigx : eigz)
                       : ((eigy < eigz) ? eigy : eigz);
        printf("  o molecular anisotropy xx:yy:zz     (scf, eigen frame)    "
                "%2.2f : %2.2f : %2.2f\n",
                eigx/min_eig,
                eigy/min_eig,
                eigz/min_eig);
        printf("  o polarizable volume (w/o 4PI/3)    (scf, eigen frame)    "
                "%3.3f A3\n", pow(eigx*eigy*eigz,1./3.)*NM3_2_A3);
        printf("  o upper polarizability tensor       (scf, input frame)    "
                "%+1.3e %+1.3e %+1.3e %+1.3e %+1.3e %+1.3e\n",
               axx*NM3_2_A3, axy*NM3_2_A3, axz*NM3_2_A3,
                             ayy*NM3_2_A3, ayz*NM3_2_A3,
                                           azz*NM3_2_A3);
    }

    return alpha;
}


int MolPolEngine::SCF_Induce(vector<APolarSite*> &poles) {

    int    maxIter  = this->_maxIter;
    double wSOR     = this->_wSOR;
    double eTOL     = this->_epsTol;

    vector< APolarSite* >::iterator pit1;
    vector< APolarSite* >::iterator pit2;

    // Induce to first order
    for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
        (*pit1)->InduceDirect();
    }

    int iter = 0;
    for ( ; iter < maxIter; ++iter) {

        // Reset induction field
        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
            (*pit1)->ResetFieldU();
        }

        // Calculate higher-order induction field
        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
            for (pit2 = pit1 + 1; pit2 < poles.end(); ++pit2) {
                _actor.FieldInduAlpha(*(*pit1),*(*pit2));
            }
        }

        // Induce dipoles
        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
            (*pit1)->Induce(wSOR);
        }

        // Check for convergence
        bool converged = true;
        double maxdU = -1;
        double avgdU = 0.0;
        int    baseN = 0;

        for (pit1 = poles.begin(); pit1 < poles.end(); ++pit1) {
             double dU = (*pit1)->HistdU();
             avgdU += dU;
             ++baseN;
             if ( dU > maxdU ) { maxdU = dU; }
             if ( dU > eTOL ) { converged = false; }
        }

        avgdU /= baseN;
        if (avgdU < eTOL/10.) { converged = true; }

        //cout << endl << "ITER" << iter << " | MAX dU " << maxdU
        //     << " | AVG dU " << avgdU;

        // Break if converged
        if (converged) { break; }
        else if (iter == maxIter - 1) {
            cout << endl << "... ... ... "
                 << "WARNING Induced multipoles did not converge to precision: "
                 << " AVG dU:U " << avgdU << flush;
            break;
        }
    }

    return iter;
}
    
    
    
    
    
    
}}