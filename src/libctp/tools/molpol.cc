#include "molpol.h"
#include <votca/ctp/polartop.h>

namespace votca { namespace ctp {
    
void MolPolTool::Initialize(Property *options) {

    string key = "options.molpol.";

    // MPS INPUT / OUTPUT    
    if (options->exists(key+"mpsfiles.input")) {
        _mps_input = options->get(key+"mpsfiles.input").as<string>();
    }
    else {
        throw std::runtime_error("No .mps input specified in options.");
    }
    if (options->exists(key+"mpsfiles.output")) {
        _mps_output = options->get(key+"mpsfiles.output").as<string>();
    }
    else {
        throw std::runtime_error("No .mps output specified in options.");
    }
    if (options->exists(key+"mpsfiles.polar")) {
        _pol_output = options->get(key+"mpsfiles.polar").as<string>();
    }
    else {
        _pol_output = "molpol." + _mps_output + ".xml";
    }

    // THOLE & INDUCTION CONVERGENCE PARAMETERS
    if ( options->exists(key+"induction.expdamp") ) {
        _aDamp = options->get(key+"induction.expdamp").as< double >();
    }
    else {
        _aDamp = 0.390;
    }
    if ( options->exists(key+"induction.wSOR") ) {
        _wSOR = options->get(key+"induction.wSOR").as< float >();
    }
    else { _wSOR = 0.25; }

    if ( options->exists(key+"induction.maxiter") ) {
        _maxIter = options->get(key+"induction.maxiter").as< int >();
    }
    else { _maxIter = 1024; }

    if ( options->exists(key+"induction.tolerance") ) {
        _epsTol = options->get(key+"induction.tolerance").as< double >();
    }
    else { _epsTol = 0.001; }
    
    // TARGET MOLECULAR POLARIZABILITY
    if (options->exists(key+"target.molpol")) {
        vector<double> u = options->get(key+"target.molpol").as< vector<double> >();
        if (u.size() != 6) {
            cout << endl << "... ... ERROR <options.molpol.target_molpol> "
                " should have this format: pxx pxy pxz pyy pyz pzz" << endl;
            throw std::runtime_error("(input error, see above)");
        }
        double A3_to_nm3 = 0.001;        
        _target = A3_to_nm3 * matrix(vec(u[0],u[1],u[2]),
                                     vec(u[1],u[3],u[4]),
                                     vec(u[2],u[4],u[5]));
    }
    if (options->exists(key+"target.optimize")) {
        _do_optimize = options->get(key+"target.optimize").as<bool>();
    }
    else {
        _do_optimize = false;
    }
    if (options->exists(key+"target.tolerance")) {
        _tolerance = options->get(key+"target.tolerance").as<double>();
    }
    else {
        _tolerance = 1e-3;
    }
}


bool MolPolTool::Evaluate() {
    
    // Load polar sites
    cout << endl << "... ... Start molpol - load polar sites" << flush;
    vector<APolarSite*> poles = APS_FROM_MPS(_mps_input, 0);
    PolarSeg pseg_input = PolarSeg(1, poles);
    PolarSeg *pseg_output = 0;
    
    if (!_do_optimize) {
        // Output only
        pseg_input.WriteMPS(_mps_output, "MOLPOL (UNSCALED)");
        matrix pmol = CalculateMolPol(pseg_input, true);
    }
    
    else {        
        
        double scale = 1.0;
        double omega = 1.0;
        
        bool converged = false;
        int loop_count = 0;
        int max_iter = 1024;
        while (true) {
            
            PolarSeg pseg_inter = PolarSeg(&pseg_input);
            
            if (tools::globals::verbose)
                cout << endl << "Relative scale s = " << scale << flush;
            for (PolarSeg::iterator pit = pseg_inter.begin();
                    pit < pseg_inter.end(); ++pit) {
                (*pit)->Pxx = omega*(*pit)->Pxx*scale + (1-omega)*(*pit)->Pxx;
                (*pit)->Pxy = omega*(*pit)->Pxy*scale + (1-omega)*(*pit)->Pxy;
                (*pit)->Pxz = omega*(*pit)->Pxz*scale + (1-omega)*(*pit)->Pxz;
                (*pit)->Pyy = omega*(*pit)->Pyy*scale + (1-omega)*(*pit)->Pyy;
                (*pit)->Pyz = omega*(*pit)->Pyz*scale + (1-omega)*(*pit)->Pyz;
                (*pit)->Pzz = omega*(*pit)->Pzz*scale + (1-omega)*(*pit)->Pzz;
            }
            
            // Polarizability tensor "inter"
            matrix p_inter = this->CalculateMolPol(pseg_inter, false);
            // ... Eigen frame
            matrix::eigensystem_t eig_inter;
            p_inter.SolveEigensystem(eig_inter);
            double eig_inter_x = eig_inter.eigenvalues[0];
            double eig_inter_y = eig_inter.eigenvalues[1];
            double eig_inter_z = eig_inter.eigenvalues[2];
            // ... Polarizable volume
            double pvol_inter = 1000*
                pow(eig_inter_x*eig_inter_y*eig_inter_z,1./3.);
            
            // Polarizability tensor "target"
            matrix p_target = _target;
            // ... Eigen frame
            matrix::eigensystem_t eig_target;
            p_target.SolveEigensystem(eig_target);
            double eig_target_x = eig_target.eigenvalues[0];
            double eig_target_y = eig_target.eigenvalues[1];
            double eig_target_z = eig_target.eigenvalues[2];
            // ... Polarizable volume
            double pvol_target = 1000*
                pow(eig_target_x*eig_target_y*eig_target_z,1./3.);
            
            // Adjust scaling parameter
            scale *= pow(pvol_target/pvol_inter,2);
            
            // Convergence test
            loop_count += 1;
            if (std::abs(1-pvol_target/pvol_inter) < _tolerance) {
                cout << endl << "... ... Iterative refinement : *CONVERGED*" << flush;
                converged = true;
                pseg_output = new PolarSeg(&pseg_inter);
                break;
            }
            if (loop_count > max_iter) {
                cout << endl << "... ... Iterative refinement : *FAILED*";
                break;
            }
        }
        
        // Output
        pseg_output->WriteMPS(_mps_output, "MOLPOL (OPTIMIZED)");
        matrix pmol = CalculateMolPol(*pseg_output, true);
    }
    
    delete pseg_output;
    pseg_output = 0;
    
    return true;
}


matrix MolPolTool::CalculateMolPol(vector<APolarSite*> &poles, bool verbose) {

    vector<APolarSite*> ::iterator pit;
    //for (pit = poles.begin(); pit < poles.end(); ++pit) {
    //    (*pit)->PrintInfo(cout);
    //}

    // Do not forget to define exp. damping factor
    _actor.SetADamp(_aDamp);


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
        (*pit)->FPx = extFx;
        (*pit)->FPy = extFy;
        (*pit)->FPz = extFz;
    }

    // Calculate induction field
    int iter_x = this->SCF_Induce(poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        siteU1x += (*pit)->U1x;
        siteU1y += (*pit)->U1y;
        siteU1z += (*pit)->U1z;
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
        (*pit)->FPx = extFx;
        (*pit)->FPy = extFy;
        (*pit)->FPz = extFz;
    }

    // Calculate induction field
    int iter_y = this->SCF_Induce(poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        siteU1x += (*pit)->U1x;
        siteU1y += (*pit)->U1y;
        siteU1z += (*pit)->U1z;
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
        (*pit)->FPx = extFx;
        (*pit)->FPy = extFy;
        (*pit)->FPz = extFz;
    }

    // Calculate induction field
    int iter_z = this->SCF_Induce(poles);

    // Add up ind. dpl.s to yield molecular U1
    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        siteU1x += (*pit)->U1x;
        siteU1y += (*pit)->U1y;
        siteU1z += (*pit)->U1z;
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

    for (pit = poles.begin(); pit < poles.end(); ++pit) {
        SUM_alpha_x += (*pit)->Pxx;
        SUM_alpha_y += (*pit)->Pyy;
        SUM_alpha_z += (*pit)->Pzz;
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






    
    if (verbose) {
        
        if (tools::globals::verbose) {
            printf("\n\n");
            printf("CONVERGED IN X %4d  Y %4d  Z %4d", iter_x, iter_y, iter_z);
            printf("\n\n");
        }
        
        printf("\n\nSummary mps='%s'\n", _mps_output.c_str());
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

    if (verbose && _pol_output != "") {

        FILE *out;
        out = fopen(_pol_output.c_str(), "w");

        vec x = EIGEN.eigenvecs[0];
        vec y = EIGEN.eigenvecs[1];
        vec z = EIGEN.eigenvecs[2];

        fprintf(out, "<polarizability mps='%s'>\n", _mps_output.c_str());
        fprintf(out, "\t<polartensor type='scf' frame='global'>\n");
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e\n", axx*NM3_2_A3, axy*NM3_2_A3, axz*NM3_2_A3);
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e\n", ayx*NM3_2_A3, ayy*NM3_2_A3, ayz*NM3_2_A3);
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e\n", azx*NM3_2_A3, azy*NM3_2_A3, azz*NM3_2_A3);
        fprintf(out, "\t</polartensor>\n");
        fprintf(out, "\t<polartensor type='scf' frame='eigen'>\n");
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", EIGEN.eigenvalues[0]*NM3_2_A3,0.,0.);
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", 0.,EIGEN.eigenvalues[1]*NM3_2_A3,0.);
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", 0.,0.,EIGEN.eigenvalues[2]*NM3_2_A3);
        fprintf(out, "\t</polartensor>\n");
        fprintf(out, "\t<principalaxes>\n");
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", x.getX(),x.getY(),x.getZ());
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", y.getX(),y.getY(),y.getZ());
        fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", z.getX(),z.getY(),z.getZ());
        fprintf(out, "\t</principalaxes>\n");
        fprintf(out, "</polarizability>");
    }

    return alpha;
}


int MolPolTool::SCF_Induce(vector<APolarSite*> &poles) {

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