#include "molpol.h"
#include <votca/xtp/polartop.h>

namespace votca { namespace xtp {
    
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
        
    // MOLPOL ENGINE SET-UP
    _molpolengine = MolPolEngine(_aDamp, _wSOR, _maxIter, _epsTol);
    
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
    if (options->exists(key+"target.pattern")) {
        string tmp = options->get(key+"target.pattern").as<string>();
        Tokenizer toker(tmp, " ,\t\n");
        toker.ToVector(_scaling_pattern);
        cout << endl << "... ... Using scaling pattern (length " 
            << _scaling_pattern.size() << ")" << flush;
    }
    else {
        ;
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
        //matrix pmol = _molpolengine.CalculateMolPol(pseg_input, true);
    }
    
    else {
        // Setup scaling pattern (Y = full scale, H = half scale, N = zero scale)
        if (_scaling_pattern.size() == 0) {
            for (PolarSeg::iterator pit = pseg_input.begin();
                pit < pseg_input.end(); ++pit) {
                _scaling_pattern.push_back("Y");
            }
        }
        else {
            assert(_scaling_pattern.size() == pseg_input.size() &&
                "<molpol> Scaling pattern from input does not match size of"
                "polar segment");
        }
        
        // Refine until converged
        double scale = 1.0;
        double omega = 1.0;
        
        bool converged = false;
        int loop_count = 0;
        int max_iter = 1024;
        while (true) {
            
            PolarSeg pseg_inter = PolarSeg(&pseg_input, true);
            
            if (tools::globals::verbose)
                cout << endl << "Relative scale s = " << scale << flush;
            vector<string>::iterator strit = _scaling_pattern.begin();
            for (PolarSeg::iterator pit = pseg_inter.begin();
                    pit < pseg_inter.end(); ++pit, ++strit) {
                double local_scale = 0.;
                if ((*strit) == "Y")
                    local_scale = scale;
                else if ((*strit) == "H") {
                    local_scale = 1+(scale-1.)*0.5;
                }
                else if ((*strit) == "N")
                    local_scale = 1.;
                else
                    assert(false && "<molpol> Unsupported scaling instruction");
                (*pit)->Pxx = omega*(*pit)->Pxx*local_scale + (1-omega)*(*pit)->Pxx;
                (*pit)->Pxy = omega*(*pit)->Pxy*local_scale + (1-omega)*(*pit)->Pxy;
                (*pit)->Pxz = omega*(*pit)->Pxz*local_scale + (1-omega)*(*pit)->Pxz;
                (*pit)->Pyy = omega*(*pit)->Pyy*local_scale + (1-omega)*(*pit)->Pyy;
                (*pit)->Pyz = omega*(*pit)->Pyz*local_scale + (1-omega)*(*pit)->Pyz;
                (*pit)->Pzz = omega*(*pit)->Pzz*local_scale + (1-omega)*(*pit)->Pzz;
            }
            
            // Polarizability tensor "inter"
            matrix p_inter = _molpolengine.CalculateMolPol(
                pseg_inter, tools::globals::verbose);
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
            if (tools::globals::verbose) 
                cout << endl << "Target volume " << pvol_target 
                    << ", current volume " << pvol_inter << flush;
            scale *= pow(pvol_target/pvol_inter,2);
            
            // Convergence test
            loop_count += 1;
            if (std::abs(1-pvol_target/pvol_inter) < _tolerance) {
                cout << endl << "... ... Iterative refinement : *CONVERGED*" << flush;
                cout << endl << "... ... Scaling coefficient  : " << scale << flush;
                converged = true;
                pseg_output = new PolarSeg(&pseg_inter, true);
                break;
            }
            if (loop_count > max_iter) {
                cout << endl << "... ... Iterative refinement : *FAILED*" << flush;
                converged = false;
                pseg_output = new PolarSeg(&pseg_inter, true);
                break;
            }
        }
        
        // Output
        if (converged)
            pseg_output->WriteMPS(_mps_output, "MOLPOL (OPTIMIZED)");
        else {
            cout << endl << "... ... ERROR Convergence not achieved. "
                 << "Check your input mps-file, target polarizability <target> "
                 << "or try decreasing <wSOR>." 
                 << flush;
        }        
        
        if (converged && _pol_output != "") {
            matrix pmol = _molpolengine.CalculateMolPol(*pseg_output, true);
            
            FILE *out;
            out = fopen(_pol_output.c_str(), "w");            
            
            double NM3_2_A3 = 1000.;
            double axx = pmol.get(0,0); double axy = pmol.get(0,1); double axz = pmol.get(0,2);
            double ayx = pmol.get(1,0); double ayy = pmol.get(1,1); double ayz = pmol.get(1,2);
            double azx = pmol.get(2,0); double azy = pmol.get(2,1); double azz = pmol.get(2,2);

            matrix::eigensystem_t pmol_eigen;
            pmol.SolveEigensystem(pmol_eigen);
            
            vec x = pmol_eigen.eigenvecs[0];
            vec y = pmol_eigen.eigenvecs[1];
            vec z = pmol_eigen.eigenvecs[2];

            fprintf(out, "<polarizability mps='%s'>\n", _mps_output.c_str());
            fprintf(out, "\t<polartensor type='scf' frame='global'>\n");
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e\n", axx*NM3_2_A3, axy*NM3_2_A3, axz*NM3_2_A3);
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e\n", ayx*NM3_2_A3, ayy*NM3_2_A3, ayz*NM3_2_A3);
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e\n", azx*NM3_2_A3, azy*NM3_2_A3, azz*NM3_2_A3);
            fprintf(out, "\t</polartensor>\n");
            fprintf(out, "\t<polartensor type='scf' frame='eigen'>\n");
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", pmol_eigen.eigenvalues[0]*NM3_2_A3,0.,0.);
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", 0.,pmol_eigen.eigenvalues[1]*NM3_2_A3,0.);
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", 0.,0.,pmol_eigen.eigenvalues[2]*NM3_2_A3);
            fprintf(out, "\t</polartensor>\n");
            fprintf(out, "\t<principalaxes>\n");
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", x.getX(),x.getY(),x.getZ());
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", y.getX(),y.getY(),y.getZ());
            fprintf(out, "\t\t%+1.7e  %+1.7e  %+1.7e \n", z.getX(),z.getY(),z.getZ());
            fprintf(out, "\t</principalaxes>\n");
            fprintf(out, "</polarizability>");
        }
    }
    
    if (tools::globals::verbose) {
        pseg_input.Coarsegrain(true);
        pseg_input.WriteMPS("molpol.coarse."+_mps_output,
            "MOLPOL (UNSCALED, COARSE)");
    }
    
    delete pseg_output;
    pseg_output = 0;
    
    return true;
}

}}
