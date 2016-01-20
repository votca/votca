#ifndef _MOLPOL_H
#define _MOLPOL_H


namespace votca { namespace xtp {

#include <votca/xtp/qmcalculator.h>



class MolPol : public QMCalculator
{
public:

    MolPol() { };
   ~MolPol() { };

    string Identify() { return "molpol"; }

    void                Initialize(Property *options);
    bool                EvaluateFrame(Topology *top);

    matrix              CalculateMolPol(vector<APolarSite*> &poles, bool verb = true);
    int                 SCF_Induce(vector<APolarSite*> &poles);


private:

    string          _mps_input;
    string          _mps_output;
    string          _pol_output;
    
    double          _aDamp;
    double          _wSOR;
    int             _maxIter;
    double          _epsTol;

    BasicInteractor _actor;

};



void MolPol::Initialize(Property *options) {

    // ===================== //
    // LOAD OPTIONS FROM XML //
    // ===================== //

    string key = "options.molpol";

    // MPS INPUT / OUTPUT

    if (options->exists(key+".mps_input")) {
        _mps_input = options->get(key+".mps_input").as<string>();
    }
    else {
        throw std::runtime_error("No .mps input specified in options.");
    }
    if (options->exists(key+".mps_output")) {
        _mps_output = options->get(key+".mps_output").as<string>();
    }
    else {
        throw std::runtime_error("No .mps output specified in options.");
    }
    if (options->exists(key+".pol_output")) {
        _pol_output = options->get(key+".pol_output").as<string>();
    }
    else {
        _pol_output = "";
    }

    // THOLE PARAMETERS

    if ( options->exists(key+".expdamp") ) {
        _aDamp = options->get(key+".expdamp").as< double >();
    }
    else {
        _aDamp = 0.390;
    }

    // CONVERGENCE PARAMETERS

    if ( options->exists(key+".wSOR") ) {
        _wSOR = options->get(key+".wSOR").as< float >();
    }
    else { _wSOR = 0.25; }

    if ( options->exists(key+".maxiter") ) {
        _maxIter = options->get(key+".maxiter").as< int >();
    }
    else { _maxIter = 1024; }

    if ( options->exists(key+".tolerance") ) {
        _epsTol = options->get(key+".tolerance").as< double >();
    }
    else { _epsTol = 0.001; }




    // ===================== //
    // LOAD POLAR SITES      //
    // ===================== //

    vector<APolarSite*> poles = APS_FROM_MPS(_mps_input, 0);



    matrix pmol = this->CalculateMolPol(poles, true);
    pmol *= 1000.0;

    printf("\nPxx Pxy Pxz %+4.7f %+4.7f %+4.7f ", pmol.get(0,0), pmol.get(0,1), pmol.get(0,2));
    printf("\nPyx Pyy Pyz %+4.7f %+4.7f %+4.7f ", pmol.get(1,0), pmol.get(1,1), pmol.get(1,2));
    printf("\nPzx Pzy Pzz %+4.7f %+4.7f %+4.7f ", pmol.get(2,0), pmol.get(2,1), pmol.get(2,2));
    
}


bool MolPol::EvaluateFrame(Topology *top) {

    cout << endl << "... ... Nothing to do here, return." << flush;
    return true;
}


matrix MolPol::CalculateMolPol(vector<APolarSite*> &poles, bool verbose) {

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

        printf("\n\n");
        printf("CONVERGED IN X %4d  Y %4d  Z %4d", iter_x, iter_y, iter_z);
        printf("\n\n");
        
        printf("SUM ATOMIC POLAR XX YY ZZ  (input frame)     "
                "%3.3f %3.3f %3.3f A3\n",
                SUM_alpha_x,
                SUM_alpha_y,
                SUM_alpha_z);
        printf("SUM ATOMIC POLAR 1/3 TRACE                   "
                "%3.3f A3\n",
                SUM_alpha_iso);
        printf("MOLECULAR  POLAR XX YY ZZ  (eigen frame)     "
                "%3.3f %3.3f %3.3f A3\n",
                EIGEN.eigenvalues[0]*NM3_2_A3,
                EIGEN.eigenvalues[1]*NM3_2_A3,
                EIGEN.eigenvalues[2]*NM3_2_A3);
        printf("MOLECULAR  POLAR 1/3 TRACE                   "
                "%3.3f A3", ISO_alpha);
        printf("\n");

        double eigx = EIGEN.eigenvalues[0];
        double eigy = EIGEN.eigenvalues[1];
        double eigz = EIGEN.eigenvalues[2];
        double min_eig =  (eigx < eigy) ?
                         ((eigx < eigz) ? eigx : eigz)
                       : ((eigy < eigz) ? eigy : eigz);
        printf("ANISOTROPY POLAR XX YY ZZ  (eigen frame)     "
                "%2.2f : %2.2f : %2.2f\n",
                eigx/min_eig,
                eigy/min_eig,
                eigz/min_eig);
        printf("POLARIZABLE VOLUME                           "
                "%3.3f A3\n", pow(eigx*eigy*eigz,1./3.)*NM3_2_A3);
        printf("UPPER %4.7f %4.7f %4.7f %4.7f %4.7f %4.7f\n",
               axx*NM3_2_A3, axy*NM3_2_A3, axz*NM3_2_A3,
                             ayy*NM3_2_A3, ayz*NM3_2_A3,
                                           azz*NM3_2_A3);
        printf("\n");


        
    }

    if (verbose && _pol_output != "") {

        FILE *out;
        out = fopen(_pol_output.c_str(), "w");

        vec x = EIGEN.eigenvecs[0];
        vec y = EIGEN.eigenvecs[1];
        vec z = EIGEN.eigenvecs[2];


        fprintf(out, "Polarizability tensor in global frame \n");
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", axx*NM3_2_A3, axy*NM3_2_A3, axz*NM3_2_A3);
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", ayx*NM3_2_A3, ayy*NM3_2_A3, ayz*NM3_2_A3);
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", azx*NM3_2_A3, azy*NM3_2_A3, azz*NM3_2_A3);
        fprintf(out, "\n");
        fprintf(out, "Polarizability tensor in eigenframe \n");
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", EIGEN.eigenvalues[0]*NM3_2_A3,0.,0.);
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", 0.,EIGEN.eigenvalues[1]*NM3_2_A3,0.);
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", 0.,0.,EIGEN.eigenvalues[2]*NM3_2_A3);
        fprintf(out, "\n");
        fprintf(out, "Polarizability tensor main axes \n");
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", x.getX(),x.getY(),x.getZ());
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", y.getX(),y.getY(),y.getZ());
        fprintf(out, "%4.7f  %4.7f  %4.7f \n", z.getX(),z.getY(),z.getZ());

    }

    return alpha;
}


int MolPol::SCF_Induce(vector<APolarSite*> &poles) {

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

#endif