#ifndef ECOULOMB_H
#define ECOULOMB_H

#include <votca/ctp/qmcalculator.h>


namespace votca { namespace ctp {
    

class ECoulomb : public QMCalculator
{
public:

    string  Identify() { return "ECoulomb"; }

    void    Initialize(Topology *top, Property *options);
    bool    EvaluateFrame(Topology *top);
    void    EvaluateSegment(Topology *top, Segment *seg, int state);
    void    Output2File(Topology *top);

private:

    double          _cutoff;
    string          _outFile;
    map< int, int > _log_seg_sphere; // <- # segs within cut-off

    int             _first;
    int             _last;
    
};


void ECoulomb::Initialize(Topology *top, Property *options) {

    string key = "options.ecoulomb";

    cout << endl << "... ... Init" << flush;    
    
    _cutoff = options->get(key+".cutoff").as< double >();
    _outFile = options->get(key+".output").as< string >();
    if (options->exists(key+".first")) {
        _first = options->get(key+".first").as< int >();
    }
    else {
        _first = 1;
    }
    if (options->exists(key+".last")) {
        _last = options->get(key+".last").as< int >();
    }
    else {
        _last = -1;
    }
}


bool ECoulomb::EvaluateFrame(Topology *top) {

    // Polar sites generated?
    if (top->isEStatified() == false) {
        cout << endl << "... ... ERROR: Estatify first. Return. " << flush;
        return 0;
    }
    else { 
        cout << endl << "... ... System is estatified. Proceed. " << flush;
    }

    vector< Segment* > ::iterator sit;
    vector< PolarSite* > ::iterator pit;

    // Reset all polar sites to neutral state
    for (sit = top->Segments().begin(); 
         sit < top->Segments().end();
         ++sit) {
    for (pit = (*sit)->PolarSites().begin();
         pit < (*sit)->PolarSites().end();
         ++pit) {
        (*pit)->Charge(0); 
    }}

    // Calculate energies for segments
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

         if ((*sit)->getId() < _first) { continue; }
         if ((*sit)->getId() == _last+1)  { break; }

         cout << endl << "... ... Evaluating site " << (*sit)->getId() << flush;

    for (int state = -1; state < 2; ++state) {
         
         this->EvaluateSegment(top, *sit, state);                 
    }}

    // Output results
    this->Output2File(top);

    return true;
}


void ECoulomb::EvaluateSegment(Topology *top, Segment *seg, int state) {    

    if (seg->hasChrgState(state)) {        

        vector< Segment*   > ::iterator ext;
        vector< PolarSite* > ::iterator pit1;
        vector< PolarSite* > ::iterator pit2;

        // ++++++++++++++++++++++++++++ //
        // Charge segment appropriately //
        // ++++++++++++++++++++++++++++ //

        for (pit1 = seg->PolarSites().begin();
             pit1 < seg->PolarSites().end();
             ++pit1) {

             (*pit1)->Charge(state);
        }

        double E_INTER   = 0.0;
        int    COUNT_EXT = 0;
        double INT2eV    = 1/(4*M_PI*8.854187817e-12)
                         * 1.602176487e-19 / 1.000e-9;

        // ++++++++++++++++++++++++++++ //
        // Calculate interaction energy //
        // ++++++++++++++++++++++++++++ //

        FILE *out;
        string sphpdb = boost::lexical_cast<string>(seg->getId())
                      + "_" + boost::lexical_cast<string>(state)
                      + "_check.pdb";
        out = fopen(sphpdb.c_str(), "w");
        vec com_shift = -1 * seg->getPos();
        for (pit1 = seg->PolarSites().begin();
             pit1 < seg->PolarSites().end();
             ++pit1) {
             (*pit1)->PrintPDB(out, com_shift);
        }

        for (ext = top->Segments().begin();
             ext < top->Segments().end();
             ++ext) {

             // Segment different from central one?
             if ((*ext)->getId() == seg->getId()) { continue; }             
             
             // Segment within cutoff?
             double r12 = abs(top->PbShortestConnect((*ext)->getPos(),
                                                       seg->getPos()));
             if (r12 > _cutoff) { continue; }

             // Check polar sites
             // (*ext)->WritePDB(out, "Multipoles", "Charges");

             ++COUNT_EXT;

             for (pit1 = (*ext)->PolarSites().begin();
                  pit1 < (*ext)->PolarSites().end();
                  ++pit1) {

                  vec dr_pbc = top->PbShortestConnect((*pit1)->getPos(),
                                                        seg->getPos());
                  vec dr_dir = seg->getPos() - (*pit1)->getPos();
                  vec pol_shift = vec(0,0,0);
                  if (abs(dr_pbc - dr_dir) > 1e-8) {
                      pol_shift = dr_dir - dr_pbc;
                  }
                  (*pit1)->PrintPDB(out, com_shift + pol_shift);

             for (pit2 = seg->PolarSites().begin();
                  pit2 < seg->PolarSites().end();
                  ++pit2) {

                  double R = abs(top->PbShortestConnect((*pit1)->getPos(),
                                                        (*pit2)->getPos()));

                  E_INTER += (*pit1)->Q00 * 1/R * (*pit2)->Q00;
             }}
        }

        fclose(out);

        // ++++++++++++++++++++++++++++++ //
        // Store result, reset to neutral //
        // ++++++++++++++++++++++++++++++ //

        seg->setEMpoles(state, INT2eV * E_INTER);
        _log_seg_sphere[seg->getId()] = COUNT_EXT;
        cout << endl << "... ... ... STATE " << state << ": "
             << " EP(~) = " << E_INTER*INT2eV << " eV"
             << " SPH(o) = " << COUNT_EXT+1 << flush;

        for (pit1 = seg->PolarSites().begin();
             pit1 < seg->PolarSites().end();
             ++pit1) {
             (*pit1)->Charge(0);
        }
    }
}


void ECoulomb::Output2File(Topology *top) {

    FILE *out;
    string topOutFile = "Frame"
                      + boost::lexical_cast<string>(top->getDatabaseId())
                      + "_" + _outFile;
    out = fopen(topOutFile.c_str(), "w");

    vector< Segment* > ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         ++sit) {

         if ((*sit)->getId() < _first) { continue; }
         else if ((*sit)->getId() == _last+1)  { break; }

        fprintf(out, "%4d ", (*sit)->getId() );
        fprintf(out, "%4s ", (*sit)->getName().c_str() );

        // Energies
        if ((*sit)->hasChrgState(0)) {
            fprintf(out, "   0 %3.8f   ", (*sit)->getEMpoles(0) );
        }
        if ((*sit)->hasChrgState(-1)) {
            fprintf(out, "  -1 %3.8f   ", (*sit)->getEMpoles(-1) );
        }

        if ((*sit)->hasChrgState(+1)) {
            fprintf(out, "  +1 %3.8f   ", (*sit)->getEMpoles(+1) );
        }

        // Iterations
        if ((*sit)->hasChrgState(0)) {
            fprintf(out, "   0 0     ");
        }
        if ((*sit)->hasChrgState(-1)) {
            fprintf(out, "  -1 0     ");
        }

        if ((*sit)->hasChrgState(+1)) {
            fprintf(out, "  +1 0     ");
        }

        // Polarizable sphere
        fprintf(out, "   EXT %4d   ",
                _log_seg_sphere[(*sit)->getId()]);

        fprintf(out, "   %4.7f %4.7f %4.7f   ",
                (*sit)->getPos().getX(),
                (*sit)->getPos().getY(),
                (*sit)->getPos().getZ());

        fprintf(out, " \n");
    }
}






}}

#endif
