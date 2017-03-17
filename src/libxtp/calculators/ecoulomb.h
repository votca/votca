/*
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_XTP_ECOULOMB_H
#define __VOTCA_XTP_ECOULOMB_H

#include <votca/xtp/qmcalculator.h>


namespace votca { namespace xtp {
    

class ECoulomb : public QMCalculator
{
public:

    string  Identify() { return "ecoulomb"; }

    void    Initialize(Property *options);
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


void ECoulomb::Initialize(Property *options) {

    // _options already has default values, update them with the supplied options
    _options.CopyValues("", *options );
    string key      = "options." + Identify();
   
    _cutoff = _options.get(key+".cutoff").as< double >();
    _outFile = _options.get(key+".output").as< string >();
    _first = _options.get(key+".first").as< int >();
    _last = _options.get(key+".last").as< int >();
}


bool ECoulomb::EvaluateFrame(Topology *top) {

    // Polar sites generated?
    if (top->isEStatified() == false) {
        cout << endl << "... ... ERROR: Estatify first by running emultipole. " << flush;
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
             vec R_pbc = top->PbShortestConnect((*ext)->getPos(), seg->getPos());
             if (abs(R_pbc) > _cutoff) { continue; }
             vec R_dir = seg->getPos() - (*ext)->getPos();
             vec shift = R_dir - R_pbc;

//             double r12 = abs(top->PbShortestConnect((*ext)->getPos(),
//                                                       seg->getPos()));
//             if (r12 > _cutoff) { continue; }

             // Check polar sites
             // (*ext)->WritePDB(out, "Multipoles", "Charges");

             ++COUNT_EXT;

             for (pit1 = (*ext)->PolarSites().begin();
                  pit1 < (*ext)->PolarSites().end();
                  ++pit1) {

//                  vec dr_pbc = top->PbShortestConnect((*pit1)->getPos(),
//                                                        seg->getPos());
//                  vec dr_dir = seg->getPos() - (*pit1)->getPos();
//                  vec pol_shift = vec(0,0,0);
//                  if (abs(dr_pbc - dr_dir) > 1e-8) {
//                      pol_shift = dr_dir - dr_pbc;
//                  }
                  (*pit1)->PrintPDB(out, com_shift + shift);

             for (pit2 = seg->PolarSites().begin();
                  pit2 < seg->PolarSites().end();
                  ++pit2) {

//                  double R = abs(top->PbShortestConnect((*pit1)->getPos(),
//                                                        (*pit2)->getPos()));

                 double R = abs((*pit1)->getPos() + shift - (*pit2)->getPos());

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
