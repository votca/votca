#ifndef _TDUMP2_H
#define _TDUMP2_H

#include <cstdlib>
#include <votca/ctp/qmcalculator2.h>


namespace votca { namespace ctp {




class TDump : public QMCalculator2 {

public:

    TDump() : _outPDBmd("MD.pdb"), _outPDBqm("QM.pdb"), 
              _framesToWrite(1),   _framesWritten(0)     { };
   ~TDump() { };

    string Identify() { return "Trajectory dump"; }

    void Initialize(Topology *top, Property *options);
    bool EvaluateFrame(Topology *top);

private:

    string _outPDBmd;
    string _outPDBqm;

    int    _framesWritten;
    int    _framesToWrite;

    
};


void TDump::Initialize(Topology *top, Property *options) {

    string key = "options.tdump";

    if (options->exists(key+".md") && options->exists(key+".qm")) {
        _outPDBmd = options->get(key+".md").as<string>();
        _outPDBqm = options->get(key+".qm").as<string>();
    }

    if (options->exists(key+".frames")) {
        _framesToWrite = options->get(key+".frames").as<int>();
    }

}

bool TDump::EvaluateFrame(Topology *top) {

    if (_framesWritten > _framesToWrite) { return 1; }

    // Rigidify system (if possible)
    bool isRigid = top->Rigidify();
    if (!isRigid) {
        return 0;
    }

    // Print coordinates
    FILE *outPDBmd = NULL;
    FILE *outPDBqm = NULL;

    outPDBmd = fopen(_outPDBmd.c_str(), "a");
    outPDBqm = fopen(_outPDBqm.c_str(), "a");

    fprintf(outPDBmd, "TITLE     VOT CAtastrophic title \n");
    fprintf(outPDBmd, "Model %8d \n" , _framesWritten);

    fprintf(outPDBqm, "TITLE     VOT CAtastrophic title \n");
    fprintf(outPDBqm, "Model %8d \n", _framesWritten);    

    vector< Segment* > ::iterator sit;
    for (sit = top->Segments().begin();
         sit < top->Segments().end();
         sit++) {

        (*sit)->WritePDB(outPDBmd, "Atoms", "MD");
        (*sit)->WritePDB(outPDBqm, "Atoms", "QM");
    }

    fprintf(outPDBmd, "TER\nENDMDL\n");
    fprintf(outPDBqm, "TER\nENDMDL\n");

    fclose(outPDBmd);
    fclose(outPDBqm);

    _framesWritten++;
}

}}


#endif
