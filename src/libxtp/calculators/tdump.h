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


#ifndef _TDUMP2_H
#define _TDUMP2_H

#include <cstdlib>
#include <votca/xtp/qmcalculator.h>


namespace votca { namespace xtp {


class TDump : public QMCalculator
{

public:

    TDump() : _outPDBmd("MD.pdb"), _outPDBqm("QM.pdb"), 
              _framesToWrite(1),   _framesWritten(0)     { };
   ~TDump() { };

    string Identify() { return "tdump"; }

    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

    string _outPDBmd;
    string _outPDBqm;

    int    _framesToWrite;
    int    _framesWritten;

    
};


void TDump::Initialize(Property *options) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( options );
 
    string key      = "options." + Identify();
    _outPDBmd = options->get(key+".md").as<string>();
    _outPDBqm = options->get(key+".qm").as<string>();
    _framesToWrite = options->get(key+".frames").as<int>();

}

bool TDump::EvaluateFrame(Topology *top) {

    if (_framesWritten > _framesToWrite) { return 1; }

    // Rigidify std::system (if possible)
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
    
    return true;
}

}}


#endif
