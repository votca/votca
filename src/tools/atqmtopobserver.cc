#include "atqmtopobserver.h"
#include <votca/csg/nblist.h>
#include <votca/csg/trajectorywriter.h>
#include <qmnblist.h>

AtQmObserver::AtQmObserver()
{}

AtQmObserver::~AtQmObserver()
{}

void AtQmObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);


    _writerCG = TrjWriterFactory().Create(".pdb");
    if(_writerCG == NULL)
        throw runtime_error("output format not supported: .pdb");

    _writerQM = TrjWriterFactory().Create(".pdb");
    if(_writerQM == NULL)
        throw runtime_error("output format not supported: .pdb");

    _writerCG->Open(string("CGtraj.pdb"));
    _writerQM->Open(string("QMtraj.pdb"));
}

void AtQmObserver::EndCG()
{
    _writerCG->Close();
    _writerQM->Close();
}

/// evaluate current conformation
void AtQmObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    _qmtop->Update(*top);
    
    _qmtop->GetJCalc();

    _writerCG->Write(top);
    ///the QM topology is more hard work:
    list<CrgUnit *> lcharges = _qmtop->crglist();
    Topology qmAtomisticTop;
    for (list<CrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++){
        _qmtop->AddAtomisticBeads(*itl,&qmAtomisticTop);
    }
    _writerQM->Write(&qmAtomisticTop);
    qmAtomisticTop.Cleanup();
}


