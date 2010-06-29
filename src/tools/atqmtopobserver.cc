#include "atqmtopobserver.h"
#include <votca/csg/nblist.h>
#include <votca/csg/trajectorywriter.h>
#include <qmnblist.h>

AtQmObserver::AtQmObserver()
{}

AtQmObserver::~AtQmObserver()
{}

void AtQmObserver::HelpText(){
    cout << "no help text available\n\n";
    cout << _op_desc << endl;
}


void AtQmObserver::AddSpecificOptions(){
    _op_desc.add_options()
    ("outCG", boost::program_options::value<string>()->default_value("CGtraj.pdb"), "  the output file for coarse grained topology")
    ("outQM", boost::program_options::value<string>()->default_value("QMtraj.pdb"), " the output file for the QM geometry")
    ;
}

void AtQmObserver::Initialize()
{
    TrajectoryWriter::RegisterPlugins();
    
    string nameCG = _op_vm["outCG"].as<string>();
    string nameQM = _op_vm["outQM"].as<string>();
    string extCG  = nameCG.substr(nameCG.length()-4,4);
    string extQM  = nameCG.substr(nameQM.length()-4,4);
    _writerCG = TrjWriterFactory().Create(nameCG);
    if(_writerCG == NULL)
        throw runtime_error(string("output format not supported: ")+ extCG);

    _writerQM = TrjWriterFactory().Create(nameQM);
    if(_writerQM == NULL)
        throw runtime_error(string("output format not supported: ")+ extQM);

    _writerCG->Open(nameCG);
    _writerQM->Open(nameQM);
}

void AtQmObserver::EndEvaluate()
{
    _writerCG->Close();
    _writerQM->Close();
}

/// evaluate current conformation
bool AtQmObserver::EvaluateFrame()
{
    _writerCG->Write(&_qmtop);
    ///the QM topology is more hard work:
    vector<CrgUnit *> lcharges = _qmtop.CrgUnits();
    Topology qmAtomisticTop;
    for (vector<CrgUnit *>::iterator itl = lcharges.begin(); itl != lcharges.end(); itl++){
        _qmtop.AddAtomisticBeads(*itl,&qmAtomisticTop);
    }
    _writerQM->Write(&qmAtomisticTop);
    qmAtomisticTop.Cleanup();
    return true;
}


