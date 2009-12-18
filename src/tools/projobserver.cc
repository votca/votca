#include "projobserver.h"
#include <votca/csg/nblist.h>
#include <votca/csg/trajectorywriter.h>
#include <qmnblist.h>

ProJObserver::ProJObserver()
{}

ProJObserver::~ProJObserver()
{}

void ProJObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

void ProJObserver::EndCG()
{}

void ProJObserver::setNNnames(string  nnnames){
    Tokenizer tok(nnnames, " ;");
    Tokenizer::iterator itok = tok.begin();
    for (; itok!= tok.end(); ++itok){
        _nnnames.push_back(*itok);
    }
}


/// evaluate current conformation



void ProJObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    _qmtop->Update(*top);
    
    BeadList list1;
    Topology *toptmp = dynamic_cast<Topology*>(_qmtop);
    list1.Generate(*toptmp, "*");
    QMNBList &nblist =(_qmtop->nblist());
    nblist.setCutoff(_cutoff);
    nblist.Generate(list1);
    TrajectoryWriter *writer;
    writer = TrjWriterFactory().Create(".pdb");
    if(writer == NULL)
        throw runtime_error("output format not supported: .pdb");
    string framedir=string("frame")+lexical_cast<string>(top->getStep())  ;

    mkdir(framedir.c_str(),0755);
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
        if(MatchNNnames(crg1, crg2)){
            Topology atoms;
            _qmtop->AddAtomisticBeads(crg1,&atoms);
            _qmtop->AddAtomisticBeads(crg2,&atoms);
        
            ///write the topo somehow now.
            string nameout =framedir + string("/")+lexical_cast<string>(crg1->GetId())+ string("and")
                + lexical_cast<string>(crg2->GetId()) + ".pdb";

        
            writer->Open(nameout);
            writer->Write(&atoms);
            writer->Close();
        }

    }
}

bool ProJObserver::MatchNNnames(CrgUnit *crg1, CrgUnit* crg2){
    vector <string>::iterator its = _nnnames.begin();
    string namecrg = crg1->GetType()->GetName()+string(":")+crg2->GetType()->GetName();
    for ( ; its!= _nnnames.end(); ++its){
        if(wildcmp(its->c_str(), namecrg.c_str()) ){
            return true;
        }
    }
    return false;
}