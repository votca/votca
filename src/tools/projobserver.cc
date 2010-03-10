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
   
    string framedir=string("frame")+lexical_cast<string>(top->getStep()) +string("/") ;
    mkdir(framedir.c_str(),0755);
    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;

        /* this next bit will be substituted by the CalcJ writer for projection
         */
         
        if(MatchNNnames(crg1, crg2)){
            _qmtop->GetJCalc().WriteProJ(*crg1, *crg2,framedir);
        }

    }
}

bool ProJObserver::MatchNNnames(CrgUnit *crg1, CrgUnit* crg2){
    vector <string>::iterator its = _nnnames.begin();
    string namecrg = crg1->getType()->GetName()+string(":")+crg2->getType()->GetName();
    for ( ; its!= _nnnames.end(); ++its){
        if(wildcmp(its->c_str(), namecrg.c_str()) ){
            return true;
        }
    }
    return false;
}