#include "easyjobserver.h"
#include <votca/csg/nblist.h>
#include <qmnblist.h>
#include <votca/tools/histogramnew.h>

EasyJObserver::EasyJObserver()
{}

EasyJObserver::~EasyJObserver()
{}

void EasyJObserver::setNNnames(string  nnnames){
    Tokenizer tok(nnnames, " ;");
    Tokenizer::iterator itok = tok.begin();
    for (; itok!= tok.end(); ++itok){
        _nnnames.push_back(*itok);
    }
}
void EasyJObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

void EasyJObserver::EndCG()
{
    HistogramNew histo;
    histo.Initialize(*(std::min_element( _Js.begin(), _Js.end() )),*(std::max_element( _Js.begin(), _Js.end() ) ) , int(_Js.size() / 10 ));
    histo.ProcessRange(_Js.begin(),_Js.end());
    cout << histo <<endl;
}

/// evaluate current conformation

void EasyJObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    _qmtop->Update(*top);
    QMNBList &nblist = _qmtop->nblist();

    BeadList list1;
    Topology *toptmp = dynamic_cast<Topology*>(_qmtop);
    list1.Generate(*toptmp, "*");

    nblist.setCutoff(_cutoff);
    nblist.Generate(list1);

    for(QMNBList::iterator iter = nblist.begin();
        iter!=nblist.end();++iter) {
        CrgUnit *crg1 = (*iter)->first;
        CrgUnit *crg2 = (*iter)->second;
         if(MatchNNnames(crg1, crg2)){
            vector <double> Js = _qmtop->GetJCalc().GetJ(*crg1, *crg2);
            //cout << crg1->GetId() << " "
            //  << crg2->GetId() << " ";
            for(int i=0; i<Js.size(); ++i) _Js.push_back(Js[i]);
            


         }
    }
}

bool EasyJObserver::MatchNNnames(CrgUnit *crg1, CrgUnit* crg2){
    vector <string>::iterator its = _nnnames.begin();
    string namecrg = crg1->GetType()->GetName()+string(":")+crg2->GetType()->GetName();
    for ( ; its!= _nnnames.end(); ++its){
        if(wildcmp(its->c_str(), namecrg.c_str()) ){
            return true;
        }
    }
    return false;
}