#include "easyjobserver.h"
#include <votca/csg/nblist.h>
#include <qmnblist.h>

EasyJObserver::EasyJObserver()
{}

EasyJObserver::~EasyJObserver()
{}

void EasyJObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

void EasyJObserver::EndCG()
{}

/// evaluate current conformation

void EasyJObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    _qmtop->Update(*top);
    QMNBList nblist;

    BeadList list1;
    Topology *toptmp = dynamic_cast<Topology*>(_qmtop);
    list1.Generate(*toptmp, "*");

    nblist.setCutoff(1.0);
    nblist.Generate(list1);

    for(NBList::iterator iter = nblist.getNBList()->begin();
        iter!=nblist.getNBList()->end();++iter) {
        QMBeadPair *pair = dynamic_cast<QMBeadPair*> (*iter);
        vector <double> Js = _qmtop->GetJCalc().GetJ(*pair->CrgUnit1(), *pair->CrgUnit2());
        cout << pair->CrgUnit1()->GetId() << " "
             << pair->CrgUnit1()->GetId() << " ";
        for(int i=0; i<Js.size(); +i)
            cout << Js[i] << " ";
        cout << endl;
    }
}

