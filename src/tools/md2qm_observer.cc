#include "md2qm_observer.h"
#include <votca/csg/nblist.h>
#include <qmnblist.h>

MD2QMObserver::MD2QMObserver()
{}


MD2QMObserver::~MD2QMObserver()
{}

void MD2QMObserver::Initialize(QMTopology &qmtop, Property &opts)
{
    _qmtop = &qmtop;
    _save.Open(qmtop, _out, 'w');
}

void MD2QMObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

/// evaluate current conformation

void MD2QMObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{

    _qmtop->Update(*top);
    QMNBList &nblist = _qmtop->nblist();

    BeadList list1;
    Topology *toptmp = dynamic_cast<Topology*>(_qmtop);
    list1.Generate(*toptmp, "*");

    nblist.setCutoff(_cutoff);
    nblist.Generate(list1);
    _save.Save();

}

void MD2QMObserver::EndCG()
{
    print_nbs_to_file(_qmtop->nblist());
    _save.Close();
}

void MD2QMObserver::print_nbs_to_file(QMNBList &nblist){
    ofstream out_nbl;
    out_nbl.open("nbl_votca.res");
    if(out_nbl!=0){
        out_nbl << "Neighbours, J(0), J_eff, rate, r_ij, abs(r_ij) [nm]" << endl;
        QMNBList::iterator iter;
        for (iter = nblist.begin(); iter != nblist.end(); ++iter) {
            out_nbl << "(" << (*iter)->first->getId() << "," << (*iter)->second->getId() << "): ";
            out_nbl << 0.0 << " " << 0.0 << " " << 0.0 << " ";
            out_nbl << (*iter)->r().getX() << " " << (*iter)->r().getY() << " " << (*iter)->r().getZ() << " ";
            out_nbl << " " << (*iter)->dist() << endl;
        }
    }
    out_nbl.close();
}

