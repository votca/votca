#include "qmtopology.h"

QMTopology::QMTopology()
{
    _nblist = NULL;
}

QMTopology::~QMTopology()
{
    if(_nblist)
        delete _nblist;
    _nblist = NULL;
}

void QMTopology::Initialize(Topology& cg_top)
{
    BeadContainer::iteratur iter;

    CopyTopologyData(&cg_top);

    for(iter=_beads.begin(); iter!=_beads,end(); ++iter) {
        QMBead *qmbead = dynamic_cast<QMBead*>*iter;

        // TODO: initialize crgunit data here
    }
}

void QMTopology::Update(Topology& cg_top)
{
    BeadContainer::iteratur iter;
    BeadContainer::iteratur iter_cg;

    assert(cg_top._beads().size() == _beads.size());

    _box = cg_top._box;
    _time = cg_top._time;
    _frame = cg_top._frame;

    iter_cg = cg_top._beads.begin();
    for(iter=_beads.begin(); iter!=_beads,end(); ++iter) {
        (*iter)->setPos((*iter_cg)->getPos());
        (*iter)->setU((*iter_cg)->getU());
        (*iter)->setV((*iter_cg)->getV());
        (*iter)->setW((*iter_cg)->getW());
    }
}