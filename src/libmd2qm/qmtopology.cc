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
    BeadContainer::iterator iter;

    CopyTopologyData(&cg_top);

    for(iter=_beads.begin(); iter!=_beads.end(); ++iter) {
        Bead *b = *iter;
        QMBead *qmbead = (QMBead *)b ;// TODO: why does dynamic_cast<QMBead*>(b); not work

        // TODO: initialize crgunit data here
    }
}

void QMTopology::Update(Topology& cg_top)
{
    BeadContainer::iterator iter;
    BeadContainer::iterator iter_cg;

    assert(cg_top.Beads().size() == _beads.size());

    _box = cg_top.getBox();
    _time = cg_top.getTime();
    _step = cg_top.getStep();

    iter_cg = cg_top.Beads().begin();
    for(iter=_beads.begin(); iter!=_beads.end(); ++iter) {
        (*iter)->setPos((*iter_cg)->getPos());
        (*iter)->setU((*iter_cg)->getU());
        (*iter)->setV((*iter_cg)->getV());
        (*iter)->setW((*iter_cg)->getW());
    }
}

void QMTopology::LoadListCharges(const string &file)
{

}