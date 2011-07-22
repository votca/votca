
#include "md2qm_observer.h"
#include <votca/csg/nblist.h>
#include <votca/ctp/qmnblist.h>

MD2QMObserver::MD2QMObserver()
{}


MD2QMObserver::~MD2QMObserver()
{}

void MD2QMObserver::Initialize(QMTopology &qmtop, Property &opts)
{
    _qmtop = &qmtop;
    _save.Open(qmtop, _out);
    if(_save.FramesInDatabase() > 0)
        throw std::runtime_error("Database already contains frame information. "
                "Appending a database is not supported yet");
}

void MD2QMObserver::BeginCG(Topology *top, Topology *top_atom)
{
    _qmtop->Initialize(*top);
}

void MD2QMObserver::EvalConfiguration(Topology *top, Topology *top_atom)
{
    if(_qmtop->getDatabaseId() != 0)
        throw std::runtime_error("writing several frames to state file not yet supported, please use --nframes=1");
    _qmtop->Update(*top);
    _save.WriteFrame();
}

void MD2QMObserver::EndCG()
{
    _save.Close();
}

