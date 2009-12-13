#include "easyjobserver.h"
#include <votca/csg/nblist.h>

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
}

