#include <votca/ctp/qmpair.h>
#include <votca/ctp/qmtopology.h>

namespace votca { namespace ctp {

QMPair::QMPair(QMCrgUnit *crg1, QMCrgUnit *crg2, QMTopology * top)
  : std::pair<QMCrgUnit *, QMCrgUnit *>(crg1,crg2), _rate_12(0.),_rate_21(0.), _in_database(false)
{
    vec crg1nm,crg2nm;
    crg1nm =  crg1->GetCom();
    crg2nm =  crg2->GetCom();
    _r = top->BCShortestConnection(crg1nm, crg2nm);
    _crg2 = second;

    // check if PBC:
    vec d = crg2nm - crg1nm;
    if (abs(d - _r) > 1e-8) {
        _ghost = new QMCrgUnit();
	_ghost->copyCrgUnit(*crg2);
        vec displ = (_r - d);
        _ghost->shift(displ);
        _crg2 = _ghost;
    }
    else {
        _ghost=NULL;
    }
}

double QMPair::calcJeff2(){
    vector<double>::iterator itj=_Js.begin();
    double j=0.;
    for (;itj!= _Js.end(); itj++){
        j+=(*itj)*(*itj);
    }
    j/= double(_Js.size());
    return j;
}

}}