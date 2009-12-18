#include "qmpair.h"
#include "qmtopology.h"


QMPair::QMPair(CrgUnit *crg1, CrgUnit *crg2, QMTopology * top):std::pair<CrgUnit *, CrgUnit *>(crg1,crg2)
{
    _r = top->BCShortestConnection(crg1->GetCom(), crg2->GetCom());
    _dist = abs(_r);

    // check if PBC:
    vec d = crg2->GetCom() - crg1->GetCom();
    if (abs(d) !=  _dist){
        _ghost = new CrgUnit();
	_ghost->copyCrgUnit(*crg2);
        vec displ = _r - d;
        _ghost->shift(displ);
        this ->second = _ghost;
    }
    else {
        _ghost=NULL;
    }
}
