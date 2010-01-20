#include "qmpair.h"
#include "qmtopology.h"


QMPair::QMPair(CrgUnit *crg1, CrgUnit *crg2, QMTopology * top):std::pair<CrgUnit *, CrgUnit *>(crg1,crg2)
{
    vec crg1nm,crg2nm;
    crg1nm =  crg1->GetCom()*(RA/10.);
    crg2nm =  crg2->GetCom()*(RA/10.);
    _r = top->BCShortestConnection(crg1nm, crg2nm);
    _dist = abs(_r);

    // check if PBC:
    vec d = crg2nm - crg1nm;
    if (abs(d) !=  _dist){
        _ghost = new CrgUnit();
	_ghost->copyCrgUnit(*crg2);
        vec displ = (_r - d)*(10./RA);
        _ghost->shift(displ);
        this ->second = _ghost;
    }
    else {
        _ghost=NULL;
    }
    /// move r and dist back to bohrs:
    _r *= 10./RA;
    _dist *= 10./RA;
}

void QMPair::setJeff(vector <double> js){
    vector<double>::iterator itj=js.begin();
    double j=0.;
    for (;itj!= js.end(); itj++){
        j+=(*itj)*(*itj);
    }
    j/= double(js.size());
    _Jeff = j;
}