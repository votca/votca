#include <votca/ctp/polarsite.h>

namespace votca { namespace ctp {


void PolarSite::ImportFrom(PolarSite *templ, string tag) {

    _pos = templ->getPos();

    if (tag == "basic") {
        _Qs[0] = templ->getQs(-1);
        _Qs[1] = templ->getQs(0);
        _Qs[2] = templ->getQs(1);
        _rank  = templ->_rank;
    }

    if (tag == "full") {
        _locX = templ->_locX;
        _locY = templ->_locY;
        _locZ = templ->_locZ;
        _top  = templ->_top;
        _seg  = templ->_seg;
        _frag = templ->_frag;
        _rank = templ->_rank; 
        _Qs   = templ->_Qs;
        _id   = templ->_id;
        _name = templ->_name;
    }
}

void PolarSite::Rotate(const matrix &rot, const vec &refPos) {

    vec dir = _pos - refPos;
    dir = rot * dir;
    _pos = refPos + dir;

    _locX = rot.getCol(0);
    _locY = rot.getCol(1);
    _locZ = rot.getCol(2);

}

void PolarSite::Translate(const vec &shift) {

    _pos += shift;

}

void PolarSite::Charge(int state) {

    int idx = state + 1;

        Q00 = _Qs[idx][0];

    if (_rank > 0) {
        Q1z = _Qs[idx][1];   // |
        Q1x = _Qs[idx][2];   // |-> NOTE: order z - x - y
        Q1y = _Qs[idx][3];   // |
    }
    if (_rank > 1) {
        Q20  = _Qs[idx][4];
        Q21c = _Qs[idx][5];
        Q21s = _Qs[idx][6];
        Q22c = _Qs[idx][7];
        Q22s = _Qs[idx][8];
    }
}


void PolarSite::PrintInfo(std::ostream &out) {

    cout << "MPOLE " << this->getId() << " " << this->getName()
         << " " << _pos.getX() << " " << _pos.getY() << " " << _pos.getZ()
         << endl;

    for (int i = -1; i < 2; ++i) {
        cout << "    STATE " << i << "   ";
        for (int j = 0; j < this->getQs(i).size(); j++) {
            cout << " " << this->getQs(i)[j];
        }
        cout << endl;
    }

    cout << "    XAXIS " << _locX.getX() << " " << _locX.getY() << " " << _locX.getZ() << endl;
    cout << "    YAXIS " << _locY.getX() << " " << _locY.getY() << " " << _locY.getZ() << endl;
    cout << "    ZAXIS " << _locZ.getX() << " " << _locZ.getY() << " " << _locZ.getZ() << endl;


}


}}

