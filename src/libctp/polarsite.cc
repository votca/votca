#include <votca/ctp/polarsite.h>

namespace votca { namespace ctp {


void PolarSite::ImportFrom(PolarSite *templ, string tag) {

    _pos = templ->getPos();

    if (tag == "basic") {
        _Qs[0] = templ->getQs(-1);
        _Qs[1] = templ->getQs(0);
        _Qs[2] = templ->getQs(1);
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

