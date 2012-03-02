#include <votca/ctp/polarsite.h>
#include <fstream>


namespace votca { namespace ctp {


void PolarSite::ImportFrom(PolarSite *templ, string tag) {

    _pos = templ->getPos();

    if (tag == "basic") {
        _Qs[0] = templ->getQs(-1);
        _Qs[1] = templ->getQs(0);
        _Qs[2] = templ->getQs(1);
        _rank  = templ->_rank;
        alpha = templ->alpha;
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
        alpha = templ->alpha;
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

//    // Rotate multipoles into global frame >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
//    for (int state = -1; state < 2; state++) {
//
//        if (_Qs[state+1].size() < 1) { continue; }
//
//        //   0     1    2    3    4    5    6    7    8    9    10   ...
//        //   Q00   Q10  Q1c  Q1s  Q20  Q21c Q21s Q22c Q22s Q30  Q31c ...
//
//        matrix R = rot;
//        matrix R_T = matrix(R.getRow(0),R.getRow(1),R.getRow(2));
//        matrix Unit = R * R_T;
//        if (this->getId() == 1) {
//            cout << endl;
//            cout << " 00 " << Unit.get(0,0);
//            cout << " 01 " << Unit.get(0,1);
//            cout << " 02 " << Unit.get(0,2);
//            cout << endl;
//            cout << " 10 " << Unit.get(1,0);
//            cout << " 11 " << Unit.get(1,1);
//            cout << " 12 " << Unit.get(1,2);
//            cout << endl;
//            cout << " 20 " << Unit.get(2,0);
//            cout << " 21 " << Unit.get(2,1);
//            cout << " 22 " << Unit.get(2,2);
//            cout << endl;
//        }
//
//        // Transform dipole moment into global frame
//        if (_Qs[state+1].size() > 1) {
//
//            double Qz = _Qs[state+1][1];
//            double Qx = _Qs[state+1][2];
//            double Qy = _Qs[state+1][3];
//
//            vec d = vec(Qx, Qy, Qz);
//            d = R * d;
//
//            _Qs[state+1][1] = d.getZ();
//            _Qs[state+1][2] = d.getX();
//            _Qs[state+1][3] = d.getY();
//        }
//
//
//        // Transform quadrupole moment into global frame
//        if (_Qs[state+1].size() > 4) {
//
//            double Qzz =      _Qs[state+1][4];
//            double Qxx = -0.5*_Qs[state+1][4] + 0.5*sqrt(3)*_Qs[state+1][7];
//            double Qyy = -0.5*_Qs[state+1][4] - 0.5*sqrt(3)*_Qs[state+1][7];
//
//            double Qxy =  0.5*sqrt(3)*_Qs[state+1][8];
//            double Qxz =  0.5*sqrt(3)*_Qs[state+1][5];
//            double Qyz =  0.5*sqrt(3)*_Qs[state+1][6];
//
//            matrix Q = matrix(vec(Qxx,Qxy,Qxz),
//                              vec(Qxy,Qyy,Qyz),
//                              vec(Qxz,Qyz,Qzz));
//
//            matrix Q_Global  = R * Q * R_T;
//
//            if (this->getId() == 1) {
//                cout << endl;
//                cout << "  " << Q_Global.get(0,0);
//                cout << "  " << Q_Global.get(0,1);
//                cout << "  " << Q_Global.get(0,2);
//                cout << endl;
//                cout << "  " << Q_Global.get(1,0);
//                cout << "  " << Q_Global.get(1,1);
//                cout << "  " << Q_Global.get(1,2);
//                cout << endl;
//                cout << "  " << Q_Global.get(2,0);
//                cout << "  " << Q_Global.get(2,1);
//                cout << "  " << Q_Global.get(2,2);
//                cout << endl;
//            }
//
//            _Qs[state+1][4] = Q_Global.get(2,2);                                  // Q20
//            _Qs[state+1][5] = 2/sqrt(3) * Q_Global.get(0,2);                      // Q21c
//            _Qs[state+1][6] = 2/sqrt(3) * Q_Global.get(1,2);                      // Q21s
//            _Qs[state+1][7] = 1/sqrt(3) * (Q_Global.get(0,0) - Q_Global.get(1,1));// Q22c
//            _Qs[state+1][8] = 2/sqrt(3) * Q_Global.get(0,1);                      // Q22s
//        }
//    }
//    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


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

void PolarSite::Induce(double wSOR) {

    U1_Hist.push_back( vec(U1x,U1y,U1z) );
    
    U1x = (1 - wSOR) * U1x + wSOR * ( - alpha * (FPx + FUx) );
    U1y = (1 - wSOR) * U1y + wSOR * ( - alpha * (FPy + FUy) );
    U1z = (1 - wSOR) * U1z + wSOR * ( - alpha * (FPz + FUz) );    
}

void PolarSite::InduceDirect() {

    U1_Hist.push_back( vec(0.,0.,0.) );
    U1x =  - alpha * FPx;
    U1y =  - alpha * FPy;
    U1z =  - alpha * FPz;
}

double PolarSite::HistdU() {

    vec dU = vec(U1x, U1y, U1z) - U1_Hist.back();
    return abs(dU)/abs(U1_Hist.back());
}



void PolarSite::Depolarize() {

    // Zero out induced moments
    U1x = U1y = U1z = 0.0;
    U1_Hist.clear();

    // Zero out fields
    FPx = FPy = FPz = 0.0;
    FUx = FUy = FUz = 0.0;
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

void PolarSite::PrintInfoInduce(std::ostream &out) {

    cout << endl;
    cout << "MPOLE " << this->getId() << " " << this->getName()
         << " " << _pos.getX() << " " << _pos.getY() << " " << _pos.getZ()
         << endl;

    cout << "    FIELD P " << FPx << " " << FPy << " " << FPz << endl;
    cout << "    FIELD U " << FUx << " " << FUy << " " << FUz << endl;
    cout << "    DIPOLE U " << U1x << " " << U1y << " " << U1z << endl;
    cout << "    Q0 " << Q00 << endl;
    cout << "    Q1 " << Q1x << " " << Q1y << " " << Q1z << endl;
    cout << "    Q2 " << Q20 << " " << Q21c << " " << Q21s << " " << Q22c << " " << Q22s << endl;
    cout << "    Alpha " << alpha << endl;
}

void PolarSite::PrintInfoVisual(FILE *out) {

    int id = this->getId();
    string el = this->getName();
    vec position = this->getPos();
    fprintf(out, "MPOLE %5d %4s %3.8f %3.8f %3.8f ", id, el.c_str(), position.getX(), position.getY(), position.getZ());

    for (int i = 0; i < this->U1_Hist.size(); ++i) {
        vec uInd = U1_Hist[i];
        fprintf(out, " %3.8f %3.8f %3.8f ", uInd.getX(), uInd.getY(), uInd.getZ() );
    }

    fprintf(out, " \n");
}


}}

