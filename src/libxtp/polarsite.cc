/*
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/polarsite.h>
#include <fstream>


namespace votca { namespace xtp {


void PolarSite::ImportFrom(PolarSite *templ, string tag) {

    _pos = templ->getPos();

    if (tag == "basic") {
        _Qs[0] = templ->getQs(-1);
        _Qs[1] = templ->getQs(0);
        _Qs[2] = templ->getQs(1);
        _rank  = templ->_rank;
        _Ps    = templ->_Ps;
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
        _Ps   = templ->_Ps;
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

    // Rotate multipoles into global frame >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    for (int state = -1; state < 2; state++) {

        // Any multipoles for this charge state available?
        if (_Qs[state+1].size() < 1) { continue; }

        //   0    1    2    3    4    5    6    7    8    9    10   ...
        //   Q00  Q10  Q1c  Q1s  Q20  Q21c Q21s Q22c Q22s Q30  Q31c ...

        matrix R = rot;
        matrix R_T = matrix(R.getRow(0),R.getRow(1),R.getRow(2));

        // Transform dipole moment into global frame
        if (_Qs[state+1].size() > 1) {

            double Qz = _Qs[state+1][1];
            double Qx = _Qs[state+1][2];
            double Qy = _Qs[state+1][3];

            vec d = vec(Qx, Qy, Qz);
            d = R * d;

            _Qs[state+1][1] = d.getZ();
            _Qs[state+1][2] = d.getX();
            _Qs[state+1][3] = d.getY();
        }

        // Transform quadrupole moment into global frame
        if (_Qs[state+1].size() > 4) {

            double Qzz =      _Qs[state+1][4];
            double Qxx = -0.5*_Qs[state+1][4] + 0.5*sqrt(3)*_Qs[state+1][7];
            double Qyy = -0.5*_Qs[state+1][4] - 0.5*sqrt(3)*_Qs[state+1][7];

            double Qxy =  0.5*sqrt(3)*_Qs[state+1][8];
            double Qxz =  0.5*sqrt(3)*_Qs[state+1][5];
            double Qyz =  0.5*sqrt(3)*_Qs[state+1][6];

            matrix Q = matrix(vec(Qxx,Qxy,Qxz),
                              vec(Qxy,Qyy,Qyz),
                              vec(Qxz,Qyz,Qzz));

            matrix Q_Global  = R * Q * R_T;

            /* if (this->getId() == 1) {
                cout << endl;
                cout << "  " << Q_Global.get(0,0);
                cout << "  " << Q_Global.get(0,1);
                cout << "  " << Q_Global.get(0,2);
                cout << endl;
                cout << "  " << Q_Global.get(1,0);
                cout << "  " << Q_Global.get(1,1);
                cout << "  " << Q_Global.get(1,2);
                cout << endl;
                cout << "  " << Q_Global.get(2,0);
                cout << "  " << Q_Global.get(2,1);
                cout << "  " << Q_Global.get(2,2);
                cout << endl;
            }                                   */

            _Qs[state+1][4] =               Q_Global.get(2,2);  // Q20
            _Qs[state+1][5] = 2 / sqrt(3) * Q_Global.get(0,2);  // Q21c
            _Qs[state+1][6] = 2 / sqrt(3) * Q_Global.get(1,2);  // Q21s
            _Qs[state+1][7] = 1 / sqrt(3) *(Q_Global.get(0,0)
                                          - Q_Global.get(1,1)); // Q22c
            _Qs[state+1][8] = 2 / sqrt(3) * Q_Global.get(0,1);  // Q22s
        }
    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


}

void PolarSite::Translate(const vec &shift) {

    _pos += shift;

}

void PolarSite::Charge(int state) {

    int idx = state + 1;

    // Adjust polarizability to charge state
        P1 = _Ps[idx];

    // Adjust multipole moments to charge state
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

void PolarSite::ChargeDelta(int state1, int state2) {

    int idx1 = state1 + 1;
    int idx2 = state2 + 1;

    // Adjust polarizability to charge state
        P1 = 0.0;

    // Adjust multipole moments to charge state
        Q00 = _Qs[idx2][0] - _Qs[idx1][0];

    if (_rank > 0) {
        Q1z = _Qs[idx2][1] - _Qs[idx1][1];   // |
        Q1x = _Qs[idx2][2] - _Qs[idx1][2];   // |-> NOTE: order z - x - y
        Q1y = _Qs[idx2][3] - _Qs[idx1][3];   // |
    }
    if (_rank > 1) {
        Q20  = _Qs[idx2][4] - _Qs[idx1][4];
        Q21c = _Qs[idx2][5] - _Qs[idx1][5];
        Q21s = _Qs[idx2][6] - _Qs[idx1][6];
        Q22c = _Qs[idx2][7] - _Qs[idx1][7];
        Q22s = _Qs[idx2][8] - _Qs[idx1][8];
    }
}

void PolarSite::Induce(double wSOR) {

    U1_Hist.push_back( vec(U1x,U1y,U1z) );
    
    U1x = (1 - wSOR) * U1x + wSOR * ( - P1 * (FPx + FUx) ); // OVERRIDE
    U1y = (1 - wSOR) * U1y + wSOR * ( - P1 * (FPy + FUy) ); // OVERRIDE
    U1z = (1 - wSOR) * U1z + wSOR * ( - P1 * (FPz + FUz) ); // OVERRIDE
}

void PolarSite::InduceDirect() {

    U1_Hist.push_back( vec(0.,0.,0.) );
    U1x =  - P1 * FPx; // OVERRIDE
    U1y =  - P1 * FPy; // OVERRIDE
    U1z =  - P1 * FPz; // OVERRIDE
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


void PolarSite::WriteXyzLine(FILE *out, vec &shift, string format) {

    //double int2ext = 1.0;

    vec pos = _pos + shift;

    /*if (format == "gaussian") {
        int2ext = 10.;
    }
    else {
        int2ext = 10.;
    }*/

    fprintf(out, "%-2s %+4.9f %+4.9f %+4.9f \n",
            _name.c_str(),
            pos.getX()*10, pos.getY()*10, pos.getZ()*10);
}


void PolarSite::WriteChkLine(FILE *out, vec &shift, bool split_dpl,
                             string format, double spacing) {

    vec pos = _pos + shift;

    string unit = "";

    if (format == "xyz") {
        unit = "angstrom";
    }
    else if (format == "gaussian") {
        unit = "angstrom";
    }

    // Take care of unit conversion
    double int2ext;

    if (unit == "nanometer") {
        int2ext = 1.;
    }
    else if (unit == "angstrom") {
        int2ext = 10.;
    }
    else if (unit == "bohr") {
        assert(false);
    } else {
        throw std::runtime_error( "Unit has to be nanometer, angstrom or bohr");
    }

    if (format == "xyz") {
        fprintf(out, "%2s ", _name.c_str());
    }

    // Print charge line
    fprintf(out, "%+4.9f %+4.9f %+4.9f %+4.7f \n",
            pos.getX()*int2ext,
            pos.getY()*int2ext,
            pos.getZ()*int2ext,
            Q00);
    

    // Split dipole moment onto charges (if desired)
    if (split_dpl) {

        vec tot_dpl = vec(U1x,U1y,U1z);

        if (_rank > 0) { tot_dpl += vec(Q1x,Q1y,Q1z); }

        matrix::eigensystem_t EIGEN;

        if (_rank == 2) {
            tot_dpl += vec(Q1x,Q1y,Q1z);
            cout << endl
                 << "WARNING: Quadrupoles are not split onto point charges."
                 << endl;

            int state = 0;

            double Qzz =      _Qs[state+1][4];
            double Qxx = -0.5*_Qs[state+1][4] + 0.5*sqrt(3)*_Qs[state+1][7];
            double Qyy = -0.5*_Qs[state+1][4] - 0.5*sqrt(3)*_Qs[state+1][7];

            double Qxy =  0.5*sqrt(3)*_Qs[state+1][8];
            double Qxz =  0.5*sqrt(3)*_Qs[state+1][5];
            double Qyz =  0.5*sqrt(3)*_Qs[state+1][6];

            matrix Q = matrix(vec(Qxx,Qxy,Qxz),
                              vec(Qxy,Qyy,Qyz),
                              vec(Qxz,Qyz,Qzz));

            
            Q.SolveEigensystem(EIGEN);


        }

        double a        = spacing;
        double mag_d    = abs(tot_dpl);
        vec    dir_d_0  = tot_dpl.normalize();
        vec    dir_d    = dir_d_0.normalize();
        vec    A        = pos + 0.5 * a * dir_d;
        vec    B        = pos - 0.5 * a * dir_d;
        double qA       = mag_d / a;
        double qB       = - qA;

        if (format == "xyz") {
            fprintf(out, " A ");
        }
        fprintf(out, "%+4.9f %+4.9f %+4.9f %+4.7f \n",
                A.getX()*int2ext,
                A.getY()*int2ext,
                A.getZ()*int2ext,
                qA);

        if (format == "xyz") {
            fprintf(out, " B ");
        }
        fprintf(out, "%+4.9f %+4.9f %+4.9f %+4.7f \n",
                B.getX()*int2ext,
                B.getY()*int2ext,
                B.getZ()*int2ext,
                qB);

        if (format == "xyz" && _rank == 2) {
            vec D1 = pos + 0.5 * a * EIGEN.eigenvecs[0];
            vec D2 = pos - 0.5 * a * EIGEN.eigenvecs[0];
            vec E1 = pos + 0.5 * a * EIGEN.eigenvecs[1];
            vec E2 = pos - 0.5 * a * EIGEN.eigenvecs[1];
            vec F1 = pos + 0.5 * a * EIGEN.eigenvecs[2];
            vec F2 = pos - 0.5 * a * EIGEN.eigenvecs[2];
            fprintf(out, " D %+4.9f %+4.9f %+4.9f \n",
                    D1.getX()*int2ext,
                    D1.getY()*int2ext,
                    D1.getZ()*int2ext);
            fprintf(out, " D %+4.9f %+4.9f %+4.9f \n",
                    D2.getX()*int2ext,
                    D2.getY()*int2ext,
                    D2.getZ()*int2ext);
            fprintf(out, " E %+4.9f %+4.9f %+4.9f \n",
                    E1.getX()*int2ext,
                    E1.getY()*int2ext,
                    E1.getZ()*int2ext);
            fprintf(out, " E %+4.9f %+4.9f %+4.9f \n",
                    E2.getX()*int2ext,
                    E2.getY()*int2ext,
                    E2.getZ()*int2ext);
            fprintf(out, " F %+4.9f %+4.9f %+4.9f \n",
                    F1.getX()*int2ext,
                    F1.getY()*int2ext,
                    F1.getZ()*int2ext);
            fprintf(out, " F %+4.9f %+4.9f %+4.9f \n",
                    F2.getX()*int2ext,
                    F2.getY()*int2ext,
                    F2.getZ()*int2ext);
        }



    }

}

void PolarSite::PrintInfo(std::ostream &out) {

    cout << "MPOLE " << this->getId() << " " << this->getName()
         << " " << _pos.getX() << " " << _pos.getY() << " " << _pos.getZ()
         << endl;

    for (int i = -1; i < 2; ++i) {
        cout << "    STATE " << i << "   ";
        for (unsigned int j = 0; j < this->getQs(i).size(); j++) {
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
    cout << "    Alpha " << P1 << endl;
}

void PolarSite::PrintInfoVisual(FILE *out) {

    int id = this->getId();
    string el = this->getName();
    vec position = this->getPos();
    fprintf(out, "MPOLE %5d %4s %3.8f %3.8f %3.8f ", id, el.c_str(), position.getX(), position.getY(), position.getZ());

    for (unsigned int i = 0; i < this->U1_Hist.size(); ++i) {
        vec uInd = U1_Hist[i];
        fprintf(out, " %3.8f %3.8f %3.8f ", uInd.getX(), uInd.getY(), uInd.getZ() );
    }
    fprintf(out, " U1 %3.8f %3.8f %3.8f ", U1x, U1y, U1z );
    fprintf(out, " FP %3.8f %3.8f %3.8f ", FPx, FPy, FPz );
    fprintf(out, " FU %3.8f %3.8f %3.8f ", FUx, FUy, FUz );

    fprintf(out, " \n");
}

void PolarSite::PrintPDB(FILE *out, vec shift) {

    int id = this->getId() % 100000;
    string name =  this->getName();
    name.resize(3);
    string resname = "RSD";
    resname.resize(3);
    int resnr = 0;
    vec position = this->getPos() + shift;

    fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   "
              "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%4.7f\n",
         id,                    // Atom serial number           %5d
         name.c_str(),          // Atom name                    %4s
         " ",                   // alternate location indicator.%1s
         resname.c_str(),       // Residue name.                %3s
         "A",                   // Chain identifier             %1s
         resnr,                 // Residue sequence number      %4d
         " ",                   // Insertion of residues.       %1s
         position.getX()*10,    // X in Angstroms               %8.3f
         position.getY()*10,    // Y in Angstroms               %8.3f
         position.getZ()*10,    // Z in Angstroms               %8.3f
         1.0,                   // Occupancy                    %6.2f
         0.0,                   // Temperature factor           %6.2f
         " ",                   // Segment identifier           %4s
         name.c_str(),          // Element symbol               %2s
         " ",                   // Charge on the atom.          %2s
         this->getQ00()
         );    
}


}}

