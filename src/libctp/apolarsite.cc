#include <votca/ctp/apolarsite.h>
#include <fstream>
#include <string>


namespace votca { namespace ctp {


void APolarSite::ImportFrom(APolarSite *templ, string tag) {

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

void APolarSite::Rotate(const matrix &rot, const vec &refPos) {

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

        // Transform polarizability tensor into global frame
        matrix P_Global = R * _Ps[state+1] * R_T;
        _Ps[state+1] = P_Global;

    }
    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


}

void APolarSite::Translate(const vec &shift) {

    _pos += shift;

}

void APolarSite::Charge(int state) {

    int idx = state + 1;

    // Adjust polarizability to charge state
        Pxx = _Ps[idx].get(0,0);
        Pxy = _Ps[idx].get(0,1);
        Pxz = _Ps[idx].get(0,2);
        Pyy = _Ps[idx].get(1,1);
        Pyz = _Ps[idx].get(1,2);
        Pzz = _Ps[idx].get(2,2);

    // Calculate principal axes
        matrix polarity = matrix( vec(Pxx,Pxy,Pxz),
                                  vec(Pxy,Pyy,Pyz),
                                  vec(Pxz,Pyz,Pzz) );
        matrix::eigensystem_t eigensystem_polarity;
        polarity.SolveEigensystem(eigensystem_polarity);

        pax         = eigensystem_polarity.eigenvecs[0];
        pay         = eigensystem_polarity.eigenvecs[1];
        paz         = eigensystem_polarity.eigenvecs[2];
        eigenpxx    = eigensystem_polarity.eigenvalues[0];
        eigenpyy    = eigensystem_polarity.eigenvalues[1];
        eigenpzz    = eigensystem_polarity.eigenvalues[2];

//        eigendamp =  (Pxx > Pyy) ?
//                         ((Pxx > Pzz) ? Pxx : Pzz)
//                       : ((Pyy > Pzz) ? Pyy : Pzz);
        eigendamp =  (eigenpxx > eigenpyy) ?
                         ((eigenpxx > eigenpzz) ? eigenpxx : eigenpzz)
                       : ((eigenpyy > eigenpzz) ? eigenpyy : eigenpzz);

        // eigendamp = 10.;

        //cout << endl << "eigensystem ...";
        //cout << endl << pax << " --- " << eigenpxx;
        //cout << endl << pay << " --- " << eigenpyy;
        //cout << endl << paz << " --- " << eigenpzz;


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

void APolarSite::ChargeDelta(int state1, int state2) {

    int idx1 = state1 + 1;
    int idx2 = state2 + 1;

    // Adjust polarizability to charge state
        Pxx = 0.0;
        Pxy = 0.0;
        Pxz = 0.0;
        Pyy = 0.0;
        Pyz = 0.0;
        Pzz = 0.0;

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


double APolarSite::getProjP(vec &dir) {


    // Correct
//    double alpha_proj = fabs(dir * pax) * fabs(dir * pax) * eigenpxx
//         + fabs(dir * pay) * fabs(dir * pay) * eigenpyy
//         + fabs(dir * paz) * fabs(dir * paz) * eigenpzz;

    // Mix
//    alpha_proj = 0.5* (1./3. * (Pxx+Pyy+Pzz)) + 0.5*alpha_proj;

    // Wrong
//    double alpha_proj = fabs(dir * pax) * eigenpxx
//         + fabs(dir * pay) * eigenpyy
//         + fabs(dir * paz) * eigenpzz;
    //cout << endl << _name << " " << alpha_proj;

    // Max
    double alpha_proj;

    if (Pxx > Pyy) {
        if (Pxx > Pzz) {
            alpha_proj = Pxx;
        }
        else {
            alpha_proj = Pzz;
        }
    }
    else {
        if (Pyy > Pzz) {
            alpha_proj = Pyy;
        }
        else {
            alpha_proj = Pzz;
        }
    }

    return this->eigendamp;
}

void APolarSite::Induce(double wSOR) {

    U1_Hist.push_back( vec(U1x,U1y,U1z) );

    U1x = (1 - wSOR) * U1x + wSOR * ( - Pxx * (FPx + FUx) - Pxy * (FPy + FUy) - Pxz * (FPz + FUz) ); // OVERRIDE
    U1y = (1 - wSOR) * U1y + wSOR * ( - Pxy * (FPx + FUx) - Pyy * (FPy + FUy) - Pyz * (FPz + FUz) ); // OVERRIDE
    U1z = (1 - wSOR) * U1z + wSOR * ( - Pxz * (FPx + FUx) - Pyz * (FPy + FUy) - Pzz * (FPz + FUz) ); // OVERRIDE
}

void APolarSite::InduceDirect() {

    U1_Hist.push_back( vec(0.,0.,0.) );
    U1x =  - Pxx * FPx - Pxy * FPy - Pxz * FPz; // OVERRIDE
    U1y =  - Pxy * FPx - Pyy * FPy - Pyz * FPz; // OVERRIDE
    U1z =  - Pxz * FPx - Pyz * FPy - Pzz * FPz; // OVERRIDE
}

double APolarSite::HistdU() {

    vec dU = vec(U1x, U1y, U1z) - U1_Hist.back();
    return abs(dU)/abs(U1_Hist.back());
}



void APolarSite::Depolarize() {

    // Zero out induced moments
    U1x = U1y = U1z = 0.0;
    U1_Hist.clear();

    // Zero out fields
    FPx = FPy = FPz = 0.0;
    FUx = FUy = FUz = 0.0;
}



void APolarSite::PrintInfo(std::ostream &out) {

    printf("\n%2s %2d POS %+2.3f %+2.3f %+2.3f "
           "POL %+2.3f %+2.3f %+2.3f %+2.3f %+2.3f %+2.3f ",
            _name.c_str(), _id, _pos.getX(), _pos.getY(), _pos.getZ(),
            Pxx*1000, Pxy*1000, Pxz*1000, Pyy*1000, Pyz*1000, Pzz*1000);
}


void APolarSite::PrintTensorPDB(FILE *out, int state) {

    this->Charge(state);

    double scale = 30.;

    vec rxyz = _pos*scale;
    
    double rx = rxyz.getX();
    double ry = rxyz.getY();
    double rz = rxyz.getZ();

    vec end_pax = rxyz + 1.0*this->pax*this->eigenpxx/this->eigendamp;
    vec end_pay = rxyz + 1.0*this->pay*this->eigenpyy/this->eigendamp;
    vec end_paz = rxyz + 1.0*this->paz*this->eigenpzz/this->eigendamp;

    fprintf(out, "%-5s %+4.7f %+4.7f %+4.7f\n",
            this->_name.c_str(),
            rx,
            ry,
            rz);
    fprintf(out, "X     %+4.7f %+4.7f %+4.7f\n",
            end_pax.getX(),
            end_pax.getY(),
            end_pax.getZ());
    fprintf(out, "Y     %+4.7f %+4.7f %+4.7f\n",
            end_pay.getX(),
            end_pay.getY(),
            end_pay.getZ());
    fprintf(out, "Z     %+4.7f %+4.7f %+4.7f\n",
            end_paz.getX(),
            end_paz.getY(),
            end_paz.getZ());

}


void APolarSite::WritePdbLine(FILE *out, const string &tag) {
    
    fprintf(out, "ATOM  %5d %4s%1s%3s %1s%4d%1s   "
              "%8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s%4.7f\n",
         _id % 100000,          // Atom serial number           %5d
         _name.c_str(),         // Atom name                    %4s
         " ",                   // alternate location indicator.%1s
         tag.c_str(),           // Residue name.                %3s
         "A",                   // Chain identifier             %1s
         _id % 10000,           // Residue sequence number      %4d
         " ",                   // Insertion of residues.       %1s
         _pos.getX()*10,        // X in Angstroms               %8.3f
         _pos.getY()*10,        // Y in Angstroms               %8.3f
         _pos.getZ()*10,        // Z in Angstroms               %8.3f
         1.0,                   // Occupancy                    %6.2f
         0.0,                   // Temperature factor           %6.2f
         " ",                   // Segment identifier           %4s
         _name.c_str(),          // Element symbol               %2s
         " ",                   // Charge on the atom.          %2s
         Q00
         );    
}


void APolarSite::WriteXyzLine(FILE *out, vec &shift, string format) {

    double int2ext = 1.0;

    vec pos = _pos + shift;

    if (format == "gaussian") {
        int2ext = 10.;
    }
    else {
        int2ext = 10.;
    }

    fprintf(out, "%-2s %+4.9f %+4.9f %+4.9f \n",
            _name.c_str(),
            pos.getX()*10, pos.getY()*10, pos.getZ()*10);
}


void APolarSite::WriteChkLine(FILE *out, vec &shift, bool split_dpl,
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
            //cout << endl
            //     << "WARNING: Quadrupoles are not split onto point charges."
            //     << endl;

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
        
        if (this->eigendamp == 0) {
            A = pos;
            B = pos;
            qA = 0;
            qB = 0;
        }

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


vector<APolarSite*> APS_FROM_MPS(string filename, int state) {

    int poleCount = 1;
    double Q0_total = 0.0;
    string units = "";
    bool useDefaultPs = true;

    vector<APolarSite*> poles;
    APolarSite *thisPole = NULL;

    vector<double> Qs; // <- multipole moments
    matrix         P1; // <- dipole polarizability

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);
            
            if ( !split.size()      ||
                  split[0] == "!"   ||
                  split[0].substr(0,1) == "!" ) { continue; }

    // ! Interesting information here, e.g.
    // ! DCV2T opt
    // ! SP        RB3LYP          6-311+G(d,p)
    // Units bohr
    //
    // C          -4.2414603400   -3.8124751600    0.0017575736    Rank  2
    //  -0.3853409355
    //  -0.0002321905   0.2401559510   0.6602334308
    //  -0.7220625314   0.0004894995  -0.0003833545   0.4526409813  -0.50937399
    //  P 1.75


            // Units used
            if ( split[0] == "Units") {
                units = split[1];
                if (units != "bohr" && units != "angstrom") {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }
            }

            // element,  position,  rank limit
            else if ( split.size() == 6 ) {

                Qs.clear();

                int id = poleCount++;  // <- starts from 1
                string name = split[0];

                double BOHR2NM = 0.0529189379;
                double ANGSTROM2NM = 0.1;
                double x, y, z;

                if (units == "bohr") {
                    x = BOHR2NM * boost::lexical_cast<double>(split[1]);
                    y = BOHR2NM * boost::lexical_cast<double>(split[2]);
                    z = BOHR2NM * boost::lexical_cast<double>(split[3]);
                }
                else if (units == "angstrom") {
                    x = ANGSTROM2NM * boost::lexical_cast<double>(split[1]);
                    y = ANGSTROM2NM * boost::lexical_cast<double>(split[2]);
                    z = ANGSTROM2NM * boost::lexical_cast<double>(split[3]);
                }
                else {
                    throw std::runtime_error( "Unit " + units + " in file "
                                            + filename + " not supported.");
                }

                vec pos = vec(x,y,z);

                int rank = boost::lexical_cast<int>(split[5]);

                APolarSite *newPole = new APolarSite(id, name);
                newPole->setRank(rank);
                newPole->setPos(pos);
                poles.push_back(newPole);
                thisPole = newPole;

            }

            // 'P', dipole polarizability
            else if ( split[0] == "P") {

                double pxx, pxy, pxz;
                double      pyy, pyz;
                double           pzz;

                if (split.size() == 7) {

                    pxx = 1e-3 * boost::lexical_cast<double>(split[1]);
                    pxy = 1e-3 * boost::lexical_cast<double>(split[2]);
                    pxz = 1e-3 * boost::lexical_cast<double>(split[3]);
                    pyy = 1e-3 * boost::lexical_cast<double>(split[4]);
                    pyz = 1e-3 * boost::lexical_cast<double>(split[5]);
                    pzz = 1e-3 * boost::lexical_cast<double>(split[6]);

                    P1 = matrix(vec(pxx,pxy,pyy),
                                vec(pxy,pyy,pyz),
                                vec(pxz,pyz,pzz));


                }

                else if (split.size() == 2) {

                    pxx = 1e-3 * boost::lexical_cast<double>(split[1]);
                    pxy = 0.0;
                    pxz = 0.0;
                    pyy = pxx;
                    pyz = 0.0;
                    pzz = pxx;

                    P1 = matrix(vec(pxx,pxy,pxz),
                                vec(pxy,pyy,pyz),
                                vec(pxz,pyz,pzz));
                }

                else {
                    throw std::runtime_error("Invalid line in " + filename
                                             + ": " + line);
                }

                thisPole->setPs(P1, state);
                useDefaultPs = false;
            }

            // multipole line
            else {

                int lineRank = int( sqrt(thisPole->getQs(state).size()) + 0.5 );

                if (lineRank == 0) {
                    Q0_total += boost::lexical_cast<double>(split[0]);
                }

                for (int i = 0; i < split.size(); i++) {

                    double qXYZ = boost::lexical_cast<double>(split[i]);

                    // Convert e*(a_0)^k to e*(nm)^k where k = rank
                    double BOHR2NM = 0.0529189379;
                    qXYZ *= pow(BOHR2NM, lineRank); // OVERRIDE

                    Qs.push_back(qXYZ);

                }
                thisPole->setQs(Qs, state);
            }
            
        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << filename << endl;
           throw runtime_error("Please supply input file.");           }


    printf("\n... ... ... Reading %-25s -> N = %2d Q0(Sum) = %+1.3f ",
                          filename.c_str(), poles.size(),  Q0_total);

    if (useDefaultPs) {

        cout << endl << "... ... ... NOTE Using default Thole polarizabilities "
             << "for charge state " << state << ". ";

        vector< APolarSite* > ::iterator pol;
        for (pol = poles.begin(); pol < poles.end(); ++pol) {
            string elem = (*pol)->getName();
            double alpha = 0.0;
            // Original set of Thole polarizabilites
            if      (elem == "C") { alpha = 1.75e-3;  } // <- conversion from
            else if (elem == "H") { alpha = 0.696e-3; } //    A³ to nm³ = 10⁻³
            else if (elem == "N") { alpha = 1.073e-3; }
            else if (elem == "O") { alpha = 0.837e-3; }
            else if (elem == "S") { alpha = 2.926e-3; }
            // Different set of Thole polarizabilities
            //if      (elem == "C") { alpha = 1.334e-3; } // <- conversion from
            //else if (elem == "H") { alpha = 0.496e-3; } //    A³ to nm³ = 10⁻³
            //else if (elem == "N") { alpha = 1.073e-3; }
            //else if (elem == "O") { alpha = 0.837e-3; }
            //else if (elem == "S") { alpha = 3.300e-3; }
            else { throw runtime_error("No polarizability given "
                                       "for polar site type " + elem + ". "); }


            P1 = matrix(vec(alpha,0,0),vec(0,alpha,0),vec(0,0,alpha));

            (*pol)->setPs(P1, state);
        }
    }

    vector< APolarSite* > ::iterator pol;
    for (pol = poles.begin(); pol < poles.end(); ++pol) {
        (*pol)->Charge(state);
    }

    return poles;
}



//inline void BasicInteractor::FieldInduAlpha(APolarSite &pol1,
//                                            APolarSite &pol2) {
//
//    // NOTE >>> Exp. damping factor a must have been set     <<< NOTE //
//    // NOTE >>> e12 points from polar site 1 to polar site 2 <<< NOTE //
//    e12  = pol2.getPos() - pol1.getPos();
//    R    = 1/abs(e12);
//    R2   = R*R;
//    R3   = R2*R;
//    R4   = R3*R;
//    R5   = R4*R;
//    e12 *= R;
//
//    // Thole damping init.
//    u3   = 1 / (R3 * sqrt(pol1.getIsoP() * pol2.getIsoP()));
//
////        rax =   pol1._locX * e12;
////        ray =   pol1._locY * e12;
////        raz =   pol1._locZ * e12;
////        rbx = - pol2._locX * e12;
////        rby = - pol2._locY * e12;
////        rbz = - pol2._locZ * e12;
//
//        rax = e12.getX();
//        ray = e12.getY();
//        raz = e12.getZ();
//        rbx = - rax;
//        rby = - ray;
//        rbz = - raz;
//
////        cxx = pol1._locX * pol2._locX;
////        cxy = pol1._locX * pol2._locY;
////        cxz = pol1._locX * pol2._locZ;
////        cyx = pol1._locY * pol2._locX;
////        cyy = pol1._locY * pol2._locY;
////        cyz = pol1._locY * pol2._locZ;
////        czx = pol1._locZ * pol2._locX;
////        czy = pol1._locZ * pol2._locY;
////        czz = pol1._locZ * pol2._locZ;
//
//        cxx = 1;
//        cxy = 0;
//        cxz = 0;
//        cyx = 0;
//        cyy = 1;
//        cyz = 0;
//        czx = 0;
//        czy = 0;
//        czz = 1;
//
//    // Fields generated by rank-1 induced m'poles
//
//    if (a*u3 < 40.0) {
//        pol1.FUx += TU1x_1x() * pol2.U1x;
//        pol1.FUx += TU1x_1y() * pol2.U1y;
//        pol1.FUx += TU1x_1z() * pol2.U1z;
//        pol1.FUy += TU1y_1x() * pol2.U1x;
//        pol1.FUy += TU1y_1y() * pol2.U1y;
//        pol1.FUy += TU1y_1z() * pol2.U1z;
//        pol1.FUz += TU1z_1x() * pol2.U1x;
//        pol1.FUz += TU1z_1y() * pol2.U1y;
//        pol1.FUz += TU1z_1z() * pol2.U1z;
//
//        pol2.FUx += TU1x_1x() * pol1.U1x;
//        pol2.FUx += TU1y_1x() * pol1.U1y;
//        pol2.FUx += TU1z_1x() * pol1.U1z;
//        pol2.FUy += TU1x_1y() * pol1.U1x;
//        pol2.FUy += TU1y_1y() * pol1.U1y;
//        pol2.FUy += TU1z_1y() * pol1.U1z;
//        pol2.FUz += TU1x_1z() * pol1.U1x;
//        pol2.FUz += TU1y_1z() * pol1.U1y;
//        pol2.FUz += TU1z_1z() * pol1.U1z;
//    }
//    else {
//        pol1.FUx += T1x_1x() * pol2.U1x;
//        pol1.FUx += T1x_1y() * pol2.U1y;
//        pol1.FUx += T1x_1z() * pol2.U1z;
//        pol1.FUy += T1y_1x() * pol2.U1x;
//        pol1.FUy += T1y_1y() * pol2.U1y;
//        pol1.FUy += T1y_1z() * pol2.U1z;
//        pol1.FUz += T1z_1x() * pol2.U1x;
//        pol1.FUz += T1z_1y() * pol2.U1y;
//        pol1.FUz += T1z_1z() * pol2.U1z;
//
//        pol2.FUx += T1x_1x() * pol1.U1x;
//        pol2.FUx += T1y_1x() * pol1.U1y;
//        pol2.FUx += T1z_1x() * pol1.U1z;
//        pol2.FUy += T1x_1y() * pol1.U1x;
//        pol2.FUy += T1y_1y() * pol1.U1y;
//        pol2.FUy += T1z_1y() * pol1.U1z;
//        pol2.FUz += T1x_1z() * pol1.U1x;
//        pol2.FUz += T1y_1z() * pol1.U1y;
//        pol2.FUz += T1z_1z() * pol1.U1z;
//    }
//}



}}