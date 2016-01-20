#ifndef VOTCA_XTP_VAVERAGE_H
#define VOTCA_XTP_VAVERAGE_H


#include <votca/xtp/qmcalculator.h>


namespace votca { namespace xtp {

    
class VAverage : public QMCalculator
{

public:

    VAverage() {};
   ~VAverage() {};

    string Identify() { return "vaverage"; }

    void Initialize(Property *opt);
    bool EvaluateFrame(Topology *top);

private:

    vector<int> _carriers;
    string _tabulate_type;
};


void VAverage::Initialize(Property *opt) {
 
    string key = "options." + this->Identify();
    if (opt->exists(key+".carriers"))
        _carriers = opt->get(key+".carriers").as< vector<int> >();
    else {
        _carriers.push_back(-1);
        _carriers.push_back(+1);
    }
    if (opt->exists(key+".tabulate")) 
        _tabulate_type = opt->get(key+".tabulate").as<string>();
    else {
        _tabulate_type = "segments";
    }
    
    if (_tabulate_type != "atoms" && _tabulate_type != "segments") {
        cout << endl;
        throw std::runtime_error("Tabulate type '" + _tabulate_type + "' "
            "not recognised (should be 'atoms' or 'segments').");
    }
    
    return;
}

bool VAverage::EvaluateFrame(Topology *top) {

    // Rigidify system (if possible)
    bool isRigid = top->Rigidify();
    if (!isRigid) {
        return 0;
    }
    
    cout << endl << "... ... Computing velocity average for all sites" << flush;
    
    for (unsigned int i = 0; i < _carriers.size(); ++i) {
        int state = _carriers[i];
        if (state == 0) continue;
        cout << endl << "... ... State = " << state << flush;
        
        // Storage for velocity averages, start from velocity zero.
        vector<vec> vavgs;
        vavgs.resize(top->Segments().size());
        vector<vec>::iterator vit;
        for (vit = vavgs.begin(); vit < vavgs.end(); ++vit) {
            (*vit) = vec(0,0,0);
        }

        // Each QM-pair (i,j) contributes to velocity (i) and velocity (j) as ...
        QMNBList &nblist = top->NBList();
        QMNBList::iterator nit;    
        for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
            QMPair *qmpair = *nit;
            Segment *seg1 = qmpair->Seg1();
            Segment *seg2 = qmpair->Seg2();
            double p1 = seg1->getOcc(state);
            double p2 = seg2->getOcc(state);
            double w12 = qmpair->getRate12(state);
            double w21 = qmpair->getRate21(state);
            vec dR12 = qmpair->R(); // points from seg. (1) to seg. (2)            
            // Split current symmetrically among sites (1) and (2)
            vec v_12_21 = 0.5*(p1*w12 - p2*w21) * dR12;
            vavgs[seg1->getId()-1] += v_12_21;
            vavgs[seg2->getId()-1] += v_12_21;
        }
        
        // Output
        string tabfile = (boost::format("vaverage.state_%1$+d.tab") % state).str();
        ofstream ofs;    
        ofs.open(tabfile.c_str(), ofstream::out);
        if (!ofs.is_open()) {
            throw runtime_error("Bad file handle: " + tabfile);
        }
        ofs << "# 1   |   2     |  3 4 5       |   6  7  8         |" << endl;
        ofs << "# ----|---------|--------------|-------------------|" << endl;
        ofs << "# ID  |   NAME  |  X Y Z [nm]  |   VX VY VZ [m/s]  |" << endl;
        vector<Segment*>::iterator sit;
        for (sit = top->Segments().begin(); sit < top->Segments().end(); ++sit) {
            vec v = vavgs[(*sit)->getId()-1]*1e-9;
            vec r = (*sit)->getPos();
            string name = (*sit)->getName();
            int id = (*sit)->getId();
            if (_tabulate_type == "atoms") {
                vector<Atom*>::iterator ait;
                for (ait = (*sit)->Atoms().begin(); ait < (*sit)->Atoms().end(); ++ait) {
                    vec r_atm = (*ait)->getPos();
                    ofs << (boost::format("%1$4d %2$-10s %3$+1.7e "
                        "%4$+1.7e %5$+1.7e %6$+1.7e %7$+1.7e %8$+1.7e")
                        % id % name % r_atm.getX() % r_atm.getY() % r_atm.getZ()
                        % v.getX() % v.getY() % v.getZ()) << endl;
                }
            }
            else if (_tabulate_type == "segments") {
                ofs << (boost::format("%1$4d %2$-10s %3$+1.7e %4$+1.7e "
                    "%5$+1.7e %6$+1.7e %7$+1.7e %8$+1.7e")
                    % id % name % r.getX() % r.getY() % r.getZ() % v.getX() 
                    % v.getY() % v.getZ()) << endl;
            }
        }
        ofs.close();
        
    }
    
    return true;
}

}}


#endif
