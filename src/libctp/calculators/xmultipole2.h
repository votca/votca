#ifndef XMULTIPOLE2_H
#define XMULTIPOLE2_H


#include <votca/ctp/qmcalculator.h>

namespace votca { namespace ctp {

class XMP : public QMCalculator
{

public:

    XMP() {};
   ~XMP() {};

    string              Identify() { return "XMultipole"; }

    void                Initialize(Topology *, Property *);

    void                Collect_JOB(string job_file, Topology *top);
    void                Collect_EMP(string emp_file, Topology *top);
    void                Collect_XMP(string xmp_file, Topology *top);
    void                Collect_MPS(Topology *top);
    vector<PolarSite*>  Parse_GDMA(string mps_file, int state);




    bool                EvaluateFrame(Topology *top);


private:

    string                         _job_file;
    string                         _emp_file;

    // Job info : JOB_ID PAIR_ID MPS_1 MPS_2 TAG
    map<int,string>                _jobId_jobTag;
    map<int,int>                   _jobId_pairId;
    map<int,string>                _jobId_mpsFile1;
    map<int,string>                _jobId_mpsFile2;
    map<int, pair<int,int> >       _pairId_seg1Id_seg2Id;

    // Allocate charges to segments (neutral, electron, hole)
    map<int,string>                 _segId_mpsFile_n;   // <=> charge state  0
    map<int,string>                 _segId_mpsFile_e;   // <=> charge state -1
    map<int,string>                 _segId_mpsFile_h;   // <=> charge state +1

    // Allocate charges to pairs (Frenkel, CT)
    map<int,string>                 _pairId_mpsFile_f1;  // <=> charge state -1
    map<int,string>                 _pairId_mpsFile_ct1; // <=> charge state +1
    map<int,string>                 _pairId_mpsFile_f2;  // <=> charge state -1
    map<int,string>                 _pairId_mpsFile_ct2; // <=> charge state +1
    

    // Store polar site container for each GDMA (i.e. .mps) file
    map<string,vector<PolarSite*> > _mpsFile_pSites;

};

void XMP::Initialize(Topology *top, Property *options) {

    string key = "options.xmultipole.control";

    _job_file = options->get(key+".job_file").as<string>();
    _emp_file = options->get(key+".emp_file").as<string>();
}


void XMP::Collect_JOB(string job_file, Topology *top) {

    QMNBList &nblist = top->NBList();

    std::string line;
    std::ifstream intt;
    intt.open(job_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

// Sample line
// # JOB_ID TAG  PAIR_ID SEG1_ID SEG1_NAME SEG1_MPS SEG2_ID SEG2_NAME SEG2_MPS
//   1      E_CT 3819    182     C60       c60.mps  392     DCV       dcv.mps

            int jobId       = boost::lexical_cast<int>(split[0]);
            string tag      = split[1];
            int pairId      = boost::lexical_cast<int>(split[2]);

            int seg1Id      = boost::lexical_cast<int>(split[3]);
            string seg1Name = split[4];
            string seg1mps  = split[5];

            int seg2Id      = boost::lexical_cast<int>(split[6]);
            string seg2Name = split[7];
            string seg2mps  = split[8];

            Segment *seg1   = top->getSegment(seg1Id);
            Segment *seg2   = top->getSegment(seg2Id);
            QMPair  *qmp    = nblist.FindPair(seg1,seg2);

            if (qmp == NULL) {
                cout << endl << "ERROR: '" << job_file
                     << "': No such pair "
                     << pairId << " " << seg1Name << " " << seg2Name << ". "
                     << flush;
                throw runtime_error("Pair specs do not match topology.");
            }

            if ( _jobId_jobTag.count(jobId) ) {
                cout << endl << "ERROR: '" << job_file
                     << "': Job-ID " << jobId << " exists more "
                        "than once. Abort." << endl;
                throw runtime_error("Rework job file.");
            }

            _jobId_jobTag[jobId]            = tag;
            _jobId_pairId[jobId]            = pairId;
            _jobId_mpsFile1[jobId]          = seg1mps;
            _jobId_mpsFile2[jobId]          = seg2mps;
            _pairId_seg1Id_seg2Id[pairId]   = pair<int,int>(seg1Id,seg2Id);


        } /* Exit loop over lines */
    }
    else { cout << endl << "ERROR: No such file " << job_file << endl;
           throw runtime_error("Please supply input file.");           }

    cout << endl << "... ... Registered " << _jobId_jobTag.size() << " jobs. "
         << flush;

}


void XMP::Collect_MPS(Topology *top) {

    map<int,string> ::iterator misit;

    // Neutral species
    for (misit = _segId_mpsFile_n.begin();
         misit != _segId_mpsFile_n.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, 0); }
    }

    // Anion
    for (misit = _segId_mpsFile_e.begin();
         misit != _segId_mpsFile_e.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, -1); }
    }

    // Cation
    for (misit = _segId_mpsFile_h.begin();
         misit != _segId_mpsFile_h.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    // Job Seg1
    for (misit = _jobId_mpsFile1.begin();
         misit != _jobId_mpsFile1.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    // Job Seg2
    for (misit = _jobId_mpsFile2.begin();
         misit != _jobId_mpsFile2.end();
         ++misit) {
        if (_mpsFile_pSites.count(misit->second) > 0 ) { continue; }
        else { _mpsFile_pSites[misit->second] = Parse_GDMA(misit->second, +1); }
    }

    cout << endl << "... ... Parsed a total of " << _mpsFile_pSites.size()
            << " multipole-definition (.mps) files. " << flush;

}


/*
void XMP::Collect_XMP(string xmp_file, Topology *top) {

    QMNBList &nblist = top->NBList();

    std::string line;
    std::ifstream intt;
    intt.open(xmp_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }


            // Sample line:
            // 3472  444 DCV5T 880 C60 f1.mps ct1.mps f2.mps ct2.mps

            assert( split.size() == 8 );                                        // <- File format correct?

            int pairId = boost::lexical_cast<int>(split[0]);
            int seg1Id = boost::lexical_cast<int>(split[1]);
            int seg2Id = boost::lexical_cast<int>(split[3]);

            int seg1Name = split[2];
            int seg2Name = split[4];

            Segment *seg1 = top->getSegment(seg1Id);
            Segment *seg2 = top->getSegment(seg2Id);
            QMPair  *qmp = nblist.FindPair(seg1,seg2);

            if (qmp == NULL) {
                cout << endl << "ERROR: '" << xmp_file
                     << "': No such pair "
                     << pairId << " " << seg1Name << " " << seg2Name << ". "
                     << flush;
                throw runtime_error("Pair specs do not match topology.");
            }

            _pairId_mpsFile_f1[pairId] = split[5];
            _pairId_mpsFile_ct1[pairId] = split[6];
            _pairId_mpsFile_f2[pairId] = split[7];
            _pairId_mpsFile_ct2[pairId] = split[8];
            _pairId_seg1Id_seg2Id[pairId] = pair<int,int>(seg1Id,seg2Id);


        } // Exit loop over lines
    }
    else { cout << endl << "ERROR: No such file " << xmp_file << endl; }
}
*/


void XMP::Collect_EMP(string emp_file, Topology *top) {

    std::string line;
    std::ifstream intt;
    intt.open(emp_file.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
            toker.ToVector(split);

            if ( !split.size()      ||
                  split[0] == "#"   ||
                  split[0].substr(0,1) == "#" ) { continue; }

            // Line sample:
            // 1022 C60 C60_n.mps C60_e.mps C60_h.mps

            int segId       = boost::lexical_cast<int>(split[0]);
            string segName  = boost::lexical_cast<string>(split[1]);

            // Some compliance checks
            assert( split.size() == 5 );                                        // <- File format correct?
            assert( top->getSegment(segId)->getName() == segName );             // <- Input matches topology?

            _segId_mpsFile_n[segId] = split[2];
            _segId_mpsFile_e[segId] = split[3];
            _segId_mpsFile_h[segId] = split[4];

        } // Exit loop over lines
    }
    else { cout << endl << "ERROR: No such file " << emp_file << endl; }

    assert( _segId_mpsFile_n.size() == top->Segments().size() );                // <- Input for all segments?

}


vector<PolarSite*> XMP::Parse_GDMA(string filename, int state) {

    int poleCount = 1;
    double Q0_total = 0.0;
    string units = "";
    bool useDefaultPs = true;

    vector<PolarSite*> poles;
    PolarSite *thisPole = NULL;

    vector<double> Qs; // <- multipole moments
    double         P1; // <- dipole polarizability

    std::string line;
    std::ifstream intt;
    intt.open(filename.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " ");
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
                P1 = -1.;

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

                PolarSite *newPole = new PolarSite(id, name);
                newPole->setRank(rank);
                newPole->setPos(pos);
                poles.push_back(newPole);
                thisPole = newPole;

            }

            // 'P', dipole polarizability
            else if ( split[0] == "P" && split.size() == 2 ) {
                P1 = 1e-3 * boost::lexical_cast<double>(split[1]);
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

    cout << endl << "... ... Reading " << filename <<
                    ": Q0(Total) = " << Q0_total << flush;

    if (useDefaultPs) {

        cout << endl << "... ... ... NOTE Using default Thole polarizabilities "
             << "for charge state " << state << ". ";

        vector< PolarSite* > ::iterator pol;
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
            (*pol)->setPs(alpha, state);
        }
    }

    return poles;
}


bool XMP::EvaluateFrame(Topology *top) {

    Collect_JOB(_job_file, top);
    Collect_EMP(_emp_file, top);
    Collect_MPS(top);

    

}


}}

#endif

