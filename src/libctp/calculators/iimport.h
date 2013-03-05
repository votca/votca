#ifndef IIMPORT_H
#define IIMPORT_H

#include <votca/ctp/qmcalculator.h>
#include <sys/stat.h>


namespace votca { namespace ctp {

class IImport : public QMCalculator
{
public:

    string      Identify() { return "IImport"; }

    void        Initialize(Topology *top, Property *options);
    bool        EvaluateFrame(Topology *top);
    void        XML2PairTI(QMPair *qmpair, string &xmlDirFile);
    void        List2PairsTI(Topology *top, string &ti_file);

private:

    bool        _importFromDirs;
    bool        _importFromList;

    string      _TI_tag;
    string      _TI_file;
};


void IImport::Initialize(Topology *top, Property *options) {

    _importFromDirs = false;
    _importFromList = false;

    string key = "options.iimport";

    if (options->exists(key+".TI_tag")) {
        _TI_tag = options->get(key+".TI_tag").as< string >();
        _importFromDirs = true;
        cout << endl << "... ... Using TI XML tag '" << _TI_tag << "'" << flush;
    }

    else if (options->exists(key+".TI_file")) {
        _TI_file = options->get(key+".TI_file").as< string >();
        _importFromList = true;
        cout << endl << "... ... Using TI table '" << _TI_file << "'" << flush;
    }
    

    
}

bool IImport::EvaluateFrame(Topology *top) {

  // Import from file
  if (_importFromList) {
      this->List2PairsTI(top, _TI_file);
  }

  // Import from DiPro folders
  else if (_importFromDirs) {
    string PARDIR = "frame" + boost::lexical_cast<string>(top->getDatabaseId());
    string SUBDIR = PARDIR + "/pair_0000_0000";

    QMNBList &nblist = top->NBList();
    QMNBList ::iterator nit;

    for (nit = nblist.begin(); nit != nblist.end(); ++nit) {
        string ID1  = boost::lexical_cast<string>((*nit)->Seg1()->getId());
        string ID2  = boost::lexical_cast<string>((*nit)->Seg2()->getId());

        SUBDIR = PARDIR + "/pair_" + ID1 + "_" + ID2 + "/TI.xml";
        this->XML2PairTI(*nit, SUBDIR);
    }
  }

}


void IImport::List2PairsTI(Topology *top, string &ti_file) {

    QMNBList &nblist = top->NBList();
    int pair_count = 0;

    std::string line;
    std::ifstream intt;
    intt.open(ti_file.c_str());

    if (intt.is_open()) {
        while (intt.good()) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t");
            toker.ToVector(split);

            if (split.size() == 0) { continue; }
            
            if (split.size() != 4 && split.size() != 6) {
                cout << endl << "... ... Invalid line: " << line << flush;
                continue;
            }

            int seg1id = boost::lexical_cast<int>(split[0]);
            int seg2id = boost::lexical_cast<int>(split[1]);

            Segment *seg1 = top->getSegment(seg1id);
            Segment *seg2 = top->getSegment(seg2id);

            QMPair *qmp = nblist.FindPair(seg1,seg2);
            if (qmp == NULL) {
                cout << endl
                     << "... ... ERROR: " << line
                     << flush;
                cout << endl
                     << "... ... Line is not compatible with neighborlist. "
                     << flush;
                throw std::runtime_error("Forgot to run -e neighborlist?");
            }            

            int e_h_1 = (split[2] == "e") ? -1 : +1;
            double j2_1 = boost::lexical_cast<double>(split[3]);

            if (split.size() == 6) {

                int e_h_2 = (split[4] == "e") ? -1 : +1;
                if (e_h_1 == e_h_2) {
                    cout << endl << "... ... Invalid line: " << line << flush;
                    continue;
                }
                double j2_2 = boost::lexical_cast<double>(split[5]);

                qmp->setJeff2(j2_2, e_h_2);
                qmp->setIsPathCarrier(true, e_h_2);
            }

            qmp->setJeff2(j2_1, e_h_1);
            qmp->setIsPathCarrier(true, e_h_1);
            ++pair_count;
            
        }
    }
    else {
        cout << endl
             << "ERROR: No such file" << ti_file << ". "
             << flush;
        throw std::runtime_error("Missing input file in IImport calculator.");
    }

    cout << endl
         << "... ... Set transfer integrals for " << pair_count << " pairs. "
         << flush;
}


void IImport::XML2PairTI(QMPair *qmpair, string &xmlDirFile) {

    printf("\n... ... Import TIs for pair %5d (ID1 %4d    ID2 %4d) ",
            qmpair->getId(), qmpair->Seg1()->getId(), qmpair->Seg2()->getId());

    string TRANSPORT = "electron_or_hole";
    int    STATE     = 0;
    
    std::string line;
    std::ifstream intt;
    intt.open(xmlDirFile.c_str());

    if (intt.is_open() ) {
        while ( intt.good() ) {

            std::getline(intt, line);
            vector<string> split;
            Tokenizer toker(line, " \t<>=\"");
            toker.ToVector(split);

            if (split.size() < 3) {
                continue;
            }

//            cout << split.size() << " " << flush;
//            vector<string> ::iterator sit;
//            for (sit = split.begin(); sit < split.end(); ++sit) {
//                cout << (*sit) << " " << flush;
//            }

            if      (split[0] == "transport") {
                assert(split[1] == "name");
                TRANSPORT = split[2];
            }
            else if (split[0] == _TI_tag) {
                assert(TRANSPORT == "electron" || TRANSPORT == "hole");
                STATE = (TRANSPORT == "electron") ? -1 : +1;
                if (_TI_tag == "J" || _TI_tag == "T_00") {
                    double J = boost::lexical_cast<double>(split[1]);
                    qmpair->setJeff2(J*J, STATE);
                    qmpair->setIsPathCarrier(1, STATE);
                    printf("\n... ... ... J2(State = %+1d) = %4.7e",
                            STATE, qmpair->getJeff2(STATE));
                }
                else if (_TI_tag == "J_sq_degen" || _TI_tag == "J_sq_boltz") {
                    double J2 = boost::lexical_cast<double>(split[1]);
                    qmpair->setJeff2(J2, STATE);
                    qmpair->setIsPathCarrier(1, STATE);
                    printf("\n... ... ... J2(State = %+1d) = %4.7e",
                            STATE, qmpair->getJeff2(STATE));
                }
            }
        } /* Exit loop over lines */
    }
    else { cout << endl
                << "... ... ... ERROR: No TI.xml in pair folder. Skip..."
                << flush;
    }

    intt.close();

}

}}

#endif
