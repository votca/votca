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

private:

    string      _TI_tag;
};


void IImport::Initialize(Topology *top, Property *options) {

    string key = "options.iimport";
    _TI_tag = options->get(key+".TI_tag").as< string >();

    cout << endl << "... ... Using TI XML tag '" << _TI_tag << "'" << flush;
}

bool IImport::EvaluateFrame(Topology *top) {

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
                    printf("\n... ... ... J2(State = %+1d) = %4.7e",
                            STATE, qmpair->getJeff2(STATE));
                }
                else if (_TI_tag == "J_sq_degen" || _TI_tag == "J_sq_boltz") {
                    double J2 = boost::lexical_cast<double>(split[1]);
                    qmpair->setJeff2(J2, STATE);
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

}

}}

#endif
