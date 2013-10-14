#ifndef _VOTCA_CTP_EIMPORT_H
#define _VOTCA_CTP_EIMPORT_H

#include <votca/ctp/qmcalculator.h>


namespace votca { namespace ctp {

class EImport : public QMCalculator
{
public:

    EImport() {};
   ~EImport() {};

   string Identify() { return "eimport"; }

   void Initialize(Property *options);
   bool EvaluateFrame(Topology *top);

private:

    string _energiesFile;
    bool   _reset;

};


void EImport::Initialize(Property *opt) {

    // update options with the VOTCASHARE defaults   
    UpdateWithDefaults( opt );
    string key = "options." + Identify();

    _energiesFile = opt->get(key + ".energies").as< string >();    

    if (opt->exists(key+".reset")) {
        int reset = opt->get(key+".reset").as<int>();
        _reset = (reset == 1) ? true : false;
        if (_reset) {
            cout << endl
                 << "... ... NOTE: Resetting site energies to zero."
                 << flush;
        }

    }
    else {
        _reset = false;
    }
}

bool EImport::EvaluateFrame(Topology *top) {

    if (_reset == true) {

        vector<Segment*> ::iterator sit;
        for (sit = top->Segments().begin();
             sit < top->Segments().end();
             ++sit) {

            (*sit)->setEMpoles( 0, 0.0);
            (*sit)->setEMpoles(-1, 0.0);
            (*sit)->setEMpoles(+1, 0.0);
        }
    }

    else {

        string filename = _energiesFile;
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


              // Sample line
              // 10 Alq3 0 0.00000000  +1 -0.49966334  0 20  +1 10  SPH 2366 ...

                int    id   = boost::lexical_cast<int>(split[0]);
                string name = boost::lexical_cast<string>(split[1]);

                Segment *seg = top->getSegment(id);

                if (seg->getName() != name) {
                    cout << endl
                         << "... ... ERROR: Structure of input file does "
                            "not match state file. "
                         << "SEG ID " << id
                         << flush;
                    continue;
                }

                int    tell = boost::lexical_cast<int>(split[6]);
                if (tell*tell == 0) {

                    int state_N = boost::lexical_cast<int>(split[2]);
                    int state_C = boost::lexical_cast<int>(split[4]);

                    assert(state_N == 0);
                    assert(state_C*state_C == 1);

                    double e_N = boost::lexical_cast<double>(split[3]);
                    double e_C = boost::lexical_cast<double>(split[5]);

                    seg->setEMpoles(state_N, e_N);
                    seg->setEMpoles(state_C, e_C);
                }

                else if (tell*tell == 1) {

                    int state_N = boost::lexical_cast<int>(split[2]);
                    int state_A = boost::lexical_cast<int>(split[4]);
                    int state_C = boost::lexical_cast<int>(split[6]);

                    assert(state_N == 0);
                    assert(state_A == -1);
                    assert(state_C == +1);

                    double e_N = boost::lexical_cast<double>(split[3]);
                    double e_A = boost::lexical_cast<double>(split[5]);
                    double e_C = boost::lexical_cast<double>(split[7]);

                    seg->setEMpoles(state_N, e_N);
                    seg->setEMpoles(state_A, e_A);
                    seg->setEMpoles(state_C, e_C);
                }

                else {
                    cout << endl << "... ... ERROR: Wrong input format. "
                         << "SEG ID " << id << flush;
                    continue;
                }
            }
        }
        else { cout << endl << "ERROR: No such file " << filename << endl; }
    }

    return 1;

}




}}


#endif
