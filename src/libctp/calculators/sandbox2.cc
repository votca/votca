#include "sandbox2.h"

namespace votca { namespace ctp {

    void Sandbox2::Initialize(Topology *top, Property *opt) {

        cout << "Initialize (Sandbox2)..." << endl;

        string key;

        key = "options.sandbox.id";
        if (opt->exists(key)) {
            _ID = opt->get(key).as< int >();
        }

        key = "options.sandbox.sec1";
        if (opt->exists(key+".p1")) {
            _p1 = opt->get(key+".p1").as< double >();
        }

        key = "options.sandbox.sec2";
        if (opt->exists(key+".p2")) {
            _p2 = opt->get(key+".p2").as< double >();
        }

        cout << "P1     " << _p1 << endl;
        cout << "P2     " << _p2 << endl;
        cout << "ID     " << _ID << endl;

    }

    bool Sandbox2::EvaluateFrame(Topology *top) {

        cout << "Calculate (Sandbox2)... " << endl;

        // Segment* seg1 = top->getSegment(1);
        // cout << "Segment 1: occ. = " << seg1->getOcc() << endl;
        // seg1->setOcc(0.5);

        top->PrintInfo(cout);

    }



}}
