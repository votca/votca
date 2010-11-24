/*
 * File:   main.cpp
 * Author: ruehle
 *
 * Created on July 6, 2010, 12:15 PM
 */

#include <stdlib.h>
#include <votca/csg/csgapplication.h>
#include <votca/tools/histogramnew.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/nblist.h>
#include <votca/csg/nblistgrid.h>

//using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgTestApp
: public CsgApplication {

    string ProgramName() {
        return "template_nblist";
    }

    void HelpText(ostream &out) {
        out << "rough template for rdf calculations";
    }

    /* why not setting a member variable? */
    bool DoThreaded() {
        return true;
    }
    void Initialize();

    bool DoTrajectory() {
        return true;
    }

    void BeginEvaluate(Topology *top, Topology *top_ref);
    void EndEvaluate();

    CsgApplication::Worker *ForkWorker(void);
    void MergeWorker(Worker *worker);
protected:
    Mutex rdfMutex;
    HistogramNew _rdf;
    double _cut_off;

};

class RDFWorker
: public CsgApplication::Worker {
public:
    /*TODO do we need con/destructor?*/
    //RDFWorker();
    ~RDFWorker();
    void EvalConfiguration(Topology *top, Topology *top_ref);
    // funktionen implementieren, unter anerem initialize
    HistogramNew _rdf;
    double _cut_off;

};

int main(int argc, char** argv) {
    CsgTestApp app;

    return app.Exec(argc, argv);
}

void CsgTestApp::Initialize() {
    CsgApplication::Initialize();
    AddProgramOptions("RDF options")
            ("c", boost::program_options::value<double>()->default_value(1.0), "the cutoff");
}

void CsgTestApp::BeginEvaluate(Topology *top, Topology *top_ref) {
    _cut_off = OptionsMap()["c"].as<double>();
    _rdf.Initialize(0, _cut_off, 50);
}

void CsgTestApp::EndEvaluate() {
    _rdf.data().y() = //_avg_vol.getAvg() * iter->second->_norm *
            element_div(_rdf.data().y(),
            element_prod(_rdf.data().x(), _rdf.data().x())
            );

    _rdf.data().Save("rdf.dat");
}

CsgApplication::Worker * CsgTestApp::ForkWorker() {
    //std::cout << "i am so forking a worker right now" << std::endl;
    RDFWorker *worker;
    worker = new RDFWorker();
    // initializiseren here?
    worker->_cut_off = OptionsMap()["c"].as<double>();
    worker->_rdf.Initialize(0, worker->_cut_off, 50);
    return worker; // evtl cast
    //i->_norm = 1. / (4. * M_PI * i->_step * beads1.size()*(beads2.size() - 1.) / 2.);
}

void CsgTestApp::MergeWorker(Worker *worker) {
    RDFWorker * myRDFWorker;
    myRDFWorker = dynamic_cast<RDFWorker*> (worker);
    rdfMutex.Lock();
    _rdf.data().y() = _rdf.data().y() + myRDFWorker->_rdf.data().y();
    //     for (int i=0; i<myRDFWorker->_rdf.data().size(); i++)
    //         std::cout << myRDFWorker->_rdf.data().x(i) << std::endl;

    rdfMutex.Unlock();
}

RDFWorker::~RDFWorker(void) {
}

void RDFWorker::EvalConfiguration(Topology *top, Topology *top_ref) {

    //std::cout << "I am so hard working! worker id: " << getId() << ", frame nr: " << top->getStep() << std::endl;
    //         << ", topology copy addr: " << top << std::endl;
    //    sleep(rand()%3);
    BeadList b;
    b.Generate(*top, "*");
    NBListGrid nb;
    nb.setCutoff(_cut_off);
    nb.Generate(b);
    NBList::iterator i;
    for (i = nb.begin(); i != nb.end(); ++i) {
        _rdf.Process((*i)->dist());
    }
    //std::cout << "      " << _rdf. << std::endl;
    //    _rdf.y() = _rdf.y() + worker.rdf.data().y() // so ungefaehr
}