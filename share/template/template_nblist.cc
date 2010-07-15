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
    : public CsgApplication
{
    string ProgramName() { return "csg_dump"; }
    void HelpText(ostream &out) { out << "Print atoms that are read from topology file to help"
        " debugging atom naming."; }

    void Initialize();

    bool DoTrajectory() {return true;}

    void BeginEvaluate(Topology *top, Topology *top_ref);
    void EvalConfiguration(Topology *top, Topology *top_ref);
    void EndEvaluate();

protected:
    HistogramNew _rdf;
    double _cut_off;

};

int main(int argc, char** argv)
{
    CsgTestApp app;

    return app.Exec(argc, argv);
}

void CsgTestApp::EvalConfiguration(Topology *top, Topology *top_ref)
{
    BeadList b;
    b.Generate(*top, "*");
    NBListGrid nb;
    nb.setCutoff(_cut_off);
    nb.Generate(b);
    NBList::iterator i;
    for(i=nb.begin(); i!=nb.end(); ++i)
        _rdf.Process((*i)->dist());
}

void CsgTestApp::Initialize()
{
    CsgApplication::Initialize();
    AddProgramOptions("RDF options")
             ("c", boost::program_options::value<double>()->default_value(1.0), "the cutoff");
}

void CsgTestApp::BeginEvaluate(Topology *top, Topology *top_ref)
{
    _cut_off = OptionsMap()["c"].as<double>();
    _rdf.Initialize(0, _cut_off, 50);
}

void CsgTestApp::EndEvaluate()
{
    _rdf.data().Save("rdf.dat");
}

