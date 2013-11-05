#ifndef VOTCA_CTP_CGPOLAR_H
#define VOTCA_CTP_CGPOLAR_H

#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/dmaspace.h>

namespace votca { 
namespace ctp {
  
class CgPolar : public QMCalculator
{
public:

    CgPolar() {};
   ~CgPolar() {};

    string Identify() { return "cgpolar"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

    Property                       *_options;
    string                         _mps_table;
    string                         _xml_file;
    XMpsMap                        _mps_mapper;
    bool                           _pdb_check;

};



void CgPolar::Initialize(Property *opt) {
    
    _options = opt;
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;    

    string key = "options.cgpolar";
        if (opt->exists(key+".multipoles")) {
            _xml_file = opt->get(key+".multipoles").as< string >();
        }
        else {
            cout << endl;
            throw std::runtime_error("No multipole mapping file provided");
        }
    key = "options.cgpolar.control";
        if (opt->exists(key+".mps_table")) {
            _mps_table = opt->get(key+".mps_table").as<string>();
        }
        else {
            _mps_table = opt->get(key+".emp_file").as<string>();
        }
        if (opt->exists(key+".pdb_check")) {
            _pdb_check = opt->get(key+".pdb_check").as<bool>();
        }
        else { _pdb_check = false; }

    return;
}




bool CgPolar::EvaluateFrame(Topology *top) {    
    
    // GENERATE BACKGROUND (= periodic bg, with empty foreground)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    PolarTop ptop(top);
    _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
    _mps_mapper.Gen_BGN(top, &ptop, NULL);
    if (_pdb_check) ptop.PrintPDB("cgpolar.fine.pdb");
    
    // COARSE-GRAIN
    vector<PolarSeg*> bgn = ptop.BGN();
    for (vector<PolarSeg*>::iterator sit=bgn.begin()+288;
        sit<bgn.end(); ++sit) {
        (*sit)->Coarsegrain();
        break;
    }
    if (_pdb_check) ptop.PrintPDB("cgpolar.coarse.pdb");
    
    
    return true;
}
    
}}

#endif