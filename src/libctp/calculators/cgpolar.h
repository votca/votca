#ifndef VOTCA_CTP_CGPOLAR_H
#define VOTCA_CTP_CGPOLAR_H

#include <votca/ctp/qmcalculator.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/dmaspace.h>
#include <votca/ctp/molpolengine.h>

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
    string                         _load_ptop_archfile;

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
        if (opt->exists(key+".load_ptop_from")) {
            _load_ptop_archfile = opt->get(key+".load_ptop_from").as<string>();
        }
        else { _load_ptop_archfile = ""; }

    return;
}




bool CgPolar::EvaluateFrame(Topology *top) {    
    
    // GENERATE POLAR TOPOLOGY
    PolarTop ptop(top);    
    if (_load_ptop_archfile == "") {
        cout << endl << "... ... Create polar topology from map" << flush;
        cout << endl << "... ... ... Initialize MPS-mapper: " << flush;
        _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
        _mps_mapper.Gen_BGN(top, &ptop, NULL);        
    }
    else {
        cout << endl << "... ... Load polar topology from hard drive: '"
            << _load_ptop_archfile << "'" << flush;
        ptop.LoadFromDrive(_load_ptop_archfile);
    }
    vector<PolarSeg*> bgn = ptop.BGN();
    
    // VERIFY INPUT: PDB, PTOP, XML
    if (_pdb_check) ptop.PrintPDB("cgpolar.fine.pdb");
    ptop.SaveToDrive("cgpolar.fine.ptop");    
    ofstream ofs;    
    ofs.open("cgpolar.fine.xml", ofstream::out);
    for (vector<PolarSeg*>::iterator sit = bgn.begin();
        sit < bgn.end(); ++sit) {
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit < (*sit)->end(); ++pit) {
            (*pit)->WriteXmlLine(ofs);
        }
    }
    ofs.close();    
    
    // COARSE-GRAIN    
    cout << endl;
    for (vector<PolarSeg*>::iterator sit=bgn.begin();
        sit<bgn.end(); ++sit) {
        //MolPolEngine engine = MolPolEngine();
        //engine.CalculateMolPol(*(*sit), true);
        //(*sit)->WriteMPS("cgpolar.fine.mps", "FINE");
        cout << "\rCoarse-grain ID = " << (*sit)->getId() << flush;
        (*sit)->Coarsegrain();
        //(*sit)->WriteMPS("cgpolar.coarse.mps", "COARSE");
        //engine.CalculateMolPol(*(*sit), true);
    }
    
    // VERIFY OUTPUT: PDB, PTOP, XML
    if (_pdb_check) ptop.PrintPDB("cgpolar.coarse.pdb");
    ptop.SaveToDrive("cgpolar.coarse.ptop");
    ofs.open("cgpolar.coarse.xml", ofstream::out);
    for (vector<PolarSeg*>::iterator sit = bgn.begin();
        sit < bgn.end(); ++sit) {
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit < (*sit)->end(); ++pit) {
            (*pit)->WriteXmlLine(ofs);
        }
    }
    ofs.close();
    
    return true;
}
    
}}

#endif