#ifndef VOTCA_XTP_CGPOLAR_H
#define VOTCA_XTP_CGPOLAR_H

#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/xmapper.h>
#include <votca/xtp/dmaspace.h>
#include <votca/xtp/molpolengine.h>
#include <boost/filesystem.hpp>

namespace votca { 
namespace xtp {
  
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
    bool                           _mps_check;
    string                         _load_ptop_archfile;
    
    bool                           _cg_anisotropic;

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
        if (opt->exists(key+".mps_check")) {
            _mps_check = opt->get(key+".mps_check").as<bool>();
        }
        else { _mps_check = false; }
        if (opt->exists(key+".load_ptop_from")) {
            _load_ptop_archfile = opt->get(key+".load_ptop_from").as<string>();
        }
        else { _load_ptop_archfile = ""; }
    
    key = "options.cgpolar.coarsegrain";
        if (opt->exists(key+".cg_anisotropic")) {
            _cg_anisotropic = opt->get(key+".cg_anisotropic").as<bool>();
        }
        else { _cg_anisotropic = false; }

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
    
    if (_mps_check) {
        string atomistic = "mps_mapped";
        if (!boost::filesystem::exists(atomistic)) {
            boost::filesystem::create_directory(atomistic);
        }
        for (vector<PolarSeg*>::iterator sit = bgn.begin();
            sit < bgn.end(); ++sit) {
            string mpsfile = atomistic + "/" + boost::lexical_cast<string>((*sit)->getId()) + ".mps";
            (*sit)->WriteMPS(mpsfile, "<cgpolar> (atomistic)");
        }
    }
    
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
    
    MolPolEngine engine = MolPolEngine(0.39, 0.30, 1024, 0.0001);
    
    for (vector<PolarSeg*>::iterator sit=bgn.begin();
        sit<bgn.end(); ++sit) {
        //matrix p0 = engine.CalculateMolPol(*(*sit), true);
        engine.CalculateMolPol(*(*sit), true);
    }
    
    // COARSE-GRAIN    
    cout << endl;
    for (vector<PolarSeg*>::iterator sit=bgn.begin();
        sit<bgn.end(); ++sit) {
        //MolPolEngine engine = MolPolEngine();
        //matrix p0 = engine.CalculateMolPol(*(*sit), true);
        engine.CalculateMolPol(*(*sit), true);
        //(*sit)->WriteMPS("cgpolar.fine.mps", "FINE");
        cout << "\rCoarse-grain ID = " << (*sit)->getId() << flush;
        (*sit)->Coarsegrain(_cg_anisotropic);        
        //(*sit)->WriteMPS("cgpolar.coarse.mps", "COARSE");
        //matrix p1 = engine.CalculateMolPol(*(*sit), true);
        //int a; cin >> a;
    }
    
    if (_mps_check) {
        string coarse = "mps_coarse";
        if (!boost::filesystem::exists(coarse)) {
            boost::filesystem::create_directory(coarse);
        }
        for (vector<PolarSeg*>::iterator sit = bgn.begin();
            sit < bgn.end(); ++sit) {
            string mpsfile = coarse + "/" + boost::lexical_cast<string>((*sit)->getId()) + ".mps";
            (*sit)->WriteMPS(mpsfile, "<cgpolar> (coarse)");
        }
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