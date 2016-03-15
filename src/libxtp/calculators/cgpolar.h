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

    std::string Identify() { return "cgpolar"; }
    void Initialize(Property *options);
    bool EvaluateFrame(Topology *top);

private:

    Property                       *_options;
    std::string                         _mps_table;
    std::string                         _xml_file;
    XMpsMap                        _mps_mapper;
    bool                           _pdb_check;
    bool                           _mps_check;
    std::string                         _load_ptop_archfile;
    
    bool                           _cg_anisotropic;

};



void CgPolar::Initialize(Property *opt) {
    
    _options = opt;
    std::cout << std::endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << std::flush;    

    std::string key = "options.cgpolar";
        if (opt->exists(key+".multipoles")) {
            _xml_file = opt->get(key+".multipoles").as< std::string >();
        }
        else {
            std::cout << std::endl;
            throw std::runtime_error("No multipole mapping file provided");
        }
    key = "options.cgpolar.control";
        if (opt->exists(key+".mps_table")) {
            _mps_table = opt->get(key+".mps_table").as<std::string>();
        }
        else {
            _mps_table = opt->get(key+".emp_file").as<std::string>();
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
            _load_ptop_archfile = opt->get(key+".load_ptop_from").as<std::string>();
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
        std::cout << std::endl << "... ... Create polar topology from map" << std::flush;
        std::cout << std::endl << "... ... ... Initialize MPS-mapper: " << std::flush;
        _mps_mapper.GenerateMap(_xml_file, _mps_table, top);
        _mps_mapper.Gen_BGN(top, &ptop, NULL);        
    }
    else {
        std::cout << std::endl << "... ... Load polar topology from hard drive: '"
            << _load_ptop_archfile << "'" << std::flush;
        ptop.LoadFromDrive(_load_ptop_archfile);
    }
    std::vector<PolarSeg*> bgn = ptop.BGN();
    
    if (_mps_check) {
        std::string atomistic = "mps_mapped";
        if (!boost::filesystem::exists(atomistic)) {
            boost::filesystem::create_directory(atomistic);
        }
        for (std::vector<PolarSeg*>::iterator sit = bgn.begin();
            sit < bgn.end(); ++sit) {
            std::string mpsfile = atomistic + "/" + boost::lexical_cast<std::string>((*sit)->getId()) + ".mps";
            (*sit)->WriteMPS(mpsfile, "<cgpolar> (atomistic)");
        }
    }
    
    // VERIFY INPUT: PDB, PTOP, XML
    if (_pdb_check) ptop.PrintPDB("cgpolar.fine.pdb");
    ptop.SaveToDrive("cgpolar.fine.ptop");    
    std::ofstream ofs;    
    ofs.open("cgpolar.fine.xml", std::ofstream::out);
    for (vector<PolarSeg*>::iterator sit = bgn.begin();
        sit < bgn.end(); ++sit) {
        for (PolarSeg::iterator pit = (*sit)->begin();
            pit < (*sit)->end(); ++pit) {
            (*pit)->WriteXmlLine(ofs);
        }
    }
    ofs.close();    
    
    MolPolEngine engine = MolPolEngine(0.39, 0.30, 1024, 0.0001);
    
    for (std::vector<PolarSeg*>::iterator sit=bgn.begin();
        sit<bgn.end(); ++sit) {
        //matrix p0 = engine.CalculateMolPol(*(*sit), true);
        engine.CalculateMolPol(*(*sit), true);
    }
    
    // COARSE-GRAIN    
    std::cout << std::endl;
    for (std::vector<PolarSeg*>::iterator sit=bgn.begin();
        sit<bgn.end(); ++sit) {
        //MolPolEngine engine = MolPolEngine();
        //matrix p0 = engine.CalculateMolPol(*(*sit), true);
        engine.CalculateMolPol(*(*sit), true);
        //(*sit)->WriteMPS("cgpolar.fine.mps", "FINE");
        std::cout << "\rCoarse-grain ID = " << (*sit)->getId() << std::flush;
        (*sit)->Coarsegrain(_cg_anisotropic);        
        //(*sit)->WriteMPS("cgpolar.coarse.mps", "COARSE");
        //matrix p1 = engine.CalculateMolPol(*(*sit), true);
        //int a; cin >> a;
    }
    
    if (_mps_check) {
        std::string coarse = "mps_coarse";
        if (!boost::filesystem::exists(coarse)) {
            boost::filesystem::create_directory(coarse);
        }
        for (std::vector<PolarSeg*>::iterator sit = bgn.begin();
            sit < bgn.end(); ++sit) {
            std::string mpsfile = coarse + "/" + boost::lexical_cast<std::string>((*sit)->getId()) + ".mps";
            (*sit)->WriteMPS(mpsfile, "<cgpolar> (coarse)");
        }
    }
    
    // VERIFY OUTPUT: PDB, PTOP, XML
    if (_pdb_check) ptop.PrintPDB("cgpolar.coarse.pdb");
    ptop.SaveToDrive("cgpolar.coarse.ptop");
    ofs.open("cgpolar.coarse.xml", std::ofstream::out);
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