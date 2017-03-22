/* 
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _VOTCA_XTP_DFT_H
#define _VOTCA_XTP_DFT_H

#include <stdio.h>

#include <votca/ctp/atom.h>
#include <votca/ctp/segment.h>

#include <votca/ctp/logger.h>
#include <votca/xtp/dftengine.h>
#include <votca/xtp/qmpackagefactory.h>


namespace votca { namespace xtp {
    using namespace std;
    
class DFT : public ctp::QMTool
{
public:

    DFT() { };
   ~DFT() { };

    string Identify() { return "dft"; }

    void   Initialize(Property *options);
    bool   Evaluate();

    DFTENGINE _dftengine;
 

private:
    
    string      _orbfile;
    string      _xyzfile;

    string      _logfile;

    string      _package;
    Property    _package_options;
    Property    _dftengine_options;
    
    string      _output_file;

    string      _reporting;    
    ctp::Logger      _log;
    
  
    void XYZ2Orbitals( Orbitals* _orbitals, string filename);
    //void Coord2Segment(Segment* _segment );
    

};




void DFT::Initialize(Property* options) {

            // update options with the VOTCASHARE defaults   
            UpdateWithDefaults( options, "xtp" );

	    /*    _do_dft_input = false;
            _do_dft_run = false;
            _do_dft_parse = false;
            _do_gwbse = false;
            _do_optimize = false;
            _restart_opt = false;
            */
            
            string key = "options." + Identify();
            _output_file = options->get(key + ".archive").as<string>();
            _reporting =  options->get(key + ".reporting").as<string> ();

            // options for dft package
            string _package_xml = options->get(key + ".package").as<string> ();
            //cout << endl << "... ... Parsing " << _package_xml << endl ;
            load_property_from_xml(_package_options, _package_xml.c_str());
            key = "package";
            _package = _package_options.get(key + ".name").as<string> ();

    

            // options for GWBSE package
key = "options." + Identify();
        _logfile = options->get(key + ".molecule.log").as<string> ();
        _orbfile = options->get(key + ".molecule.orbitals").as<string> ();
            string _dftengine_xml = options->get(key + ".dftengine").as<string> ();
            //cout << endl << "... ... Parsing " << _package_xml << endl ;
            load_property_from_xml(_dftengine_options, _dftengine_xml.c_str());

            
            // job tasks
            /*string _tasks_string = options->get(key + ".tasks").as<string> ();
            if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
            if (_tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
            if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;
            if (_tasks_string.find("optimize") != std::string::npos) _do_optimize = true;
	    */           
 
            
            // something unique
            // _package = options->get(key + ".package").as<string> ();

            // initial coordinates
	    _xyzfile = options->get(key + ".xyz").as<string>();

            
            // get the path to the shared folders with xml files
            char *votca_share = getenv("VOTCASHARE");
            if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
            string xmlFile = string(getenv("VOTCASHARE")) + string("/xtp/packages/") + _package + string("_idft_pair.xml");
            //load_property_from_xml(_package_options, xmlFile);

            // register all QM packages (Gaussian, TURBOMOLE, etc)
            QMPackageFactory::RegisterAll();


           
            
}

bool DFT::Evaluate() {

    
    if ( _reporting == "silent")   _log.setReportLevel( ctp::logERROR ); // only output ERRORS, GEOOPT info, and excited state info for trial geometry
    if ( _reporting == "noisy")    _log.setReportLevel( ctp::logDEBUG ); // OUTPUT ALL THE THINGS
    if ( _reporting == "default")  _log.setReportLevel( ctp::logINFO );  // 
    
    _log.setMultithreading( true );
    _log.setPreface(ctp::logINFO,    "\n... ...");
    _log.setPreface(ctp::logERROR,   "\n... ...");
    _log.setPreface(ctp::logWARNING, "\n... ...");
    _log.setPreface(ctp::logDEBUG,   "\n... ..."); 

    //TLogLevel _ReportLevel = _log.getReportLevel( ); // backup report level

    // Create new orbitals object and fill with atom coordinates
    Orbitals _orbitals;
    XYZ2Orbitals( &_orbitals, _xyzfile );

    vector <ctp::Segment* > _segments;
    // Create a new segment
    ctp::Segment _segment(0, "mol");
                
    //    if (_do_dft_input) {
    //       ReadXYZ( &_segment, _xyzfile );

                ifstream in;
                double x, y, z;
                //int natoms, id;
                int id;
                string label, type;
                vec pos;

                LOG(ctp::logDEBUG,_log) << " Reading molecular coordinates from " << _xyzfile << flush;
                in.open(_xyzfile.c_str(), ios::in);
                if (!in) throw runtime_error(string("Error reading coordinates from: ")
                        + _xyzfile);

                id = 1;
                while (in.good()) { // keep reading until end-of-file
                    in >> type;
                    in >> x;
                    in >> y;
                    in >> z;
                    if (in.eof()) break;
                    // cout << type << ":" << x << ":" << y << ":" << z << endl;

                    // creating atoms and adding them to the molecule
                    ctp::Atom *pAtom = new ctp::Atom(id++, type);
                    vec position(x / 10, y / 10, z / 10); // xyz has Angstrom, votca stores nm
                    pAtom->setPos(position);
                    pAtom->setQMPart(id, position);
                    pAtom->setElement(type);
                    _segment.AddAtom(pAtom);

                }
                in.close();



       _segments.push_back(&_segment);
       //}


    //cout << "here" << endl;
    // get the corresponding object from the QMPackageFactory
        QMPackage *_qmpackage =  QMPackages().Create( _package );
    _qmpackage->setLog( &_log );       
    _qmpackage->Initialize( &_package_options );
    // set the run dir 
    _qmpackage->setRunDir(".");
          _qmpackage->WriteInputFile( _segments, &_orbitals );

	  //       if ( _do_dft_run ){
          _qmpackage->Run();
	  //}

      // parse DFT data, if required
      //if ( _do_dft_parse ){
        LOG(ctp::logDEBUG,_log) << "Parsing DFT data " << _output_file << flush;
        _qmpackage->setOrbitalsFileName( _orbfile );
        //int _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( &_orbitals );
        _qmpackage->ParseOrbitalsFile( &_orbitals );
        _qmpackage->setLogFileName( _logfile );
        _qmpackage->ParseLogFile( &_orbitals );
        //int _parse_log_status = _qmpackage->ParseLogFile( &_orbitals );
        _orbitals.setDFTbasis(_qmpackage->getBasisSetName());
 
        




    // initialize the DFTENGINE
    DFTENGINE _dft;
    _dft.Initialize( &_dftengine_options );
    _dft.setLogger(&_log);
    // RUN
    _dft.Evaluate( &_orbitals );

            
       LOG(ctp::logDEBUG,_log) << "Saving data to " << _output_file << flush;
       std::ofstream ofs( ( _output_file).c_str() );
       boost::archive::binary_oarchive oa( ofs );


       
       oa << _orbitals;
       ofs.close();

     
     
     
     
    Property _summary; 
    _summary.add("output","");
    //Property *_job_output = &_summary.add("output","");
    votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
     
    //ofs (_output_file.c_str(), std::ofstream::out);
    //ofs << *_job_output;    
    //ofs.close();
    
    return true;
}








void DFT::XYZ2Orbitals(Orbitals* _orbitals, string filename){
 
    
               
                ifstream in;
                double x, y, z;
                //int natoms, id;
                string label, type;
                vec pos;

                LOG(ctp::logDEBUG,_log) << " Reading molecular coordinates from " << _xyzfile << flush;
                in.open(_xyzfile.c_str(), ios::in);
                if (!in) throw runtime_error(string("Error reading coordinates from: ")
                        + _xyzfile);

                //id = 1;
                while (in.good()) { // keep reading until end-of-file
                    in >> type;
                    in >> x;
                    in >> y;
                    in >> z;
                    if (in.eof()) break;
                    // cout << type << ":" << x << ":" << y << ":" << z << endl;

                    // creating atoms and adding them to the molecule
                    _orbitals->AddAtom( type, x, y, z );


                }
                in.close();



                return;
}

  
  

}}


#endif
