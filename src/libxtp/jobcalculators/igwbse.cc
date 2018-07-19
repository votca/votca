/*
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include "iqm.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/split.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/constants.h>

using boost::format;
using namespace boost::filesystem;




namespace votca {
  namespace xtp {

   
    void IQM::Initialize(votca::tools::Property* options) {

      // tasks to be done by IBSE: dft_input, dft_run, dft_parse, mgbft, bse_coupling
      _do_dft_input = false;
      _do_dft_run = false;
      _do_dft_parse = false;
      _do_gwbse = false;
      _do_coupling = false;

      _store_dft = false;
      _store_singlets = false;
      _store_triplets = false;
      _store_ehint = false;
      _write_orbfile = false;



      ParseOptionsXML(options);

      // register all QM packages (Gaussian, turbomole, etc))
      QMPackageFactory::RegisterAll();
      return;
    }

    void IGWBSE::ParseOptionsXML(votca::tools::Property *opt) {


      // parsing general ibse options
      string key = "options." + Identify();
      // _energy_difference = opt->get( key + ".degeneracy" ).as< double > ();


      // job tasks
      string _tasks_string = opt->get(key + ".tasks").as<string> ();
      if (_tasks_string.find("input") != std::string::npos) _do_dft_input = true;
      if (_tasks_string.find("dft") != std::string::npos) _do_dft_run = true;
      if (_tasks_string.find("parse") != std::string::npos) _do_dft_parse = true;
      if (_tasks_string.find("gwbse") != std::string::npos) _do_gwbse = true;
      if (_tasks_string.find("coupling") != std::string::npos) _do_coupling = true;

      // storage options
      string _store_string = opt->get(key + ".store").as<string> ();
      if (_store_string.find("dft") != std::string::npos) _store_dft = true;
      if (_store_string.find("singlets") != std::string::npos) _store_singlets = true;
      if (_store_string.find("triplets") != std::string::npos) _store_triplets = true;
      if (_store_string.find("ehint") != std::string::npos) _store_ehint = true;
      if (_store_dft || _store_singlets || _store_triplets || _store_ehint) {
        _write_orbfile = true;
      }
      // options for gwbse
      string _gwbse_xml = opt->get(key + ".gwbse_options").as<string> ();
      //cout << endl << "... ... Parsing " << _package_xml << endl ;
      load_property_from_xml(_gwbse_options, _gwbse_xml.c_str());
      string _coupling_xml = opt->get(key + ".bsecoupling_options").as<string>();
      load_property_from_xml(_coupling_options, _coupling_xml.c_str());

      // options for dft package
      string _package_xml = opt->get(key + ".dftpackage").as<string> ();
      //cout << endl << "... ... Parsing " << _package_xml << endl ;
      load_property_from_xml(_package_options, _package_xml.c_str());
      key = "package";
      _package = _package_options.get(key + ".name").as<string> ();

      // job file specification
      key = "options." + Identify();

      if (opt->exists(key + ".job_file")) {
        _jobfile = opt->get(key + ".job_file").as<string>();
      } else {
        throw std::runtime_error("Job-file not set. Abort.");
      }

      //options for parsing data into sql file   
      key = "options." + Identify() + ".readjobfile";
      if (opt->exists(key + ".singlets")) {
        string _parse_string_s = opt->get(key + ".singlets").as<string> ();
        _singlet_levels = FillParseMaps(_parse_string_s);
      }
      if (opt->exists(key + ".triplets")) {
        string _parse_string_t = opt->get(key + ".triplets").as<string> ();
        _triplet_levels = FillParseMaps(_parse_string_t);
      }

      return;
    }

    std::map<std::string, int> IGWBSE::FillParseMaps(string Mapstring) {
      Tokenizer tok_cleanup(Mapstring, ", \t\n");
      std::vector <std::string> strings_vec;
      tok_cleanup.ToVector(strings_vec);
      std::vector<string>::iterator sit;
      std::map<std::string, int> type2level;
      for (sit = strings_vec.begin(); sit < strings_vec.end(); ++sit) {
        //std::vector<string>temp;
        std::vector<string> segmentpnumber;
        string temp = *sit;
        boost::trim(temp);
        boost::algorithm::split(segmentpnumber, temp, boost::is_any_of(": "), boost::token_compress_on);
        if (segmentpnumber.size() != 2) {

          throw runtime_error("Parser igwbse: Segment and exciton labels are not separated properly");
        }
        if (segmentpnumber[1].size() != 2) {
          throw runtime_error("State identifier " + segmentpnumber[1] + " unknown, right now only states up to number 9 are parsed. s1,s2,t1, etc..");
        }

        //boost::algorithm::split( temp, (*sit), boost::is_any_of(""),boost::token_compress_on );
        int number = boost::lexical_cast<int>(segmentpnumber[1].at(1)) - 1;
        string type = boost::lexical_cast<string>(segmentpnumber[0]);
        type2level[type] = number;
      }
      return type2level;
    }

    void IGWBSE::LoadOrbitals(string file_name, Orbitals& orbitals, ctp::Logger *log) {

      CTP_LOG(ctp::logDEBUG, *log) << "Loading " << file_name << flush;
      try {
        orbitals.ReadFromCpt(file_name);
      } catch (std::runtime_error& error) {
        CTP_LOG(ctp::logERROR, *log) << "Failed loading orbitals from " << file_name << flush;
      }

    }

    ctp::Job::JobResult IGWBSE::EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *opThread) {

      // report back to the progress observer
      ctp::Job::JobResult jres = ctp::Job::JobResult();

      string igwbse_work_dir = "OR_FILES";
      string egwbse_work_dir = "OR_FILES";
      string frame_dir = "frame_" + boost::lexical_cast<string>(top->getDatabaseId());

      bool _run_dft_status = false;
      bool _parse_log_status = false;
      bool _parse_orbitals_status = false;
      bool _calculate_integrals = false;
      stringstream sout;
      string output;


      // get the logger from the thread
      ctp::Logger* pLog = opThread->getLogger();

      // get the information about the job executed by the thread
      int _job_ID = job->getId();
      Property _job_input = job->getInput();
      list<Property*> segment_list = _job_input.Select("segment");
      int ID_A = segment_list.front()->getAttribute<int>("id");
      string type_A = segment_list.front()->getAttribute<string>("type");
      int ID_B = segment_list.back()->getAttribute<int>("id");
      string type_B = segment_list.back()->getAttribute<string>("type");

      // set the folders 
      string _pair_dir = (format("%1%%2%%3%%4%%5%") % "pair" % "_" % ID_A % "_" % ID_B).str();

      path arg_path, arg_pathA, arg_pathB, arg_pathAB;

      string orbFileA = (arg_pathA / egwbse_work_dir / "molecules_gwbse" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_A % ".orb").str()).c_str();
      string orbFileB = (arg_pathB / egwbse_work_dir / "molecules_gwbse" / frame_dir / (format("%1%_%2%%3%") % "molecule" % ID_B % ".orb").str()).c_str();
      string orbFileAB = (arg_pathAB / igwbse_work_dir / "pairs_gwbse" / frame_dir / (format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B % ".orb").str()).c_str();
      string _orb_dir = (arg_path / igwbse_work_dir / "pairs_gwbse" / frame_dir).c_str();

      ctp::Segment *seg_A = top->getSegment(ID_A);
      assert(seg_A->getName() == type_A);

      ctp::Segment *seg_B = top->getSegment(ID_B);
      assert(seg_B->getName() == type_B);

      vector < ctp::Segment* > segments;
      segments.push_back(seg_A);
      segments.push_back(seg_B);


      CTP_LOG(ctp::logINFO, *pLog) << ctp::TimeStamp() << " Evaluating pair "
              << _job_ID << " [" << ID_A << ":" << ID_B << "] out of " <<
              (top->NBList()).size() << flush;

      string _package_append = _package + "_gwbse";

      string _qmpackage_work_dir = (arg_path / igwbse_work_dir / _package_append / frame_dir / _pair_dir).c_str();
      // get the corresponding object from the QMPackageFactory
      QMPackage *_qmpackage = QMPackages().Create(_package);
      // set a log file for the package
      _qmpackage->setLog(pLog);
      // set the run dir 
      _qmpackage->setRunDir(_qmpackage_work_dir);
      // get the package options
      _qmpackage->Initialize(&_package_options);

      // if asked, prepare the input files
      if (_do_dft_input) {
        boost::filesystem::create_directories(_qmpackage_work_dir);
        Orbitals *_orbitalsAB = NULL;
        if (_qmpackage->GuessRequested()) { // do not want to do an SCF loop for a dimer
          CTP_LOG(ctp::logINFO, *pLog) << "Guess requested, reading molecular orbitals" << flush;
          Orbitals _orbitalsA, _orbitalsB;
          _orbitalsAB = new Orbitals();
          // load the corresponding monomer orbitals and prepare the dimer guess 

          // failed to load; wrap-up and finish current job
          try {
            _orbitalsA.ReadFromCpt(orbFileA);
          } catch (std::runtime_error& error) {
            CTP_LOG(ctp::logERROR, *pLog) << "Do input: failed loading orbitals from " << orbFileA << flush;
            cout << *pLog;
            output += "failed on " + orbFileA;
            jres.setOutput(output);
            jres.setStatus(ctp::Job::FAILED);
            delete _qmpackage;
            return jres;
          }
 
          try {
            _orbitalsB.ReadFromCpt(orbFileB);
          } catch (std::runtime_error& error) {
            CTP_LOG(ctp::logERROR, *pLog) << "Do input: failed loading orbitals from " << orbFileB << flush;
            cout << *pLog;
            output += "failed on " + orbFileB;
            jres.setOutput(output);
            jres.setStatus(ctp::Job::FAILED);
            delete _qmpackage;
            return jres;
          }
          CTP_LOG(ctp::logDEBUG, *pLog) << "Constructing the guess for dimer orbitals" << flush;
          Orbitals::PrepareGuess(_orbitalsA, _orbitalsB, *_orbitalsAB);
        }



        // if a pair object is available, take into account PBC, otherwise write as is
        ctp::QMNBList* nblist = &top->NBList();
        ctp::QMPair* pair = nblist->FindPair(seg_A, seg_B);

        if (pair == NULL) {
          vector < ctp::Segment* > segments;
          segments.push_back(seg_A);
          segments.push_back(seg_B);
          CTP_LOG(ctp::logWARNING, *pLog) << "PBCs are not taken into account when writing the coordinate file!" << flush;
          _qmpackage->WriteInputFile(segments, _orbitalsAB);
        } else {
          _qmpackage->WriteInputFilePBC(pair, _orbitalsAB);
        }





        delete _orbitalsAB;
      } // end of the input



      // run the executable
      if (_do_dft_run) {
        _run_dft_status = _qmpackage->Run();
        if (!_run_dft_status) {
          output += "run failed; ";
          CTP_LOG(ctp::logERROR, *pLog) << _qmpackage->getPackageName() << " run failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(ctp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }
      } // end of the run


      // This will be later used to write orbitals of the dimer to a file 
      // SOMETHING TO CLEANUP
      Orbitals _orbitalsAB;
      // parse the log/orbitals files
      if (_do_dft_parse) {
        _parse_log_status = _qmpackage->ParseLogFile(_orbitalsAB);

        if (!_parse_log_status) {
          output += "log incomplete; ";
          CTP_LOG(ctp::logERROR, *pLog) << "LOG parsing failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(ctp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }

        _parse_orbitals_status = _qmpackage->ParseOrbitalsFile(_orbitalsAB);

        if (!_parse_orbitals_status) {
          output += "Orbitals parsing failed; ";
          CTP_LOG(ctp::logERROR, *pLog) << "Orbitals parsing failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(ctp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }




      }// end of the parse orbitals/log
      else {
        LoadOrbitals(orbFileAB, _orbitalsAB, pLog);
      }
      BSECoupling _bsecoupling;
      // do excited states calculation
      if (_do_gwbse) {
        GWBSE _gwbse = GWBSE(_orbitalsAB);
        _gwbse.setLogger(pLog);
        _gwbse.Initialize(&_gwbse_options);
        _gwbse.Evaluate();
    
      } // end of excited state calculation, exciton data is in _orbitalsAB
  
      // calculate the coupling
      Property _job_summary;
      Orbitals _orbitalsA, _orbitalsB;
      if (_do_coupling) {
        // orbitals must be loaded from a file
        if (!_do_gwbse) LoadOrbitals(orbFileAB, _orbitalsAB, pLog);



        try {
          _orbitalsA.ReadFromCpt(orbFileA);
        } catch (std::runtime_error& error) {
          CTP_LOG(ctp::logERROR, *pLog) << "Failed loading orbitals from " << orbFileA << flush;
          cout << *pLog;
          output += "failed on " + orbFileA;
          jres.setOutput(output);
          jres.setStatus(ctp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }

        try {
          _orbitalsB.ReadFromCpt(orbFileB);
        } catch (std::runtime_error& error) {
          CTP_LOG(ctp::logERROR, *pLog) << "Failed loading orbitals from " << orbFileB << flush;
          cout << *pLog;
          output += "failed on " + orbFileB;
          jres.setOutput(output);
          jres.setStatus(ctp::Job::FAILED);
          delete _qmpackage;
          return jres;
        }


        _bsecoupling.setLogger(pLog);
        _bsecoupling.Initialize(&_coupling_options);
        _calculate_integrals = _bsecoupling.CalculateCouplings(&_orbitalsA, &_orbitalsB, &_orbitalsAB);
        // std::cout << _log;

        if (!_calculate_integrals) {
          output += "integrals failed; ";
          CTP_LOG(ctp::logERROR, *pLog) << "Calculating integrals failed" << flush;
          cout << *pLog;
          jres.setOutput(output);
          jres.setStatus(ctp::Job::FAILED);
          return jres;
        }

      }
      Property *_job_output = &_job_summary.add("output", "");
      if (_calculate_integrals) {
        // adding coupling elements 
        _orbitalsAB.setSingletCouplings(_bsecoupling.getJAB_singletstorage());
        _orbitalsAB.setTripletCouplings(_bsecoupling.getJAB_tripletstorage());
        Property *_pair_summary = &_job_output->add("pair", "");
        Property *_type_summary = &_pair_summary->add("type", "");
        _bsecoupling.addoutput(_type_summary, &_orbitalsA, & _orbitalsB);
      }


      votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
      sout << iomXML << _job_summary;
      CTP_LOG(ctp::logINFO, *pLog) << ctp::TimeStamp() << " Finished evaluating pair " << ID_A << ":" << ID_B << flush;
      if (_write_orbfile) {
        // save orbitals 
        boost::filesystem::create_directories(_orb_dir);

        CTP_LOG(ctp::logDEBUG, *pLog) << "Saving orbitals to " << orbFileAB << flush;
        



        if (!_store_dft) {
          _orbitalsAB.AOVxc().resize(0, 0);

          _orbitalsAB.MOCoefficients().resize(0, 0);
        }

        if (!_store_singlets) {
          _orbitalsAB.BSESingletCoefficients().resize(0, 0);
          _orbitalsAB.BSESingletEnergies().resize(0, 0);
        }

        if (!_store_triplets) {
          _orbitalsAB.BSETripletCoefficients().resize(0, 0);
          _orbitalsAB.BSETripletEnergies().resize(0, 0);
        }

        // serialization of electron-hole interaction only if explicitly requested
        if (!_store_ehint) {
          _orbitalsAB.eh_t().resize(0, 0);
          _orbitalsAB.eh_s().resize(0, 0);
        }

        _orbitalsAB.WriteToCpt(orbFileAB);
      } else {
        CTP_LOG(ctp::logDEBUG, *pLog) << "Orb file is not saved according to options " << flush;
      }

      // cleanup whatever is not needed
      _qmpackage->CleanUp();
      delete _qmpackage;

      jres.setOutput(_job_summary);
      jres.setStatus(ctp::Job::COMPLETE);

      return jres;
    }

    void IGWBSE::WriteJobFile(ctp::Topology *top) {

      cout << endl << "... ... Writing job file " << flush;
      std::ofstream ofs;
      ofs.open(_jobfile.c_str(), std::ofstream::out);
      if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);


      ctp::QMNBList::iterator pit;
      ctp::QMNBList &nblist = top->NBList();

      int jobCount = 0;
      if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
      }

      ofs << "<jobs>" << endl;
      string tag = "";

      for (pit = nblist.begin(); pit != nblist.end(); ++pit) {
        if ((*pit)->getType() == ctp::QMPair::Excitoncl) {
          continue;
        }
        //if ((*pit)->HasGhost()){ // Used to only produce jobs concerned with pbcs
        int id1 = (*pit)->Seg1()->getId();
        string name1 = (*pit)->Seg1()->getName();
        int id2 = (*pit)->Seg2()->getId();
        string name2 = (*pit)->Seg2()->getName();

        int id = ++jobCount;

        Property Input;
        Property *pInput = &Input.add("input", "");
        Property *pSegment = &pInput->add("segment", boost::lexical_cast<string>(id1));
        pSegment->setAttribute<string>("type", name1);
        pSegment->setAttribute<int>("id", id1);

        pSegment = &pInput->add("segment", boost::lexical_cast<string>(id2));
        pSegment->setAttribute<string>("type", name2);
        pSegment->setAttribute<int>("id", id2);

        ctp::Job job(id, tag, Input, ctp::Job::AVAILABLE);
        job.ToStream(ofs, "xml");
        //}
      }

      // CLOSE STREAM
      ofs << "</jobs>" << endl;
      ofs.close();

      cout << endl << "... ... In total " << jobCount << " jobs" << flush;
      return;
    }

    /** 
     * Imports electronic couplings with superexchange
     */

    void IGWBSE::ReadJobFile(ctp::Topology *top) {

      Property xml;

      vector<Property*> records;

      // gets the neighborlist from the topology
      ctp::QMNBList &nblist = top->NBList();
      int _number_of_pairs = nblist.size();
      int _current_pairs = 0;
      int _incomplete_jobs = 0;

      // output using logger
      ctp::Logger _log;
      _log.setReportLevel(ctp::logINFO);


      // load the QC results in a vector indexed by the pair ID
      load_property_from_xml(xml, _jobfile);
      list<Property*> jobProps = xml.Select("jobs.job");
      records.resize(nblist.size() + 1);
      //to skip pairs which are not in the jobfile
      for (unsigned i = 0; i < records.size(); i++) {
        records[i] = NULL;
      }
      // loop over all jobs = pair records in the job file
      for (list<Property*> ::iterator it = jobProps.begin(); it != jobProps.end(); ++it) {

        //int level_segA=1;
        //int level_segB=1;

        // if job produced an output, then continue with analysis
        if ((*it)->exists("output") && (*it)->exists("output.pair")) {

          // get the output records
          Property poutput = (*it)->get("output.pair");
          // job file is stupid, because segment ids are only in input have to get them out l
          list<Property*> pinput = (*it)->Select("input.segment");
          vector<int> id;
          for (list<Property*> ::iterator iit = pinput.begin(); iit != pinput.end(); ++iit) {
            id.push_back((*iit)->getAttribute<int>("id"));
          }
          if (id.size() != 2) throw std::runtime_error("Getting pair ids from jobfile failed, check jobfile.");

          double idA = id[0];
          double idB = id[1];

          // segments which correspond to these ids           
          ctp::Segment *segA = top->getSegment(idA);
          ctp::Segment *segB = top->getSegment(idB);
          // pair that corresponds to the two segments
          ctp::QMPair *qmp = nblist.FindPair(segA, segB);

          if (qmp == NULL) { // there is no pair in the neighbor list with this name
            CTP_LOG_SAVE(ctp::logINFO, _log) << "No pair " << idA << ":" << idB << " found in the neighbor list. Ignoring" << flush;
          } else {
            //cout << "Store in record: " <<  idA << ":" << idB << flush; 
            records[qmp->getId()] = &((*it)->get("output.pair.type"));

          }
        } else {
          throw runtime_error("\nERROR: Job file incomplete.\n Check your job file for FAIL, AVAILABLE, or ASSIGNED. Exiting\n");
        }
      } // finished loading from the file


      // loop over all pairs in the neighbor list
      std::cout << "Neighborlist size " << top->NBList().size() << std::endl;
      for (ctp::QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {

        ctp::QMPair *pair = *ipair;
        if (records[ pair->getId() ] == NULL) continue; //skip pairs which are not in the jobfile

        ctp::Segment* segmentA = pair->Seg1();
        ctp::Segment* segmentB = pair->Seg2();



        //cout << "Processing pair " << segmentA->getId() << ":" << segmentB->getId() << flush;

        ctp::QMPair::PairType _ptype = pair->getType();
        Property* pair_property = records[ pair->getId() ];



        // If a pair is of a direct type 
        if (_ptype == ctp::QMPair::PairType::Hopping || _ptype == ctp::QMPair::PairType::SuperExchangeAndHopping) {
          bool foundsinglet = false;
          bool foundtriplet = false;

          if (pair_property->exists("singlets")) {

            //bool found=false;
            double coupling;
            list<Property*> singlets = pair_property->Select("singlets.coupling");
            int stateA = _singlet_levels[segmentA->getName()];
            int stateB = _singlet_levels[segmentB->getName()];
            for (list<Property*> ::iterator iit = singlets.begin(); iit != singlets.end(); ++iit) {
              int state1 = (*iit)->getAttribute<int>("excitonA");
              int state2 = (*iit)->getAttribute<int>("excitonB");
              if (state1 == stateA && state2 == stateB) {
                coupling = boost::lexical_cast<double>((*iit)->value());
                pair->setJeff2(coupling*coupling, 2);
                pair->setIsPathCarrier(true, 2);
                foundsinglet = true;
                break;
              }
            }
          }
          if (pair_property->exists("triplets")) {

            //bool found=false;
            double coupling;
            list<Property*> triplets = pair_property->Select("triplets.coupling");
            int stateA = _triplet_levels[segmentA->getName()];
            int stateB = _triplet_levels[segmentB->getName()];
            for (list<Property*> ::iterator iit = triplets.begin(); iit != triplets.end(); ++iit) {
              int state1 = (*iit)->getAttribute<int>("excitonA");
              int state2 = (*iit)->getAttribute<int>("excitonB");
              if (state1 == stateA && state2 == stateB) {
                coupling = boost::lexical_cast<double>((*iit)->value());
                pair->setJeff2(coupling*coupling, 3);
                pair->setIsPathCarrier(true, 3);
                foundtriplet = true;
                break;
              }
            }
          }
          if (foundsinglet || foundtriplet) {
            _current_pairs++;
          }
        } else {
          cout << "WARNING Pair " << pair->getId() << " is not of any of the Hopping or SuperExchangeAndHopping type, what did you do to the jobfile?" << flush;
        }

        //cout << endl;

      }

      CTP_LOG_SAVE(ctp::logINFO, _log) << "Pairs [total:updated] " << _number_of_pairs << ":" << _current_pairs << " Incomplete jobs: " << _incomplete_jobs << flush;
      cout << _log;
      return;
    }

  }
};
