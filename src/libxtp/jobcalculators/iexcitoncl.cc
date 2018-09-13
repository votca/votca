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


#include "iexcitoncl.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/constants.h>
#include <votca/tools/propertyiomanipulator.h>

#include <votca/xtp/apolarsite.h>
#include <votca/xtp/polarseg.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/xinteractor.h>

using namespace boost::filesystem;
using namespace votca::tools;

namespace votca {
    namespace xtp {

        // +++++++++++++++++++++++++++++ //
        // IEXCITON MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        void IEXCITON::Initialize(tools::Property* options) {

            cout << endl
                    << "... ... Initialized with " << _nThreads << " threads. "
                    << flush;

            _maverick = (_nThreads == 1) ? true : false;

            _induce = false;
            _epsilon = 1;
            _cutoff = -1;

            string key = "options." + Identify();

            if (options->exists(key + ".job_file")) {
                _jobfile = options->get(key + ".job_file").as<string>();
            } else {
                throw std::runtime_error("Job-file not set. Abort.");
            }
            key = "options." + Identify();
            if (options->exists(key + ".mapping")) {
                _xml_file = options->get(key + ".mapping").as<string>();
            } else {
                throw std::runtime_error("Mapping-file not set. Abort.");
            }
            if (options->exists(key + ".emp_file")) {
                _emp_file = options->get(key + ".emp_file").as<string>();
            } else {
                throw std::runtime_error("Emp-file not set. Abort.");
            }
            
            
            if (options->exists(key + ".states")) {
                string parse_string = options->get(key + ".states").as<string>();
                _statemap = FillParseMaps(parse_string);    
            } 
            if (options->exists(key + ".epsilon")) {
                _epsilon = options->get(key + ".epsilon").as<double>();
            } else {
                _epsilon = 1;
            }
            if (options->exists(key + ".cutoff")) {
                _cutoff = options->get(key + ".cutoff").as<double>();
            } else {
                _cutoff = -1;
            }

            if (options->exists(key + ".induce")) {
                _induce = options->get(key + ".induce").as<bool>();
            }

            cout << "done" << endl;
        }
        
        
      std::map<std::string, QMState> IEXCITON::FillParseMaps(const string& Mapstring) {
      Tokenizer split_options(Mapstring, ", \t\n");
      std::map<std::string, QMState> type2level;
      for (const string& substring : split_options) {
        std::vector<string> segmentpnumber;
        Tokenizer tok(substring, ":");
        tok.ToVector(segmentpnumber);
        if (segmentpnumber.size() != 2) {
          throw runtime_error("Parser iqm: Segment and exciton labels:" + substring + "are not separated properly");
        }
        QMState state=QMState(segmentpnumber[1]);
        if(!state.isTransition()){
            throw std::runtime_error("State to calculate must be a transition state");
        }
        string segmentname = segmentpnumber[0];
        type2level[segmentname] = state;
      }
      return type2level;
    }

        void IEXCITON::PreProcess(xtp::Topology *top) {

            // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
            cout << endl << "... ... Initialize MPS-mapper: " << flush;

            _mps_mapper.GenerateMap(_xml_file, _emp_file, top);
        }

       
        xtp::Job::JobResult IEXCITON::EvalJob(xtp::Topology *top, xtp::Job *job, xtp::QMThread *opThread) {

            // report back to the progress observer
            xtp::Job::JobResult jres = xtp::Job::JobResult();

            // get the logger from the thread
            xtp::Logger* pLog = opThread->getLogger();

            // get the information about the job executed by the thread
            int job_ID = job->getId();
            Property job_input = job->getInput();
            list<Property*> segment_list = job_input.Select("segment");
            int ID_A = segment_list.front()->getAttribute<int>("id");
            string type_A = segment_list.front()->getAttribute<string>("type");
            string mps_fileA = segment_list.front()->getAttribute<string>("mps_file");
            int ID_B = segment_list.back()->getAttribute<int>("id");
            string type_B = segment_list.back()->getAttribute<string>("type");
            string mps_fileB = segment_list.back()->getAttribute<string>("mps_file");

            xtp::Segment *seg_A = top->getSegment(ID_A);
            xtp::Segment *seg_B = top->getSegment(ID_B);

            XTP_LOG(xtp::logINFO, *pLog) << xtp::TimeStamp() << " Evaluating pair "
                    << job_ID << " [" << ID_A << ":" << ID_B << "]" << flush;

            vector<xtp::APolarSite*> seg_A_raw = xtp::APS_FROM_MPS(mps_fileA, 0, opThread);
            vector<xtp::APolarSite*> seg_B_raw = xtp::APS_FROM_MPS(mps_fileB, 0, opThread);

            xtp::PolarSeg seg_A_polar = *(_mps_mapper.MapPolSitesToSeg(seg_A_raw, seg_A));
            xtp::PolarSeg seg_B_polar = *(_mps_mapper.MapPolSitesToSeg(seg_B_raw, seg_B));

            double JAB = EvaluatePair(top, &seg_A_polar, &seg_B_polar, pLog);

            for (xtp::APolarSite* site : seg_A_raw) {
                delete site;
            }
            seg_A_raw.clear();
            for (xtp::APolarSite* site : seg_B_raw) {
                delete site;
            }
            seg_B_raw.clear();

            Property job_summary;
            Property& job_output = job_summary.add("output", "");
            Property& pair_summary = job_output.add("pair", "");
            string nameA = seg_A->getName();
            string nameB = seg_B->getName();
            pair_summary.setAttribute("idA", ID_A);
            pair_summary.setAttribute("idB", ID_B);
            pair_summary.setAttribute("typeA", nameA);
            pair_summary.setAttribute("typeB", nameB);
            Property& coupling_summary = pair_summary.add("Coupling", "");
            coupling_summary.setAttribute("jABstatic", JAB);

            jres.setOutput(job_summary);
            jres.setStatus(xtp::Job::COMPLETE);

            return jres;
        }
        
     QMState IEXCITON::GetElementFromMap(const std::string& elementname )const{
      QMState state;
      try{
        state = _statemap.at(elementname);
      }
      catch (std::out_of_range& error) {
        std::string errormessage="Map does not have segment of type: "+elementname;
        errormessage+="\n segments in map are:";
        for(const auto& s:_statemap){
         errormessage+="\n\t"+s.first;
        }
        throw std::runtime_error(errormessage);
      }
      return state;
    }

        double IEXCITON::EvaluatePair(xtp::Topology *top, xtp::PolarSeg* Seg1, xtp::PolarSeg* Seg2, xtp::Logger* pLog) {

            xtp::XInteractor actor;
            actor.ResetEnergy();
            Seg1->CalcPos();
            Seg2->CalcPos();
            vec s = top->PbShortestConnect(Seg1->getPos(), Seg2->getPos()) + Seg1->getPos() - Seg2->getPos();

            double E = 0.0;
            for (xtp::APolarSite* site1 : *Seg1) {
                for (xtp::APolarSite* site2 : *Seg2) {
                    actor.BiasIndu(*site1, *site2, s);
                    site1->Depolarize();
                    site2->Depolarize();
                    E += actor.E_f(*site1, *site2);
                }
            }

            if (_cutoff >= 0) {
                if (abs(s) > _cutoff) {
                    E = E / _epsilon;
                }
            }
            return E * conv::int2eV;
        }

        void IEXCITON::WriteJobFile(xtp::Topology *top) {

            cout << endl << "... ... Writing job file " << flush;
            std::ofstream ofs;
            ofs.open(_jobfile.c_str(), std::ofstream::out);
            if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);
            xtp::QMNBList &nblist = top->NBList();
            int jobCount = 0;
            if (nblist.size() == 0) {
                cout << endl << "... ... No pairs in neighbor list, skip." << flush;
                return;
            }

            ofs << "<jobs>" << endl;
            string tag = "";

            for (xtp::QMPair* pair:nblist) {
                if (pair->getType() == 3) {
                    int id1 = pair->Seg1()->getId();
                    string name1 = pair->Seg1()->getName();
                    int id2 = pair->Seg2()->getId();
                    string name2 = pair->Seg2()->getName();
                    int id = ++jobCount;
                    QMState state1=GetElementFromMap(name1);
                    QMState state2=GetElementFromMap(name2);
                    
                    string mps_file1 = (boost::format("MP_FILES/%s_%s.mps") % name1 % state1.ToString()).str();
                    string mps_file2 = (boost::format("MP_FILES/%s_%s.mps") % name1 % state2.ToString()).str();

                    Property Input;
                    Property& pInput = Input.add("input", "");
                    Property& pSegment1 = pInput.add("segment", boost::lexical_cast<string>(id1));
                    pSegment1.setAttribute<string>("type", name1);
                    pSegment1.setAttribute<int>("id", id1);
                    pSegment1.setAttribute<string>("mps_file", mps_file1);
                    Property& pSegment2 = pInput.add("segment", boost::lexical_cast<string>(id2));
                    pSegment2.setAttribute<string>("type", name2);
                    pSegment2.setAttribute<int>("id", id2);
                    pSegment2.setAttribute<string>("mps_file", mps_file2);

                    xtp::Job job(id, tag, Input, xtp::Job::AVAILABLE);
                    job.ToStream(ofs, "xml");
                }
            }

            // CLOSE STREAM
            ofs << "</jobs>" << endl;
            ofs.close();

            cout << endl << "... ... In total " << jobCount << " jobs" << flush;

        }

        void IEXCITON::ReadJobFile(xtp::Topology *top) {

            Property xml;
            vector<Property*> records;
            // gets the neighborlist from the topology
            xtp::QMNBList &nblist = top->NBList();
            int number_of_pairs = nblist.size();
            int current_pairs = 0;
            xtp::Logger log;
            log.setReportLevel(xtp::logINFO);

            // load the QC results in a vector indexed by the pair ID
            load_property_from_xml(xml, _jobfile);
            list<Property*> jobProps = xml.Select("jobs.job");
            records.resize(number_of_pairs + 1);

            //to skip pairs which are not in the jobfile
            for (unsigned i = 0; i < records.size(); i++) {
                records[i] = NULL;
            }
            // loop over all jobs = pair records in the job file
            for (Property* prop:jobProps) {
                // if job produced an output, then continue with analysis
                if (prop->exists("output") && prop->exists("output.pair")) {
                    current_pairs++;
                    Property& poutput = prop->get("output.pair");
                    int idA = poutput.getAttribute<int>("idA");
                    int idB = poutput.getAttribute<int>("idB");      
                    xtp::Segment *segA = top->getSegment(idA);
                    xtp::Segment *segB = top->getSegment(idB);
                    xtp::QMPair *qmp = nblist.FindPair(segA, segB);

                    if (qmp == NULL) { 
                        XTP_LOG_SAVE(xtp::logINFO, log) << "No pair " << idA << ":" << idB << " found in the neighbor list. Ignoring" << flush;
                    } else {
                        records[qmp->getId()] = &(prop->get("output.pair"));
                    }
                } else {
                    Property thebadone = prop->get("id");
                    throw runtime_error("\nERROR: Job file incomplete.\n Job with id " + thebadone.as<string>() + " is not finished. Check your job file for FAIL, AVAILABLE, or ASSIGNED. Exiting\n");
                }
            } // finished loading from the file


            // loop over all pairs in the neighbor list
            XTP_LOG_SAVE(xtp::logINFO, log) << "Neighborlist size " << top->NBList().size() << flush;
            for (xtp::QMPair *pair: top->NBList()) {

                if (records[ pair->getId() ] == NULL) continue; //skip pairs which are not in the jobfile
                double Jeff2 = 0.0;
                double jAB = 0.0;
                if (pair->getType() == xtp::QMPair::Excitoncl) {
                    Property* pair_property = records[ pair->getId() ];
                    list<Property*> pCoupling = pair_property->Select("Coupling");
                    for (Property* coup:pCoupling) {
                        jAB = coup->getAttribute<double>("jABstatic");
                    }
                    Jeff2 = jAB*jAB;
                    pair->setJeff2(Jeff2, 2);
                    pair->setIsPathCarrier(true, 2);
                }
            }
            XTP_LOG_SAVE(xtp::logINFO, log) << "Pairs [total:updated] " << number_of_pairs << ":" << current_pairs << flush;
            cout << log;
        }

    }
};
