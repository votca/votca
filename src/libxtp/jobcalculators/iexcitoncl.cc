/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include "votca/xtp/segmentmapper.h"

#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <votca/tools/constants.h>
#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/eeinteractor.h>
#include <votca/xtp/logger.h>

using namespace std;
using namespace boost::filesystem;
using namespace votca::tools;

namespace votca {
namespace xtp {

// +++++++++++++++++++++++++++++ //
// IEXCITON MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IEXCITON::Initialize(tools::Property& options) {

  string key = "options." + Identify();
  ParseCommonOptions(options);

  if (options.exists(key + ".states")) {
    string parse_string = options.get(key + ".states").as<string>();
    _statemap = FillParseMaps(parse_string);
  }

  cout << "done" << endl;
}

std::map<std::string, QMState> IEXCITON::FillParseMaps(
    const string& Mapstring) {
  Tokenizer split_options(Mapstring, ", \t\n");
  std::map<std::string, QMState> type2level;
  for (const string& substring : split_options) {
    std::vector<string> segmentpnumber;
    Tokenizer tok(substring, ":");
    tok.ToVector(segmentpnumber);
    if (segmentpnumber.size() != 2) {
      throw runtime_error("Parser iqm: Segment and exciton labels:" +
                          substring + "are not separated properly");
    }
    QMState state = QMState(segmentpnumber[1]);
    if (!state.isTransition()) {
      throw std::runtime_error("State to calculate must be a transition state");
    }
    string segmentname = segmentpnumber[0];
    type2level[segmentname] = state;
  }
  return type2level;
}

Job::JobResult IEXCITON::EvalJob(Topology& top, Job& job, QMThread& opThread) {

  // report back to the progress observer
  Job::JobResult jres = Job::JobResult();

  // get the logger from the thread
  Logger& pLog = opThread.getLogger();

  // get the information about the job executed by the thread
  int job_ID = job.getId();
  Property job_input = job.getInput();
  vector<Property*> segment_list = job_input.Select("segment");
  int ID_A = segment_list.front()->getAttribute<int>("id");
  string type_A = segment_list.front()->getAttribute<string>("type");
  string mps_fileA = segment_list.front()->getAttribute<string>("mps_file");
  int ID_B = segment_list.back()->getAttribute<int>("id");
  string type_B = segment_list.back()->getAttribute<string>("type");
  string mps_fileB = segment_list.back()->getAttribute<string>("mps_file");

  Segment& seg_A = top.getSegment(ID_A);
  Segment& seg_B = top.getSegment(ID_B);
  QMNBList& nblist = top.NBList();
  QMPair* pair = nblist.FindPair(&seg_A, &seg_B);
  if (pair == nullptr) {
    throw std::runtime_error(
        "pair between segments " + std::to_string(seg_A.getId()) + ":" +
        std::to_string(seg_B.getId()) + " not found in neighborlist ");
  }

  XTP_LOG_SAVE(logINFO, pLog) << TimeStamp() << " Evaluating pair " << job_ID
                              << " [" << ID_A << ":" << ID_B << "]" << flush;

  StaticMapper map(pLog);
  map.LoadMappingFile(_mapfile);
  StaticSegment seg1 = map.map(*(pair->Seg1()), mps_fileA);
  StaticSegment seg2 = map.map(*(pair->Seg2PbCopy()), mps_fileB);
  eeInteractor actor;
  double JAB = actor.CalcStaticEnergy(seg1, seg2);
  _cutoff = 0;
  Property job_summary;
  Property& job_output = job_summary.add("output", "");
  Property& pair_summary = job_output.add("pair", "");
  string nameA = seg_A.getName();
  string nameB = seg_B.getName();
  pair_summary.setAttribute("idA", ID_A);
  pair_summary.setAttribute("idB", ID_B);
  pair_summary.setAttribute("typeA", nameA);
  pair_summary.setAttribute("typeB", nameB);
  Property& coupling_summary = pair_summary.add("Coupling", "");
  coupling_summary.setAttribute("jABstatic", JAB * tools::conv::hrt2ev);

  jres.setOutput(job_summary);
  jres.setStatus(Job::COMPLETE);

  return jres;
}

QMState IEXCITON::GetElementFromMap(const std::string& elementname) const {
  QMState state;
  try {
    state = _statemap.at(elementname);
  } catch (std::out_of_range& error) {
    std::string errormessage =
        "Map does not have segment of type: " + elementname;
    errormessage += "\n segments in map are:";
    for (const auto& s : _statemap) {
      errormessage += "\n\t" + s.first;
    }
    throw std::runtime_error(errormessage);
  }
  return state;
}

void IEXCITON::WriteJobFile(Topology& top) {

  cout << endl << "... ... Writing job file " << flush;
  std::ofstream ofs;
  ofs.open(_jobfile, std::ofstream::out);
  if (!ofs.is_open())
    throw runtime_error("\nERROR: bad file handle: " + _jobfile);
  QMNBList& nblist = top.NBList();
  int jobCount = 0;
  if (nblist.size() == 0) {
    cout << endl << "... ... No pairs in neighbor list, skip." << flush;
    return;
  }

  ofs << "<jobs>" << endl;
  string tag = "";

  for (QMPair* pair : nblist) {
    if (pair->getType() == QMPair::PairType::Excitoncl) {
      int id1 = pair->Seg1()->getId();
      string name1 = pair->Seg1()->getName();
      int id2 = pair->Seg2()->getId();
      string name2 = pair->Seg2()->getName();
      int id = jobCount;
      QMState state1 = GetElementFromMap(name1);
      QMState state2 = GetElementFromMap(name2);

      string mps_file1 =
          (boost::format("MP_FILES/%s_%s.mps") % name1 % state1.ToString())
              .str();
      string mps_file2 =
          (boost::format("MP_FILES/%s_%s.mps") % name1 % state2.ToString())
              .str();

      Property Input;
      Property& pInput = Input.add("input", "");
      Property& pSegment1 =
          pInput.add("segment", boost::lexical_cast<string>(id1));
      pSegment1.setAttribute<string>("type", name1);
      pSegment1.setAttribute<int>("id", id1);
      pSegment1.setAttribute<string>("mps_file", mps_file1);
      Property& pSegment2 =
          pInput.add("segment", boost::lexical_cast<string>(id2));
      pSegment2.setAttribute<string>("type", name2);
      pSegment2.setAttribute<int>("id", id2);
      pSegment2.setAttribute<string>("mps_file", mps_file2);

      Job job(id, tag, Input, Job::AVAILABLE);
      job.ToStream(ofs, "xml");
      jobCount++;
    }
  }

  // CLOSE STREAM
  ofs << "</jobs>" << endl;
  ofs.close();

  cout << endl << "... ... In total " << jobCount << " jobs" << flush;
}

void IEXCITON::ReadJobFile(Topology& top) {

  Property xml;
  vector<Property*> records;
  // gets the neighborlist from the topology
  QMNBList& nblist = top.NBList();
  int number_of_pairs = nblist.size();
  int current_pairs = 0;
  Logger log;
  log.setReportLevel(logINFO);

  // load the QC results in a vector indexed by the pair ID
  load_property_from_xml(xml, _jobfile);
  vector<Property*> jobProps = xml.Select("jobs.job");
  records.resize(number_of_pairs + 1);

  // to skip pairs which are not in the jobfile
  for (unsigned i = 0; i < records.size(); i++) {
    records[i] = NULL;
  }
  // loop over all jobs = pair records in the job file
  for (Property* prop : jobProps) {
    // if job produced an output, then continue with analysis
    if (prop->exists("output") && prop->exists("output.pair")) {
      current_pairs++;
      Property& poutput = prop->get("output.pair");
      int idA = poutput.getAttribute<int>("idA");
      int idB = poutput.getAttribute<int>("idB");
      Segment& segA = top.getSegment(idA);
      Segment& segB = top.getSegment(idB);
      QMPair* qmp = nblist.FindPair(&segA, &segB);

      if (qmp == NULL) {
        XTP_LOG_SAVE(logINFO, log)
            << "No pair " << idA << ":" << idB
            << " found in the neighbor list. Ignoring" << flush;
      } else {
        records[qmp->getId()] = &(prop->get("output.pair"));
      }
    } else {
      Property thebadone = prop->get("id");
      throw runtime_error("\nERROR: Job file incomplete.\n Job with id " +
                          thebadone.as<string>() +
                          " is not finished. Check your job file for FAIL, "
                          "AVAILABLE, or ASSIGNED. Exiting\n");
    }
  }  // finished loading from the file

  // loop over all pairs in the neighbor list
  XTP_LOG_SAVE(logINFO, log)
      << "Neighborlist size " << top.NBList().size() << flush;
  for (QMPair* pair : top.NBList()) {

    if (records[pair->getId()] == NULL)
      continue;  // skip pairs which are not in the jobfile
    double Jeff2 = 0.0;
    double jAB = 0.0;
    if (pair->getType() == QMPair::Excitoncl) {
      Property* pair_property = records[pair->getId()];
      vector<Property*> pCoupling = pair_property->Select("Coupling");
      for (Property* coup : pCoupling) {
        jAB = coup->getAttribute<double>("jABstatic");
      }
      Jeff2 = jAB * jAB * tools::conv::ev2hrt * tools::conv::ev2hrt;
      pair->setJeff2(Jeff2, QMStateType::Singlet);
    }
  }
  XTP_LOG_SAVE(logINFO, log) << "Pairs [total:updated] " << number_of_pairs
                             << ":" << current_pairs << flush;
  cout << log;
}

}  // namespace xtp
};  // namespace votca
