/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Third party includes
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/propertyiomanipulator.h>

// Local VOTCA includes
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/segmentmapper.h"

// Local private VOTCA includes
#include "iexcitoncl.h"

using namespace boost::filesystem;
using namespace votca::tools;

namespace votca {
namespace xtp {

// +++++++++++++++++++++++++++++ //
// IEXCITON MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IEXCITON::ParseSpecificOptions(const tools::Property& options) {

  if (options.exists("states")) {
    std::string parse_string = options.get(".states").as<std::string>();
    statemap_ = FillParseMaps(parse_string);
  }
}

std::map<std::string, QMState> IEXCITON::FillParseMaps(
    const std::string& Mapstring) {
  Tokenizer split_options(Mapstring, ", \t\n");
  std::map<std::string, QMState> type2level;
  for (const std::string& substring : split_options) {
    Tokenizer tok(substring, ":");
    std::vector<std::string> segmentpnumber = tok.ToVector();
    if (segmentpnumber.size() != 2) {
      throw std::runtime_error("Parser iqm: Segment and exciton labels:" +
                               substring + "are not separated properly");
    }
    QMState state = QMState(segmentpnumber[1]);
    if (!state.isTransition()) {
      throw std::runtime_error("State to calculate must be a transition state");
    }
    std::string segmentname = segmentpnumber[0];
    type2level[segmentname] = state;
  }
  return type2level;
}

Job::JobResult IEXCITON::EvalJob(const Topology& top, Job& job,
                                 QMThread& opThread) {

  // report back to the progress observer
  Job::JobResult jres = Job::JobResult();

  // get the logger from the thread
  Logger& pLog = opThread.getLogger();

  // get the information about the job executed by the thread
  Index job_ID = job.getId();
  Property job_input = job.getInput();
  std::vector<Property*> segment_list = job_input.Select("segment");
  Index ID_A = segment_list.front()->getAttribute<Index>("id");
  std::string type_A = segment_list.front()->getAttribute<std::string>("type");
  std::string mps_fileA =
      segment_list.front()->getAttribute<std::string>("mps_file");
  Index ID_B = segment_list.back()->getAttribute<Index>("id");
  std::string type_B = segment_list.back()->getAttribute<std::string>("type");
  std::string mps_fileB =
      segment_list.back()->getAttribute<std::string>("mps_file");

  const Segment& seg_A = top.getSegment(ID_A);
  if (type_A != seg_A.getType()) {
    throw std::runtime_error("SegmentA: type " + seg_A.getType() +
                             " and type in jobfile " + type_A +
                             " do not agree for ID:" + std::to_string(ID_A));
  }
  const Segment& seg_B = top.getSegment(ID_B);
  if (type_B != seg_B.getType()) {
    throw std::runtime_error("SegmentB: type " + seg_B.getType() +
                             " and type in jobfile " + type_B +
                             " do not agree for ID:" + std::to_string(ID_B));
  }
  const QMNBList& nblist = top.NBList();
  const QMPair* pair = nblist.FindPair(&seg_A, &seg_B);
  if (pair == nullptr) {
    throw std::runtime_error(
        "pair between segments " + std::to_string(seg_A.getId()) + ":" +
        std::to_string(seg_B.getId()) + " not found in neighborlist ");
  }

  XTP_LOG(Log::error, pLog) << TimeStamp() << " Evaluating pair " << job_ID
                            << " [" << ID_A << ":" << ID_B << "]" << std::flush;

  StaticMapper map(pLog);
  map.LoadMappingFile(mapfile_);
  StaticSegment seg1 = map.map(*(pair->Seg1()), mps_fileA);
  StaticSegment seg2 = map.map(pair->Seg2PbCopy(), mps_fileB);
  eeInteractor actor;
  double JAB = actor.CalcStaticEnergy(seg1, seg2);
  cutoff_ = 0;
  Property job_summary;
  Property& job_output = job_summary.add("output", "");
  Property& pair_summary = job_output.add("pair", "");
  std::string nameA = seg_A.getType();
  std::string nameB = seg_B.getType();
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
    state = statemap_.at(elementname);
  } catch (std::out_of_range&) {
    std::string errormessage =
        "Map does not have segment of type: " + elementname;
    errormessage += "\n segments in map are:";
    for (const auto& s : statemap_) {
      errormessage += "\n\t" + s.first;
    }
    throw std::runtime_error(errormessage);
  }
  return state;
}

void IEXCITON::WriteJobFile(const Topology& top) {

  std::cout << "\n... ... Writing job file " << jobfile_ << std::flush;
  std::ofstream ofs;
  ofs.open(jobfile_, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("\nERROR: bad file handle: " + jobfile_);
  }
  const QMNBList& nblist = top.NBList();
  Index jobCount = 0;
  if (nblist.size() == 0) {
    std::cout << "\n... ... No pairs in neighbor list, skip." << std::flush;
    return;
  }

  ofs << "<jobs>\n";
  std::string tag = "";

  for (const QMPair* pair : nblist) {
    if (pair->getType() == QMPair::PairType::Excitoncl) {
      Index id1 = pair->Seg1()->getId();
      std::string name1 = pair->Seg1()->getType();
      Index id2 = pair->Seg2()->getId();
      std::string name2 = pair->Seg2()->getType();
      Index id = jobCount;
      QMState state1 = GetElementFromMap(name1);
      QMState state2 = GetElementFromMap(name2);

      std::string mps_file1 =
          (boost::format("MP_FILES/%s_%s.mps") % name1 % state1.ToString())
              .str();
      std::string mps_file2 =
          (boost::format("MP_FILES/%s_%s.mps") % name2 % state2.ToString())
              .str();

      Property Input;
      Property& pInput = Input.add("input", "");
      Property& pSegment1 =
          pInput.add("segment", boost::lexical_cast<std::string>(id1));
      pSegment1.setAttribute<std::string>("type", name1);
      pSegment1.setAttribute<Index>("id", id1);
      pSegment1.setAttribute<std::string>("mps_file", mps_file1);
      Property& pSegment2 =
          pInput.add("segment", boost::lexical_cast<std::string>(id2));
      pSegment2.setAttribute<std::string>("type", name2);
      pSegment2.setAttribute<Index>("id", id2);
      pSegment2.setAttribute<std::string>("mps_file", mps_file2);

      Job job(id, tag, Input, Job::AVAILABLE);
      job.ToStream(ofs);
      jobCount++;
    }
  }

  // CLOSE STREAM
  ofs << "</jobs>\n";
  ofs.close();

  std::cout << "\n... ... In total " << jobCount << " jobs" << std::flush;
}

void IEXCITON::ReadJobFile(Topology& top) {

  Logger log;
  log.setReportLevel(Log::current_level);
  Property xml;
  xml.LoadFromXML(jobfile_);
  std::vector<Property*> jobProps = xml.Select("jobs.job");

  Index updated_jobs = 0;
  Index incomplete_jobs = 0;

  // loop over all jobs = pair records in the job file
  for (const Property* prop : jobProps) {
    // if job produced an output, then continue with analysis
    if (prop->exists("output") && prop->exists("output.pair")) {
      const Property& poutput = prop->get("output.pair");
      Index idA = poutput.getAttribute<Index>("idA");
      Index idB = poutput.getAttribute<Index>("idB");
      Segment& segA = top.getSegment(idA);
      Segment& segB = top.getSegment(idB);
      QMPair* qmp = top.NBList().FindPair(&segA, &segB);

      if (qmp == nullptr) {
        XTP_LOG(Log::error, log)
            << "No pair " << idA << ":" << idB
            << " found in the neighbor list. Ignoring" << std::flush;
      } else if (qmp->getType() == QMPair::Excitoncl) {
        updated_jobs++;
        const Property& pCoupling = poutput.get("Coupling");
        double jAB = pCoupling.getAttribute<double>("jABstatic");
        double Jeff2 = jAB * jAB * tools::conv::ev2hrt * tools::conv::ev2hrt;
        qmp->setJeff2(Jeff2, QMStateType::Singlet);
      }
    } else {
      incomplete_jobs++;
    }
  }

  XTP_LOG(Log::error, log) << "Neighborlist size " << top.NBList().size()
                           << std::flush;

  XTP_LOG(Log::error, log) << "Pairs in jobfile [total:updated:incomplete] "
                           << jobProps.size() << ":" << updated_jobs << ":"
                           << incomplete_jobs << std::flush;
  std::cout << log;
}

}  // namespace xtp
}  // namespace votca
