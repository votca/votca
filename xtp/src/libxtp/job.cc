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
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

// Third party includes
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/propertyiomanipulator.h>

// Local VOTCA includes
#include "votca/xtp/job.h"

namespace votca {
namespace xtp {

Job::Job(const tools::Property &prop) {

  // DEFINED BY USER
  id_ = prop.get("id").as<Index>();
  tag_ = prop.get("tag").as<std::string>();
  input_ = prop.get("input");
  if (prop.exists("status")) {
    status_ = ConvertStatus(prop.get("status").as<std::string>());
  } else {
    status_ = AVAILABLE;
  }

  // GENERATED DURING RUNTIME
  if (prop.exists("host")) {
    host_ = prop.get("host").as<std::string>();
    has_host_ = true;
  }
  if (prop.exists("time")) {
    time_ = prop.get("time").as<std::string>();
    has_time_ = true;
  }
  if (prop.exists("output")) {
    output_ = prop.get("output");
    has_output_ = true;
  }
  if (prop.exists("error")) {
    error_ = prop.get("error").as<std::string>();
    has_error_ = true;
  }
}

Job::Job(Index id, const std::string &tag, const tools::Property &input,
         JobStatus status) {

  id_ = id;
  tag_ = tag;
  input_ = input.get("input");
  status_ = status;
}

std::string Job::ConvertStatus(JobStatus status) const {

  std::string converted;
  switch (status) {
    case AVAILABLE:
      converted = "AVAILABLE";
      break;
    case ASSIGNED:
      converted = "ASSIGNED";
      break;
    case FAILED:
      converted = "FAILED";
      break;
    case COMPLETE:
      converted = "COMPLETE";
      break;
    default:
      throw std::runtime_error("Incomprehensible status (enum)");
  }
  return converted;
}

Job::JobStatus Job::ConvertStatus(std::string status) const {
  JobStatus converted;
  if (status == "AVAILABLE") {
    converted = AVAILABLE;
  } else if (status == "ASSIGNED") {
    converted = ASSIGNED;
  } else if (status == "FAILED") {
    converted = FAILED;
  } else if (status == "COMPLETE") {
    converted = COMPLETE;
  } else {
    throw std::runtime_error("Incomprehensible status: " + status);
  }
  return converted;
}

void Job::Reset() {
  output_ = tools::Property();
  has_output_ = false;
  error_ = "";
  has_error_ = false;
  return;
}

void Job::ToStream(std::ofstream &ofs) const {

  tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 0,
                                      "\t\t");
  std::string tab = "\t";
  ofs << tab << "<job>\n";
  ofs << tab << tab << (boost::format("<id>%1$d</id>\n") % id_).str();
  ofs << tab << tab << (boost::format("<tag>%1$s</tag>\n") % tag_).str();
  ofs << iomXML << input_;
  ofs << tab << tab
      << (boost::format("<status>%1$s</status>\n") % ConvertStatus(status_)).str();

  if (has_host_) {
    ofs << tab << tab << (boost::format("<host>%1$s</host>\n") % host_).str();
  }
  if (has_time_) {
    ofs << tab << tab << (boost::format("<time>%1$s</time>\n") % time_).str();
  }
  if (has_output_) {
    ofs << iomXML << output_;
  }
  if (has_error_) {
    ofs << tab << tab << (boost::format("<error>%1$s</error>\n") % error_).str();
  }
  ofs << tab << "</job>\n";
  return;
}

void Job::UpdateFrom(const Job &ext) {
  status_ = ext.getStatus();
  if (ext.hasHost()) {
    has_host_ = true;
    host_ = ext.getHost();
  }
  if (ext.hasTime()) {
    has_time_ = true;
    time_ = ext.getTime();
  }
  if (ext.hasOutput()) {
    has_output_ = true;
    output_ = ext.getOutput();
  }
  if (ext.hasError()) {
    has_error_ = true;
    error_ = ext.getError();
  }
  return;
}

void Job::UpdateFromResult(const JobResult &res) {
  status_ = res.getStatus();
  if (res.hasOutput()) {
    output_ = res.getOutput();
    has_output_ = true;
  }
  if (res.hasError()) {
    error_ = res.getError();
    has_error_ = true;
  }
  attemptsCount_++;
  return;
}

std::vector<Job> LOAD_JOBS(const std::string &job_file) {

  tools::Property xml;
  xml.LoadFromXML(job_file);

  std::vector<tools::Property *> jobProps = xml.Select("jobs.job");
  std::vector<Job> jobs;
  jobs.reserve(jobProps.size());
  for (tools::Property *prop : jobProps) {
    jobs.push_back(Job(*prop));
  }

  return jobs;
}

void WRITE_JOBS(const std::vector<Job> &jobs, const std::string &job_file) {
  std::ofstream ofs;
  ofs.open(job_file, std::ofstream::out);
  if (!ofs.is_open()) {
    throw std::runtime_error("Bad file handle: " + job_file);
  }
  ofs << "<jobs>" << std::endl;
  for (auto &job : jobs) {
    job.ToStream(ofs);
  }
  ofs << "</jobs>" << std::endl;

  ofs.close();
  return;
}

void UPDATE_JOBS(const std::vector<Job> &from, std::vector<Job> &to,
                 const std::string &thisHost) {
  std::vector<Job>::iterator it_int;
  std::vector<Job>::const_iterator it_ext;

  if (to.size() != from.size()) {
    throw std::runtime_error("Progress file out of sync (::size), abort.");
  }

  for (it_int = to.begin(), it_ext = from.begin(); it_int != to.end();
       ++it_int, ++it_ext) {
    Job &job_int = *it_int;
    const Job &job_ext = *it_ext;
    if (job_int.getId() != job_ext.getId()) {
      throw std::runtime_error("Progress file out of sync (::id), abort.");
    }
    if (job_ext.hasHost() && job_ext.getHost() != thisHost) {
      job_int.UpdateFrom(job_ext);
    }
  }

  return;
}

}  // namespace xtp
}  // namespace votca
