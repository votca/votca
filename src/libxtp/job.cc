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
/// For an earlier history see ctp repo commit
/// 77795ea591b29e664153f9404c8655ba28dc14e9

#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <votca/tools/propertyiomanipulator.h>
#include <votca/xtp/job.h>

namespace votca {
namespace xtp {
using boost::format;
Job::Job(const tools::Property &prop) {

  // DEFINED BY USER
  _id = prop.get("id").as<int>();
  _tag = prop.get("tag").as<std::string>();
  _input = prop.get("input");
  if (prop.exists("status"))
    _status = ConvertStatus(prop.get("status").as<std::string>());
  else
    _status = AVAILABLE;

  // GENERATED DURING RUNTIME
  if (prop.exists("host")) {
    _host = prop.get("host").as<std::string>();
    _has_host = true;
  }
  if (prop.exists("time")) {
    _time = prop.get("time").as<std::string>();
    _has_time = true;
  }
  if (prop.exists("output")) {
    _output = prop.get("output");
    _has_output = true;
  }
  if (prop.exists("error")) {
    _error = prop.get("error").as<std::string>();
    _has_error = true;
  }
}

Job::Job(int id, std::string &tag, std::string &inputstr, std::string status) {

  _id = id;
  _tag = tag;
  tools::Property input("input", inputstr, "");
  _input = input;
  _status = ConvertStatus(status);
}

Job::Job(int id, std::string &tag, tools::Property &input, JobStatus status) {

  _id = id;
  _tag = tag;
  _input = input.get("input");
  _status = status;
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
  if (status == "AVAILABLE")
    converted = AVAILABLE;
  else if (status == "ASSIGNED")
    converted = ASSIGNED;
  else if (status == "FAILED")
    converted = FAILED;
  else if (status == "COMPLETE")
    converted = COMPLETE;
  else
    throw std::runtime_error("Incomprehensible status: " + status);
  return converted;
}

void Job::Reset() {
  _output = tools::Property();
  _has_output = false;
  _error = "";
  _has_error = false;
  return;
}

void Job::ToStream(std::ofstream &ofs) const {

  tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 0,
                                      "\t\t");
  std::string tab = "\t";
  ofs << tab << "<job>\n";
  ofs << tab << tab << (format("<id>%1$d</id>\n") % _id).str();
  ofs << tab << tab << (format("<tag>%1$s</tag>\n") % _tag).str();
  ofs << iomXML << _input;
  ofs << tab << tab
      << (format("<status>%1$s</status>\n") % ConvertStatus(_status)).str();

  if (_has_host)
    ofs << tab << tab << (format("<host>%1$s</host>\n") % _host).str();
  if (_has_time)
    ofs << tab << tab << (format("<time>%1$s</time>\n") % _time).str();
  if (_has_output) ofs << iomXML << _output;
  if (_has_error)
    ofs << tab << tab << (format("<error>%1$s</error>\n") % _error).str();
  ofs << tab << "</job>\n";
  return;
}

void Job::UpdateFrom(const Job &ext) {
  _status = ext.getStatus();
  if (ext.hasHost()) {
    _has_host = true;
    _host = ext.getHost();
  }
  if (ext.hasTime()) {
    _has_time = true;
    _time = ext.getTime();
  }
  if (ext.hasOutput()) {
    _has_output = true;
    _output = ext.getOutput();
  }
  if (ext.hasError()) {
    _has_error = true;
    _error = ext.getError();
  }
  return;
}

void Job::UpdateFromResult(const JobResult &res) {
  _status = res.getStatus();
  if (res.hasOutput()) {
    _output = res.getOutput();
    _has_output = true;
  }
  if (res.hasError()) {
    _error = res.getError();
    _has_error = true;
  }
  _attemptsCount++;
  return;
}

std::vector<Job> LOAD_JOBS(const std::string &job_file) {

  tools::Property xml;
  load_property_from_xml(xml, job_file);

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

  if (to.size() != from.size())
    throw std::runtime_error("Progress file out of sync (::size), abort.");

  for (it_int = to.begin(), it_ext = from.begin(); it_int != to.end();
       ++it_int, ++it_ext) {
    Job &job_int = *it_int;
    const Job &job_ext = *it_ext;
    if (job_int.getId() != job_ext.getId())
      throw std::runtime_error("Progress file out of sync (::id), abort.");
    if (job_ext.hasHost() && job_ext.getHost() != thisHost)
      job_int.UpdateFrom(job_ext);
  }

  return;
}

}  // namespace xtp
}  // namespace votca
