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

using boost::format;
namespace votca {
namespace xtp {

Job::Job(tools::Property *prop)
    : _has_host(false),
      _has_time(false),
      _has_error(false),
      _has_output(false),
      _has_sqlcmd(false) {

  // DEFINED BY USER
  _id = prop->get("id").as<int>();
  _tag = prop->get("tag").as<std::string>();
  _input = prop->get("input");
  _attemptsCount = 0;

  if (prop->exists("status"))
    _status = ConvertStatus(prop->get("status").as<std::string>());
  else
    _status = AVAILABLE;

  if (prop->exists("sqlcmd")) {
    _sqlcmd = prop->get("sqlcmd").as<std::string>();
    _has_sqlcmd = true;
  }

  // GENERATED DURING RUNTIME
  if (prop->exists("host")) {
    _host = prop->get("host").as<std::string>();
    _has_host = true;
  }
  if (prop->exists("time")) {
    _time = prop->get("time").as<std::string>();
    _has_time = true;
  }
  if (prop->exists("output")) {
    _output = prop->get("output");
    _has_output = true;
  }
  if (prop->exists("error")) {
    _error = prop->get("error").as<std::string>();
    _has_error = true;
  }
}

Job::Job(int id, std::string &tag, std::string &inputstr, std::string status)
    : _has_host(false),
      _has_time(false),
      _has_error(false),
      _has_output(false),
      _has_sqlcmd(false) {

  _id = id;
  _tag = tag;
  tools::Property input("input", inputstr, "");
  _input = input;
  _status = ConvertStatus(status);
  _attemptsCount = 0;
}

Job::Job(int id, std::string &tag, tools::Property &input, JobStatus status)
    : _has_host(false),
      _has_time(false),
      _has_error(false),
      _has_output(false),
      _has_sqlcmd(false) {

  _id = id;
  _tag = tag;
  _input = input.get("input");
  _status = status;
  _attemptsCount = 0;
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

void Job::ToStream(std::ofstream &ofs, std::string fileformat) {

  tools::PropertyIOManipulator iomXML(tools::PropertyIOManipulator::XML, 0,
                                      "\t\t");

  if (fileformat == "xml") {
    std::string tab = "\t";

    ofs << tab << "<job>\n";
    ofs << tab << tab << (format("<id>%1$d</id>\n") % _id).str();
    ofs << tab << tab << (format("<tag>%1$s</tag>\n") % _tag).str();
    ofs << iomXML << _input;
    ofs << tab << tab
        << (format("<status>%1$s</status>\n") % ConvertStatus(_status)).str();

    if (_has_sqlcmd)
      ofs << tab << tab << (format("<sqlcmd>%1$s</sqlcmd>\n") % _sqlcmd).str();
    if (_has_host)
      ofs << tab << tab << (format("<host>%1$s</host>\n") % _host).str();
    if (_has_time)
      ofs << tab << tab << (format("<time>%1$s</time>\n") % _time).str();
    if (_has_output) ofs << iomXML << _output;
    if (_has_error)
      ofs << tab << tab << (format("<error>%1$s</error>\n") % _error).str();
    ofs << tab << "</job>\n";
  } else if (fileformat == "tab") {
    std::string time = _time;
    if (!_has_time) time = "__:__";
    std::string host = _host;
    if (!_has_host) host = "__:__";
    std::string status = ConvertStatus(_status);

    std::stringstream input;
    std::stringstream output;

    input << iomXML << _input;
    output << iomXML << _output;
    ofs << (format("%4$10s %5$20s %6$10s %1$5d %2$10s %3$30s %7$s %8$s\n") %
            _id % _tag % input.str() % status % host % time % _error %
            output.str())
               .str();
  } else {
    assert(false);
  }

  return;
}

void Job::UpdateFrom(Job *ext) {

  _status = ext->getStatus();
  if (ext->hasHost()) {
    _has_host = true;
    _host = ext->getHost();
  }
  if (ext->hasTime()) {
    _has_time = true;
    _time = ext->getTime();
  }
  if (ext->hasOutput()) {
    _has_output = true;
    _output = ext->getOutput();
  }
  if (ext->hasError()) {
    _has_error = true;
    _error = ext->getError();
  }

  return;
}

void Job::SaveResults(JobResult *res) {

  _status = res->_status;
  if (res->_has_output) {
    _output = res->_output;
    _has_output = true;
  }
  if (res->_has_error) {
    _error = res->_error;
    _has_error = true;
  }

  _attemptsCount += 1;

  return;
}

}  // namespace xtp
}  // namespace votca
