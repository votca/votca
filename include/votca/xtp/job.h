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

#ifndef VOTCA_XTP_JOB_H
#define VOTCA_XTP_JOB_H

#include <fstream>
#include <iostream>
#include <string>
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class Job {

 public:
  enum JobStatus { AVAILABLE, ASSIGNED, FAILED, COMPLETE };

  Job(const tools::Property &prop);
  Job(int id, std::string &tag, std::string &input, std::string stat);
  Job(int id, std::string &tag, tools::Property &input, JobStatus stat);

  std::string ConvertStatus(JobStatus) const;
  JobStatus ConvertStatus(std::string) const;

  class JobResult {
   public:
    void setStatus(JobStatus stat) { _status = stat; }
    void setOutput(std::string output) {
      _has_output = true;
      _output = tools::Property().add("output", output);
    }
    void setOutput(tools::Property &output) {
      _has_output = true;
      _output = output.get("output");
    }

    JobStatus getStatus() const { return _status; }
    bool hasOutput() const { return _has_output; }
    const tools::Property &getOutput() const { return _output; }
    bool hasError() const { return _has_error; }
    const std::string &getError() const { return _error; }

    void setError(std::string error) {
      _has_error = true;
      _error = error;
    }

   private:
    JobStatus _status;
    tools::Property _output;
    bool _has_output = false;
    std::string _error;
    bool _has_error = false;
  };

  void Reset();
  void ToStream(std::ofstream &ofs, std::string fileformat) const;
  void UpdateFrom(const Job &ext);
  void UpdateFromResult(const JobResult &res);

  int getId() const { return _id; }
  std::string getTag() const { return _tag; }
  tools::Property &getInput() { return _input; }
  const tools::Property &getInput() const { return _input; }
  const JobStatus &getStatus() const { return _status; }
  std::string getStatusStr() const { return ConvertStatus(_status); }

  bool hasHost() const { return _has_host; }
  bool hasTime() const { return _has_time; }
  bool hasOutput() const { return _has_output; }
  bool hasError() const { return _has_error; }

  bool isAvailable() const { return (_status == AVAILABLE) ? true : false; }
  bool isAssigned() const { return (_status == ASSIGNED) ? true : false; }
  bool isFailed() const { return (_status == FAILED) ? true : false; }
  bool isComplete() const { return (_status == COMPLETE) ? true : false; }
  bool isFresh() const { return (_attemptsCount < 1) ? true : false; }

  void setStatus(JobStatus stat) { _status = stat; }
  void setStatus(std::string stat) { _status = ConvertStatus(stat); }
  void setTime(std::string time) {
    _time = time;
    _has_time = true;
  }
  void setHost(std::string host) {
    _host = host;
    _has_host = true;
  }
  void setOutput(std::string output) {
    _output = tools::Property().add("output", output);
    _has_output = true;
  }

  const std::string &getHost() const {
    assert(_has_host && "Job has no host");
    return _host;
  }
  const std::string &getTime() const {
    assert(_has_time && "Job has no time");
    return _time;
  }
  const tools::Property &getOutput() const {
    assert(_has_output && "Job has no output");
    return _output;
  }
  const std::string &getError() const {
    assert(_has_error && "Job has no error");
    return _error;
  }

 private:
  // Defined by user
  int _id;
  std::string _tag;
  JobStatus _status;
  int _attemptsCount = 0;
  tools::Property _input;

  // Generated during runtime
  std::string _host;
  bool _has_host = false;
  std::string _time;
  bool _has_time = false;
  tools::Property _output;
  bool _has_error = false;
  bool _has_output = false;
  std::string _error;
};  // namespace xtp

std::vector<Job> LOAD_JOBS(const std::string &xml_file);
void WRITE_JOBS(const std::vector<Job> &jobs, const std::string &job_file,
                std::string fileformat);
void UPDATE_JOBS(const std::vector<Job> &from, std::vector<Job> &to,
                 const std::string &thisHost);
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOB_H
