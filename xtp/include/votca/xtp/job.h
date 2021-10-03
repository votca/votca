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

#pragma once
#ifndef VOTCA_XTP_JOB_H
#define VOTCA_XTP_JOB_H

// Standard includes
#include <fstream>
#include <iostream>
#include <string>

// VOTCA includes
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class Job {

 public:
  enum JobStatus { AVAILABLE, ASSIGNED, FAILED, COMPLETE };

  Job(const tools::Property &prop);
  Job(Index id, const std::string &tag, const tools::Property &input,
      JobStatus status);

  std::string ConvertStatus(JobStatus) const;
  JobStatus ConvertStatus(std::string) const;

  class JobResult {
   public:
    void setStatus(JobStatus stat) { status_ = stat; }
    void setOutput(std::string output) {
      has_output_ = true;
      output_ = tools::Property().add("output", output);
    }
    void setOutput(tools::Property &output) {
      has_output_ = true;
      output_ = output.get("output");
    }

    JobStatus getStatus() const { return status_; }
    bool hasOutput() const { return has_output_; }
    const tools::Property &getOutput() const { return output_; }

    tools::Property &getOutput() { return output_; }
    bool hasError() const { return has_error_; }
    const std::string &getError() const { return error_; }

    void setError(std::string error) {
      has_error_ = true;
      error_ = error;
    }

   private:
    JobStatus status_;
    tools::Property output_;
    bool has_output_ = false;
    std::string error_;
    bool has_error_ = false;
  };

  void Reset();
  void ToStream(std::ofstream &ofs) const;
  void UpdateFrom(const Job &ext);
  void UpdateFromResult(const JobResult &res);

  Index getId() const { return id_; }
  std::string getTag() const { return tag_; }
  tools::Property &getInput() { return input_; }
  const tools::Property &getInput() const { return input_; }
  const JobStatus &getStatus() const { return status_; }
  std::string getStatusStr() const { return ConvertStatus(status_); }

  bool hasHost() const { return has_host_; }
  bool hasTime() const { return has_time_; }
  bool hasOutput() const { return has_output_; }
  bool hasError() const { return has_error_; }

  bool isAvailable() const { return (status_ == AVAILABLE) ? true : false; }
  bool isAssigned() const { return (status_ == ASSIGNED) ? true : false; }
  bool isFailed() const { return (status_ == FAILED) ? true : false; }
  bool isComplete() const { return (status_ == COMPLETE) ? true : false; }
  bool isFresh() const { return (attemptsCount_ < 1) ? true : false; }

  void setStatus(JobStatus stat) { status_ = stat; }
  void setStatus(std::string stat) { status_ = ConvertStatus(stat); }
  void setTime(std::string time) {
    time_ = time;
    has_time_ = true;
  }
  void setHost(std::string host) {
    host_ = host;
    has_host_ = true;
  }
  void setOutput(std::string output) {
    output_ = tools::Property().add("output", output);
    has_output_ = true;
  }

  const std::string &getHost() const {
    assert(has_host_ && "Job has no host");
    return host_;
  }
  const std::string &getTime() const {
    assert(has_time_ && "Job has no time");
    return time_;
  }
  const tools::Property &getOutput() const {
    assert(has_output_ && "Job has no output");
    return output_;
  }
  const std::string &getError() const {
    assert(has_error_ && "Job has no error");
    return error_;
  }

 private:
  // Defined by user
  Index id_;
  std::string tag_;
  JobStatus status_;
  Index attemptsCount_ = 0;
  tools::Property input_;

  // Generated during runtime
  std::string host_;
  bool has_host_ = false;
  std::string time_;
  bool has_time_ = false;
  tools::Property output_;
  bool has_error_ = false;
  bool has_output_ = false;
  std::string error_;
};  // namespace xtp

std::vector<Job> LOAD_JOBS(const std::string &xml_file);
void WRITE_JOBS(const std::vector<Job> &jobs, const std::string &job_file);
void UPDATE_JOBS(const std::vector<Job> &from, std::vector<Job> &to,
                 const std::string &thisHost);
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOB_H
