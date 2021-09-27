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

#pragma once
#ifndef VOTCA_XTP_LOGGER_H
#define VOTCA_XTP_LOGGER_H

// Standard includes
#include <chrono>
#include <iostream>
#include <sstream>

// VOTCA includes
#include <votca/tools/globals.h>

namespace votca {
namespace xtp {

/*
 * Macros to use the Logger: XTP_LOG(level,logger) << message
 */
#define XTP_LOG(level, log)           \
  if (level > (log).getReportLevel()) \
    ;                                 \
  else                                \
    (log)(level)

/*
 * Custom buffer to store messages
 */
class LogBuffer final : public std::stringbuf {

 public:
  LogBuffer() : std::stringbuf() { LogLevel_ = Log::current_level; }

  // sets the log level (needed for output)
  void setLogLevel(Log::Level LogLevel) { LogLevel_ = LogLevel; }

  // sets Multithreading (buffering required)
  void setMultithreading(bool maverick) { maverick_ = maverick; }

  // sets preface strings for Log::error, Log::warning, ...
  void setPreface(Log::Level level, std::string preface) {
    switch (level) {

      case Log::Level::error:
        errorPreface_ = preface;
        break;
      case Log::Level::warning:
        warnPreface_ = preface;
        break;
      case Log::Level::info:
        infoPreface_ = preface;
        break;
      case Log::Level::debug:
        dbgPreface_ = preface;
        break;
    }
  }

  void EnablePreface() { writePreface_ = true; }
  void DisablePreface() { writePreface_ = false; }

  // flushes all collected messages
  void FlushBuffer() {
    std::cout << stringStream_.str();
    stringStream_.str("");
  }

  // returns the pointer to the collected messages
  std::string Messages() {
    std::string messages_ = stringStream_.str();
    stringStream_.str("");
    return messages_;
  }

 private:
  // Log Level (WARNING, INFO, etc)
  Log::Level LogLevel_ = Log::Level::error;

  // temporary buffer to store messages
  std::ostringstream stringStream_;

  // Multithreading
  bool maverick_ = true;

  std::string errorPreface_ = "\n ERROR   ";
  std::string warnPreface_ = "\n WARNING ";
  std::string infoPreface_ = "\n         ";
  std::string dbgPreface_ = "\n DEBUG   ";
  bool writePreface_ = true;

 protected:
  int sync() {

    std::ostringstream message_;

    if (writePreface_) {
      switch (LogLevel_) {
        case Log::Level::error:
          message_ << errorPreface_;
          break;
        case Log::Level::warning:
          message_ << warnPreface_;
          break;
        case Log::Level::info:
          message_ << infoPreface_;
          break;
        case Log::Level::debug:
          message_ << dbgPreface_;
          break;
      }
    }

    if (!maverick_) {
      // collect all messages of one thread
      stringStream_ << message_.str() << " " << str();
    } else {
      // if only one thread outputs, flush immediately
      std::cout << message_.str() << " " << str() << std::flush;
    }
    message_.str("");
    str("");
    return 0;
  }
};

/** \class Logger
 *   \brief Logger is used for thread-safe output of messages
 *
 *  Logger writes messages into LogBuffer.
 *  Inheritance from ostream allows to use overloaded << and >> for writing.
 *  Example:
 *
 *  \code
 *  #include <logger.h>
 *  Logger log; // create a logger object
 *  log.setReportLevel(Log::error); // output only log messages starting from a
 *  level XTP_LOG(Log::error,*log) << "Error detected" << flush; // write to
 * the logger at an ERROR level cout << log; // output logger content to
 * standard output \endcode
 *
 *  Logger has four predefined log levels: error, warning, info,
 * debug.
 */
class Logger final : public std::ostream {

  friend std::ostream &operator<<(std::ostream &log_out, Logger &logger) {
    log_out << logger.Messages();
    return log_out;
  }

 public:
  Logger() : std::ostream(&buffer_), ReportLevel_(Log::current_level) {
    setMultithreading(maverick_);
  }
  Logger(Log::Level ReportLevel)
      : std::ostream(&buffer_), ReportLevel_(ReportLevel) {
    setMultithreading(maverick_);
  }

  Logger &operator()(Log::Level LogLevel) {
    buffer_.setLogLevel(LogLevel);
    return *this;
  }

  void setReportLevel(Log::Level ReportLevel) { ReportLevel_ = ReportLevel; }
  void setMultithreading(bool maverick) {
    maverick_ = maverick;
    buffer_.setMultithreading(maverick_);
  }
  bool isMaverick() const { return maverick_; }

  Log::Level getReportLevel() const { return ReportLevel_; }

  void setPreface(Log::Level level, const std::string &preface) {
    buffer_.setPreface(level, preface);
  }

  void setCommonPreface(const std::string &preface) {
    setPreface(Log::info, preface);
    setPreface(Log::error, preface);
    setPreface(Log::debug, preface);
    setPreface(Log::warning, preface);
  }

  void EnablePreface() { buffer_.EnablePreface(); }

  void DisablePreface() { buffer_.DisablePreface(); }

 private:
  LogBuffer buffer_;
  // at what level of detail output messages
  Log::Level ReportLevel_ = Log::error;

  // if true, only a single processor job is executed
  bool maverick_ = false;

  std::string Messages() { return buffer_.Messages(); }
};

/**
 *   \brief Timestamp returns the current time as a string
 *  Example: cout << TimeStamp()
 */
class TimeStamp {
 public:
  friend std::ostream &operator<<(std::ostream &os, const TimeStamp &) {
    std::time_t now_time =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm *timeinfo = std::localtime(&now_time);
    os << timeinfo->tm_year + 1900 << "-" << timeinfo->tm_mon + 1 << "-"
       << timeinfo->tm_mday << " " << timeinfo->tm_hour << ":"
       << timeinfo->tm_min << ":" << timeinfo->tm_sec;
    return os;
  }

  explicit TimeStamp() = default;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LOGGER_H
