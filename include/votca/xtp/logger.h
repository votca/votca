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

#pragma once
#ifndef VOTCA_XTP_LOG_H
#define VOTCA_XTP_LOG_H

#include <chrono>
#include <iostream>
#include <sstream>
#include <votca/tools/globals.h>
namespace votca {
namespace xtp {

/*
 * Macros to use the Logger: XTP_LOG(level,logger) << message
 */
#define XTP_LOG(level, log)                              \
  if (&log != nullptr && level > (log).getReportLevel()) \
    ;                                                    \
  else                                                   \
    (log)(level)

#define XTP_LOG_SAVE(level, log)      \
  if (level > (log).getReportLevel()) \
    ;                                 \
  else                                \
    (log)(level)

/*
 * Custom buffer to store messages
 */
class LogBuffer : public std::stringbuf {

 public:
  LogBuffer() : std::stringbuf() { _LogLevel = Log::current_level; }

  // sets the log level (needed for output)
  void setLogLevel(Log::Level LogLevel) { _LogLevel = LogLevel; }

  // sets Multithreading (buffering required)
  void setMultithreading(bool maverick) { _maverick = maverick; }

  // sets preface strings for Log::error, Log::warning, ...
  void setPreface(Log::Level level, std::string preface) {
    switch (level) {

      case Log::Level::error:
        _errorPreface = preface;
        break;
      case Log::Level::warning:
        _warnPreface = preface;
        break;
      case Log::Level::info:
        _infoPreface = preface;
        break;
      case Log::Level::debug:
        _dbgPreface = preface;
        break;
    }
  }

  void EnablePreface() { _writePreface = true; }
  void DisablePreface() { _writePreface = false; }

  // flushes all collected messages
  void FlushBuffer() {
    std::cout << _stringStream.str();
    _stringStream.str("");
  }

  // returns the pointer to the collected messages
  std::string Messages() {
    std::string _messages = _stringStream.str();
    _stringStream.str("");
    return _messages;
  }

 private:
  // Log Level (WARNING, INFO, etc)
  Log::Level _LogLevel = Log::Level::error;

  // temporary buffer to store messages
  std::ostringstream _stringStream;

  // Multithreading
  bool _maverick;

  std::string _errorPreface = "\n ERROR   ";
  std::string _warnPreface = "\n WARNING ";
  std::string _infoPreface = "\n         ";
  std::string _dbgPreface = "\n DEBUG   ";
  bool _writePreface = true;

 protected:
  int sync() override {

    std::ostringstream _message;

    if (_writePreface) {
      switch (_LogLevel) {
        case Log::Level::error:
          _message << _errorPreface;
          break;
        case Log::Level::warning:
          _message << _warnPreface;
          break;
        case Log::Level::info:
          _message << _infoPreface;
          break;
        case Log::Level::debug:
          _message << _dbgPreface;
          break;
      }
    }

    if (!_maverick) {
      // collect all messages of one thread
      _stringStream << _message.str() << " " << str();
    } else {
      // if only one thread outputs, flush immediately
      std::cout << _message.str() << " " << str() << std::flush;
    }
    _message.str("");
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
 *  #include <votca/xtp/logger.h>
 *  Logger log; // create a logger object
 *  log.setReportLevel(Log::error); // output only log messages starting from a
 *  level XTP_LOG(Log::error,*log) << "Error detected" << flush; // write to
 * the logger at an ERROR level cout << log; // output logger content to
 * standard output \endcode
 *
 *  Logger has four predefined log levels: error, warning, info,
 * debug.
 */
class Logger : public std::ostream {

  friend std::ostream &operator<<(std::ostream &log_out, Logger &logger) {
    log_out << logger.Messages();
    return log_out;
  }

 public:
  Logger() : std::ostream(new LogBuffer()), _ReportLevel(Log::current_level){};
  Logger(Log::Level ReportLevel)
      : std::ostream(new LogBuffer()), _ReportLevel(ReportLevel) {}

  ~Logger() final {
    delete rdbuf();
    rdbuf(nullptr);
  }

  Logger &operator()(Log::Level LogLevel) {
    dynamic_cast<LogBuffer *>(rdbuf())->setLogLevel(LogLevel);
    return *this;
  }

  void setReportLevel(Log::Level ReportLevel) { _ReportLevel = ReportLevel; }
  void setMultithreading(bool maverick) {
    _maverick = maverick;
    dynamic_cast<LogBuffer *>(rdbuf())->setMultithreading(_maverick);
  }
  bool isMaverick() const { return _maverick; }

  Log::Level getReportLevel() const { return _ReportLevel; }

  void setPreface(Log::Level level, const std::string &preface) {
    dynamic_cast<LogBuffer *>(rdbuf())->setPreface(level, preface);
  }

  void setCommonPreface(const std::string &preface) {
    setPreface(Log::info, preface);
    setPreface(Log::error, preface);
    setPreface(Log::debug, preface);
    setPreface(Log::warning, preface);
  }

  void EnablePreface() { dynamic_cast<LogBuffer *>(rdbuf())->EnablePreface(); }

  void DisablePreface() {
    dynamic_cast<LogBuffer *>(rdbuf())->DisablePreface();
  }

 private:
  // at what level of detail output messages
  Log::Level _ReportLevel = Log::error;

  // if true, only a single processor job is executed
  bool _maverick = false;

  std::string Messages() {
    return dynamic_cast<LogBuffer *>(rdbuf())->Messages();
  }
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

#endif  // VOTCA_XTP_LOG_H
