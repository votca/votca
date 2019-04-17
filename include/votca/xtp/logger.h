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

#ifndef VOTCA_XTP_LOG_H
#define VOTCA_XTP_LOG_H

#include <iostream>
#include <sstream>

namespace votca {
namespace xtp {

enum TLogLevel { logERROR, logWARNING, logINFO, logDEBUG };

/*
 * Macros to use the Logger: XTP_LOG(level,logger) << message
 */
#define XTP_LOG(level, log)                           \
  if (&log != NULL && level > (log).getReportLevel()) \
    ;                                                 \
  else                                                \
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
  LogBuffer()
      : std::stringbuf(),
        _errorPreface(" ERROR   "),
        _warnPreface(" WARNING "),
        _infoPreface("         "),
        _dbgPreface(" DEBUG   "),
        _writePreface(true) {}

  // sets the log level (needed for output)
  void setLogLevel(TLogLevel LogLevel) { _LogLevel = LogLevel; }

  // sets Multithreading (buffering required)
  void setMultithreading(bool maverick) { _maverick = maverick; }

  // sets preface strings for logERROR, logWARNING, ...
  void setPreface(TLogLevel level, std::string preface) {
    switch (level) {

      case logERROR:
        _errorPreface = preface;
        break;
      case logWARNING:
        _warnPreface = preface;
        break;
      case logINFO:
        _infoPreface = preface;
        break;
      case logDEBUG:
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
  TLogLevel _LogLevel = TLogLevel::logDEBUG;
  ;

  // temporary buffer to store messages
  std::ostringstream _stringStream;

  // Multithreading
  bool _maverick;

  std::string _errorPreface;
  std::string _warnPreface;
  std::string _infoPreface;
  std::string _dbgPreface;
  std::string _timePreface;
  bool _writePreface;

 protected:
  virtual int sync() {

    std::ostringstream _message;

    if (_writePreface) {
      switch (_LogLevel) {
        case logERROR:
          _message << _errorPreface;
          break;
        case logWARNING:
          _message << _warnPreface;
          break;
        case logINFO:
          _message << _infoPreface;
          break;
        case logDEBUG:
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
 *  Logger* log = new Logger(); // create a logger object
 *  log->setReportLevel(logDEBUG); // output only log messages starting from a
 * DEBUG level XTP_LOG(logERROR,*log) << "Error detected" << flush; // write to
 * the logger at an ERROR level cout << log; // output logger content to
 * standard output \endcode
 *
 *  Logger has four predefined log levels: logERROR, logWARNING, logINFO,
 * logDEBUG.
 */
class Logger : public std::ostream {

  friend std::ostream &operator<<(std::ostream &out, Logger &logger) {
    out << logger.Messages();
    return out;
  }

 public:
  Logger(TLogLevel ReportLevel = logWARNING) : std::ostream(new LogBuffer()) {
    _ReportLevel = ReportLevel;
    _maverick = false;
  }

  ~Logger() {
    // dynamic_cast<LogBuffer *>( rdbuf())->FlushBuffer();
    delete rdbuf();
    rdbuf(NULL);
  }

  Logger &operator()(TLogLevel LogLevel) {
    // rdbuf()->pubsync();
    dynamic_cast<LogBuffer *>(rdbuf())->setLogLevel(LogLevel);
    return *this;
  }

  void setReportLevel(TLogLevel ReportLevel) { _ReportLevel = ReportLevel; }
  void setMultithreading(bool maverick) {
    _maverick = maverick;
    dynamic_cast<LogBuffer *>(rdbuf())->setMultithreading(_maverick);
  }
  bool isMaverick() { return _maverick; }

  TLogLevel getReportLevel() { return _ReportLevel; }

  void setPreface(TLogLevel level, std::string preface) {
    dynamic_cast<LogBuffer *>(rdbuf())->setPreface(level, preface);
  }

  void EnablePreface() { dynamic_cast<LogBuffer *>(rdbuf())->EnablePreface(); }

  void DisablePreface() {
    dynamic_cast<LogBuffer *>(rdbuf())->DisablePreface();
  }

 private:
  // at what level of detail output messages
  TLogLevel _ReportLevel = TLogLevel::logDEBUG;

  // if true, only a single processor job is executed
  bool _maverick;

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
  friend std::ostream &operator<<(std::ostream &os, const TimeStamp &ts) {
    time_t rawtime;
    tm *timeinfo;
    time(&rawtime);
    timeinfo = localtime(&rawtime);
    os << (timeinfo->tm_year) + 1900 << "-" << timeinfo->tm_mon + 1 << "-"
       << timeinfo->tm_mday << " " << timeinfo->tm_hour << ":"
       << timeinfo->tm_min << ":" << timeinfo->tm_sec;
    return os;
  }

  explicit TimeStamp(){};
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_LOG_H
