/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

#ifndef __VOTCA_CTP_LOG_H
#define	__VOTCA_CTP_LOG_H

#include <sstream>
#include <iostream>

namespace votca { namespace ctp {

enum TLogLevel {logERROR, logWARNING, logINFO, logDEBUG};
 
/*
 * Macros to use the Logger (level,logger) << message
 */
#define LOG(level, log) \
if ( &log != NULL && level > (log).getReportLevel() ) ; \
else (log)(level)

/*
 * Custom buffer to store messages
 */
class LogBuffer : public std::stringbuf {

public:
	LogBuffer() : std::stringbuf() {}
        
        // sets the log level (needed for output)
	void setLogLevel(TLogLevel LogLevel) { _LogLevel = LogLevel; }
        
        // sets Multithreading (buffering required)
        void setMultithreading( bool maverick ) { _maverick = maverick; }
        
        // flushes all collected messages
        void FlushBuffer(){ std::cout << _stringStream.str(); _stringStream.str(""); }
        
        // returns the pointer to the collected messages
        std::string Messages() { 
            string _messages = _stringStream.str(); 
            _stringStream.str("");
            return _messages; 
        }

private:
  
  // Log Level (WARNING, INFO, etc) 
  TLogLevel _LogLevel;
  
  // temporary buffer to store messages
  std::ostringstream _stringStream;
  
  // Multithreading
  bool _maverick;

protected:
	virtual int sync() {
            
            std::ostringstream _message;

            switch ( _LogLevel )
            {
                case logERROR: 
                    _message << " ERROR   ";
                    break;
                case logWARNING:
                    _message << " WARNING ";
                    break;      
                case logINFO:
                    _message << " ";
                    break;      
                case logDEBUG:
                    _message << " DEBUG   ";
                    break;      
            }
            
            if ( !_maverick ) {
                // collect all messages of one thread
                _stringStream << _message.str()  << " " << str();
            } else {
                // if only one thread outputs, flash immediately
                std::cout << _message.str() << " " << str() << std::flush;
            }
            _message.str("");
	    str("");
	    return 0;
	}
      
/*        
        // timestamp 
        std::string timestamp() {
                std::ostringstream stream;   
                time_t rawtime;
                tm * timeinfo;
 
                time(&rawtime);
                timeinfo = localtime( &rawtime );
 
                stream  //<< (timeinfo->tm_year)+1900
                        //<< "-" << timeinfo->tm_mon + 1
                        //<< "-" << timeinfo->tm_mday 
                        << " " << timeinfo->tm_hour
                        << ":" << timeinfo->tm_min 
                        << ":"  << timeinfo->tm_sec
                        ;
                return stream.str();  
        }        
 */
};


/**
*   \brief Logger for thread-safe output of messages
*  Logger writes messages into LogBuffer
*  Inheritance from ostream allows to use << and >> for writing    
*/
class Logger : public std::ostream {

       friend std::ostream& operator<<( std::ostream& out, Logger&  logger ) {
           out << logger.Messages();
           return out;
       }
    
public:
	Logger( TLogLevel ReportLevel =  logWARNING) : std::ostream(new LogBuffer()) { 
            _ReportLevel = ReportLevel; 
            _maverick = false;
        }
        
	~Logger() {
            //dynamic_cast<LogBuffer *>( rdbuf())->FlushBuffer();
            delete rdbuf();
            rdbuf(NULL);
	}
        
	Logger &operator()( TLogLevel LogLevel ) {
		//rdbuf()->pubsync();
		dynamic_cast<LogBuffer *>( rdbuf() )->setLogLevel(LogLevel);
		return *this;
	}
        
        void setReportLevel( TLogLevel ReportLevel ) { _ReportLevel = ReportLevel; }
        void setMulithreading( bool maverick ) { 
            _maverick = maverick;
            dynamic_cast<LogBuffer *>( rdbuf() )->setMultithreading( _maverick );
        }
        
        TLogLevel getReportLevel( ) { return _ReportLevel; }
        
private:
    // at what level of detail output messages
    TLogLevel _ReportLevel;
    
    // if true, only a single processor job is executed
    bool      _maverick;
    
    std::string Messages() {
        return dynamic_cast<LogBuffer *>( rdbuf() )->Messages();
    }
        
};

/**
*   \brief Timestamp returns the current time as a string
*  Example: cout << TimeStamp()
*/
class TimeStamp 
{
  public:
    friend std::ostream & operator<<(std::ostream &os, const TimeStamp& ts)
    {
        time_t rawtime;
        tm * timeinfo;
        time(&rawtime);
        timeinfo = localtime( &rawtime );
        os  << (timeinfo->tm_year)+1900
            << "-" << timeinfo->tm_mon + 1
            << "-" << timeinfo->tm_mday 
            << " " << timeinfo->tm_hour
            << ":" << timeinfo->tm_min 
            << ":"  << timeinfo->tm_sec;
         return os;    
    }
    
    explicit TimeStamp() {};
    
};


}}

#endif /* __VOTCA_CTP_LOG_H */