#ifndef _CTP_QMTHREAD_H
#define _CTP_QMTHREAD_H

#include "votca/tools/globals.h"
#include "votca/tools/thread.h"
#include "votca/xtp/logger.h"
#include <ctime>
#include <iostream>
#include <string>

namespace votca {
namespace xtp {

// ++++++++++++++++++++++++++++++++++++++ //
// Thread class with local string stream //
// ++++++++++++++++++++++++++++++++++++++ //

class QMThread : public tools::Thread {

  /*
  friend ostream& operator<<( ostream& out, QMThread&  t ) {
     out << (t._ss).str();
     t._ss.str( "" );
     return out;
 }


  template <class Traits>
  friend QMThread& operator<<( QMThread &t,  Traits& inp ) {

      if ( tools::globals::verbose ) {
          if ( t._maverick ) { std::cout << inp ; }
          else { t._ss << inp ; }
      }
      return t;
  }

  template <class Traits>
  friend QMThread& operator>>( Traits& inp, QMThread & t ) {

      if ( tools::globals::verbose ) {
          if ( t._maverick ) { std::cout << inp ; }
          else { t._ss << inp ; }
      }
      return t;
  }
  */

 public:
  QMThread() { _maverick = false; }
  QMThread(bool maverick) { _maverick = maverick; };
  ~QMThread(){};

  int getId() { return _id; }
  void setId(int id) { _id = id; }
  bool isMaverick() { return _maverick; }

  Logger* getLogger() { return &_logger; }
  virtual void Run(void) { ; }

 protected:
  int _id;
  std::stringstream _ss;
  bool _maverick;
  Logger _logger;
};

}  // namespace xtp
}  // namespace votca

#endif /* _CTP_QMTHREAD_H */
