/*
 *            Copyright 2009-2018 The VOTCA Development Team
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
/// For an earlier history see ctp repo commit 77795ea591b29e664153f9404c8655ba28dc14e9

#ifndef VOTCA_XTP_QMTHREAD_H
#define	VOTCA_XTP_QMTHREAD_H

#include "votca/tools/thread.h"
#include "votca/tools/globals.h"
#include "votca/xtp/logger.h"
#include <iostream>
#include <string>
#include <ctime>

namespace TOOLS = votca::tools;

namespace votca { namespace xtp {

  // ++++++++++++++++++++++++++++++++++++++ //
  // Thread class with local string stream //
  // ++++++++++++++++++++++++++++++++++++++ //

  class QMThread : public TOOLS::Thread
  {
    public:
      QMThread() { _maverick = false; }
      QMThread(bool maverick) { _maverick = maverick; }; 
      ~QMThread() {};

      int  getId() { return _id; }
      void setId(int id) { _id = id; }
      bool isMaverick() { return _maverick; }

      Logger* getLogger() { return &_logger; }
      virtual void Run(void) { ; }

    protected:

      int              _id;
      std::stringstream     _ss;
      bool             _maverick;
      Logger           _logger;

  };


}}

#endif // VOTCA_XTP_QMTHREAD_H
