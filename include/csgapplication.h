/*
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_CSGAPPLICATION_H
#define	__VOTCA_CSGAPPLICATION_H

#include <votca/tools/application.h>
#include "topology.h"

namespace votca { namespace csg {
using namespace votca::tools;

class CsgApplication
    : public Application
{
public:
    CsgApplication();
    ~CsgApplication();

    void Initialize();
    bool EvaluateOptions();
    void Run(void);

    void ShowHelpText(std::ostream &out);

    /// \brief overload and return true to enable mapping command line options
    virtual bool DoMapping(void) { return false; }
    /// \brief overload and return true to enable trajectory command line options
    virtual bool DoTrajectory(void) { return false; }

    /// \brief called after topology was loaded
    virtual void EvaluateTopology(Topology *top) {}

    /// \brief called before the first frame
    virtual void BeginEvaluate(Topology *top, Topology *top_ref = 0);
    /// \brief called after the last frame
    virtual void EndEvaluate();
    // \brief called for each frame which is mapped
    virtual void EvalConfiguration(Topology *top, Topology *top_ref = 0);
};

}}

#endif	/* __VOTCA_CSGAPPLICATION_H */

