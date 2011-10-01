/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _QMCALCULATOR_H
#define	_QMCALCULATOR_H

#include "qmtopology.h"

namespace votca { namespace ctp {

/// the idea of this class is to make QMApplications more flexible

class QMCalculator{
public:
    QMCalculator() {};
    virtual ~QMCalculator() {};

    //virtual const char *Description() = 0;
    const char *Description( const char *name ) {
        // loading the documentation xml file from VOTCASHARE
        char *votca_share = getenv("VOTCASHARE");
        if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
        string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/")+name+string(".xml");
        try {
            Property options;
            load_property_from_xml(options, xmlFile);
            return (options.get(name+string(".description")).as<string>()).c_str();
        } catch(std::exception &error) {
            return (string("XML file or description tag missing: ")+xmlFile).c_str();
        }
    }

    virtual void Initialize(QMTopology *top, Property *options) {}
    virtual bool EvaluateFrame(QMTopology *top) { return true; }
    virtual void EndEvaluate(QMTopology *top) {}
protected:
};

}}

#endif	/* _QMCALCULATOR_H */

