/*
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_CALCULATOR_H
#define	__VOTCA_KMC_CALCULATOR_H

#include <votca/tools/property.h>

namespace votca { namespace xtp {

    using namespace tools;

class KMCCalculator
{
public:
    
    KMCCalculator() {};
    virtual     ~KMCCalculator() {};
    
    virtual void Initialize(const char *filename, Property *options, const char *outputfile) {};
    virtual bool EvaluateFrame() { return true; }
    virtual void EndEvaluate() {}    
    
    void         setnThreads(int nThreads) { _nThreads = nThreads; }
    
protected:
    
    int _nThreads;
};

}}

#endif	/* __VOTCA_KMC_CALCULATOR_H */

