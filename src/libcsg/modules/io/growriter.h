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
// 
// File:   growriter.h
// Author: victor
//
// Created on 9. January 2008, 18:11
//

#ifndef _GROWRITER_H
#define	_GROWRITER_H

#include <stdio.h>
#include "topology.h"
#include "trajectorywriter.h"

using namespace std;

class GROWriter
: public TrajectoryWriter
{
public:
    
    void Open(string file, bool bAppend = false);
    void Close();
    
    void Write(Topology *conf);

    private:
        FILE *_out;
};

#endif	/* _GROWRITER_H */

