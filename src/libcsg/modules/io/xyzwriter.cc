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

#include <stdio.h>
#include <string>

#include <votca/csg/xyzwriter.h>

namespace votca { namespace csg {

using namespace std;

void XYZWriter::Open(string file, bool bAppend)
{
    if(bAppend){
    _out.open(file, std::ios_base::app);
    }else{
    _out.open(file);
    }
}

void XYZWriter::Close()
{
    _out.close();
}


void XYZWriter::Write(Topology *conf)
{
    std::string header=(boost::format("frame: %1$d time: %2$f\n")
                         % (conf->getStep()+1) %conf->getTime()).str();
    Write<Topology>(*conf,header);  
}






}}
