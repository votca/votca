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

#include <votca/xtp/charges.h>

namespace votca { namespace xtp {

int multipoles::read_crg_eps(const char * namefile ){
        ifstream in(namefile);
        string word;
	double mpl;

        int count = 0;
        while( in  ) {
                count ++;
                in >> word;
                if (count %3 == 0){
                    sscanf(word.c_str(), "%lf", &mpl);
                    mpls.push_back(mpl);
                }
        }
        in.close();
        return 0;
}

}}
