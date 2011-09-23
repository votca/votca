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

#ifndef __VOTCA_MD2QM_CUSTOMFIELDS_H
#define	__VOTCA_MD2QM_CUSTOMFIELDS_H

#include <map>
#include <string>

namespace votca { namespace ctp {

class CustomFields
{
public:
    void setDouble(std::string key, double value) { _double_values[key] = value; }
    double getDouble(std::string key) { return _double_values[key]; }
    std::map<std::string, double> &DoubleValues() { return _double_values; }
    bool DoubleExists(std::string key) {
   return _double_values.find(key) != _double_values.end();
} 
protected:
    std::map<std::string, double> _double_values;
};

}}

#endif	/* CUSTOMFIELDS_H */

