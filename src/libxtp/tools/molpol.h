/* 
 *            Copyright 2009-2016 The VOTCA Development Team
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

#ifndef MOLPOLTOOL_H
#define MOLPOLTOOL_H


#include <votca/xtp/qmtool.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/molpolengine.h>


namespace votca { namespace xtp {


class MolPolTool : public QMTool
{
public:

    MolPolTool() { };
   ~MolPolTool() { };

    string Identify() { return "molpol"; }

    void   Initialize(Property *options);
    bool   Evaluate();

private:

    // FILE I/O OPTIONS
    string          _mps_input;
    string          _mps_output;
    string          _pol_output;
    // MOLPOL OPTIMIZATION
    bool            _do_optimize;
    matrix          _target;
    double          _tolerance;
    vector<string>  _scaling_pattern;
    // SCF MOLPOL OPTIONS
    double          _aDamp;
    double          _wSOR;
    int             _maxIter;
    double          _epsTol;
    // MOLPOL ENGINE
    MolPolEngine    _molpolengine;

};


}}

#endif