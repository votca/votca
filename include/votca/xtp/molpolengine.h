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

#ifndef VOTCA_XTP_MOLPOLENGINE_H
#define VOTCA_XTP_MOLPOLENGINE_H

#include <votca/xtp/polartop.h>

namespace votca {
namespace xtp {
    
    
class MolPolEngine
{
public:
    
    MolPolEngine() 
        : _aDamp(0.390), _wSOR(0.30), _maxIter(1024), _epsTol(0.0001)
        { _actor.SetADamp(_aDamp); }
    MolPolEngine(double aDamp, double wSOR, int maxIter, double epsTol)
        : _aDamp(aDamp), _wSOR(wSOR), _maxIter(maxIter), _epsTol(epsTol)
        { _actor.SetADamp(_aDamp); }
   ~MolPolEngine() {}
    
    matrix CalculateMolPol(std::vector<APolarSite*> &poles, bool verbose = true);
    int    SCF_Induce(std::vector<APolarSite*> &poles);
    
private:
    
    BasicInteractor _actor;
    double _aDamp;
    double _wSOR;
    int _maxIter;
    double _epsTol;
};

}}

#endif // VOTCA_XTP_MOLPOLENGINE_H
