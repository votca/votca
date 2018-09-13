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

#ifndef VOTCA_XTP_POLARFRAG_H
#define VOTCA_XTP_POLARFRAG_H

#include <votca/xtp/apolarsite.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/base_object.hpp>

namespace votca { namespace xtp {
    
    
// Forward declarations
class PolarSeg;

class PolarFrag : public std::vector<APolarSite*>
{    
public:    
    
   ~PolarFrag();
    int getId() const { return _id; }
    const std::string &getName() const { return _name; }
    const votca::tools::vec &CalcPosCenterOfGeom();
    const votca::tools::vec &CalcPosPolarWeights();
    const votca::tools::vec &getPos() { return _pos; }
    
    void GeneratePermInduCgSite(bool do_cg_polarizabilities);
    APolarSite *getInduCgSite() { return _indu_cg_site; }
    APolarSite *getPermCgSite() { return _perm_cg_site; }
    
    // Serialization interface
    template<class Archive>
    void serialize(Archive &arch, const unsigned int version) {
        arch & boost::serialization::base_object< std::vector<APolarSite*> >(*this);
        arch & _id;
        arch & _name;
        arch & _pseg;
        arch & _pos;
        return;
    }
    
    friend class PolarSeg;
    friend class boost::serialization::access;
    
    
private:
    
    PolarFrag() 
        : _id(-1), _name(""), _pseg(NULL), _pos(vec(0,0,0)),
          _indu_cg_site(NULL), _perm_cg_site(NULL) {}
    PolarFrag(PolarSeg *pseg, int id, std::string name)
        : _id(id), _name(name), _pseg(pseg), _pos(vec(0,0,0)),
          _indu_cg_site(NULL), _perm_cg_site(NULL) {}
    
    int _id;
    std::string _name;
    PolarSeg *_pseg;
    votca::tools::vec _pos;
    
    APolarSite *_indu_cg_site;
    APolarSite *_perm_cg_site;
};

}}

#endif // VOTCA_XTP_POLARFRAG_H
