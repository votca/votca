#ifndef VOTCA_CTP_POLARFRAG_H
#define VOTCA_CTP_POLARFRAG_H

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <votca/xtp/ewald/apolarsite.h>

namespace votca {
namespace xtp {

typedef Eigen::Vector3d vec;
typedef Eigen::Matrix3d matrix;

// Forward declarations
class PolarSeg;

class PolarFrag : public std::vector<APolarSite *> {
 public:
  ~PolarFrag();
  int getId() const { return _id; }
  const std::string &getName() const { return _name; }
  const vec &CalcPosCenterOfGeom();
  const vec &CalcPosPolarWeights();
  const vec &getPos() { return _pos; }

  void GeneratePermInduCgSite(bool do_cg_polarizabilities);
  APolarSite *getInduCgSite() { return _indu_cg_site; }
  APolarSite *getPermCgSite() { return _perm_cg_site; }

  // Serialization interface
  template <class Archive>
  void serialize(Archive &arch, const unsigned int version) {
    arch &boost::serialization::base_object<vector<APolarSite *> >(*this);
    arch &_id;
    arch &_name;
    arch &_pseg;
    arch &_pos;
    if (version > 5) {
      return;
    }  // to avoid warning at compile time

    return;
  }

  friend class PolarSeg;
  friend class boost::serialization::access;

 private:
  PolarFrag()
      : _id(-1),
        _name(""),
        _pseg(NULL),
        _pos(vec(0, 0, 0)),
        _indu_cg_site(NULL),
        _perm_cg_site(NULL) {}
  PolarFrag(PolarSeg *pseg, int id, std::string name)
      : _id(id),
        _name(name),
        _pseg(pseg),
        _pos(vec(0, 0, 0)),
        _indu_cg_site(NULL),
        _perm_cg_site(NULL) {}

  int _id;
  std::string _name;
  PolarSeg *_pseg;
  vec _pos;

  APolarSite *_indu_cg_site;
  APolarSite *_perm_cg_site;
};

}  // namespace xtp
}  // namespace votca

#endif