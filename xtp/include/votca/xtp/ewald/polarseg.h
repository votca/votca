#ifndef __POLARSEG__H
#define __POLARSEG__H

// #include <votca/tools/vec.h>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <votca/xtp/ewald/apolarsite.h>
#include <votca/xtp/ewald/polarfrag.h>
#include <votca/xtp/vxc_grid.h>

namespace votca {
namespace xtp {

typedef Eigen::Vector3d vec;
typedef Eigen::Matrix3d matrix;

class PolarNb;

class PolarSeg : public std::vector<APolarSite *> {

 public:
  PolarSeg()
      : _id(-1),
        _pos(vec(0, 0, 0)),
        _is_charged(true),
        _is_polarizable(true),
        _indu_cg_site(NULL),
        _perm_cg_site(NULL), _grid() {}
  PolarSeg(int id, std::vector<APolarSite *> &psites);
  PolarSeg(PolarSeg *templ, bool do_depolarize);
  explicit PolarSeg(int id)
      : _id(id),
        _pos(vec(0, 0, 0)),
        _is_charged(true),
        _is_polarizable(true),
        _indu_cg_site(NULL),
        _perm_cg_site(NULL), _grid() {}
  ~PolarSeg();

  const int &getId() { return _id; }
  const vec &getPos() { return _pos; }
  void setId(int id) { _id = id; }

  // Polar fragments
  PolarFrag *AddFragment(std::string name);
  vector<PolarFrag *> &PolarFrags() { return _pfrags; }
  // Local neighbor-list
  vector<PolarNb *> &PolarNbs() { return _nbs; }
  void ReservePolarNbs(int nbsize) { _nbs.reserve(nbsize); }
  PolarNb *AddNewPolarNb(PolarSeg *pseg);
  void AddPolarNb(PolarNb *nb) { _nbs.push_back(nb); }
  void ClearPolarNbs();
  // Position & total charge
  void Translate(const vec &shift);
  void CalcPos();
  double CalcTotQ();
  vec CalcTotD();
  void Coarsegrain(bool cg_anisotropic);
  void GeneratePermInduCgSite(bool do_cg_polarizabilities);
  APolarSite *getInduCgSite() { return _indu_cg_site; }
  APolarSite *getPermCgSite() { return _perm_cg_site; }
  // Evaluates to "true" if ANY contained polar site has charge != 0
  void CalcIsCharged();
  bool IsCharged() { return _is_charged; }
  // Evaluates to "true" if ANY contained polar site has polarizability > 0
  void CalcIsPolarizable();
  bool IsPolarizable() { return _is_polarizable; }

  // setting and accessing a numerical integration grid
  void setGrid(const Vxc_Grid& grid) {_grid = grid; _has_grid = true;}
  bool hasGrid() {return _has_grid;}
  Vxc_Grid &getGrid() {return _grid;}

  // File output methods
  void PrintPolarNbPDB(std::string outfile);
  void WriteMPS(std::string mpsfile, std::string tag = "");
  // Serialization interface
  template <class Archive>
  void serialize(Archive &arch, const unsigned int version) {
    arch &boost::serialization::base_object<std::vector<APolarSite *> >(*this);
    arch &_pfrags;
    arch &_id;
    arch &_pos;
    arch &_is_charged;
    arch &_is_polarizable;
    if (version > 5) return;
    return;
  }

 private:
  int _id;
  vec _pos;
  bool _is_charged;
  bool _is_polarizable;
  vector<PolarFrag *> _pfrags;
  vector<PolarNb *> _nbs;
  APolarSite *_indu_cg_site;
  APolarSite *_perm_cg_site;
  // Integration grid for QM-Ewald
  Vxc_Grid _grid;
  bool _has_grid = false;
};

class PolarNb {
 public:
  // s22x = dr(nb,shift;pbc) - dr(nb,shift) + imageboxvector
  // Use s22x to obtain effective connection vector
  // dreff = top->ShortestConnect(ref,nb) + _pbcshift
  PolarNb(PolarSeg *nb, vec &dr12, vec &s22x)
      : _nb(nb), _dr12(dr12), _s22x(s22x){};
  explicit PolarNb(PolarSeg *nb)
      : _nb(nb), _dr12(vec(0, 0, 0)), _s22x(vec(0, 0, 0)){};
  ~PolarNb() { _nb = 0; }
  PolarSeg *getNb() { return _nb; }
  vec &getR() { return _dr12; }
  vec &getS() { return _s22x; }

 private:
  PolarSeg *_nb;
  vec _dr12;
  vec _s22x;
};

}  // namespace xtp
}  // namespace votca

#endif