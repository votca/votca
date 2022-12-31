#ifndef VOTCA_CTP_EWALD3D_H
#define VOTCA_CTP_EWALD3D_H

#include <votca/csg/boundarycondition.h>
#include <votca/xtp/ewald/ewaldnd.h>
#include <votca/xtp/ewald/polartop.h>
#include <votca/xtp/ewald/xjob.h>
#include <votca/xtp/qmthread.h>

namespace CSG = votca::csg;

namespace votca {
namespace xtp {

// NOTE: This is not a conventional 3D Ewald summation, so use carefully
//       (tuned for the purpose of site-energy calculations)
// NOTE: PolarTop should be set-up with three containers: FGC, FGN, BGN
//       MGN is set-up in constructor using the real-space c/o (from input)
//       The topology is used to retrieve information on the sim. box (PB)
//       All polar segments should be positioned as nearest images of the
//       foreground charge density (FGC, FGN).
//       All polar segments should be appropriately charged (Q00).
// NOTE: The k-shell grouping algorithm can fail for strongly skewed boxes.
//

class Ewald3D3D : public Ewald3DnD {

 public:
  Ewald3D3D(const Topology *top, PolarTop *ptop, tools::Property *opt,
            Logger *log);
  ~Ewald3D3D();

  std::string IdentifyMethod() { return "3D x 3D"; }
  EWD::triple<> ConvergeReciprocalSpaceSum(std::vector<PolarSeg *> &target);
  EWD::triple<> CalculateShapeCorrection(std::vector<PolarSeg *> &target);
  double CalculateSq2(vec &k);

 private:
};

}  // namespace xtp
}  // namespace votca

#endif