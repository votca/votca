#ifndef VOTCA_XTP_PEWALD3D_H
#define VOTCA_XTP_PEWALD3D_H

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

class PEwald3D3D : public Ewald3DnD {

 public:
  PEwald3D3D(const Topology *top, PolarTop *ptop, tools::Property *opt,
             Logger *log);
  ~PEwald3D3D();

  std::string IdentifyMethod() { return "Polar 3D x 3D"; }
  void GenerateKVectors(std::vector<PolarSeg *> &ps1,
                        std::vector<PolarSeg *> &ps2);

  EWD::triple<> ConvergeRealSpaceSum(std::vector<PolarSeg *> &target);
  EWD::triple<> ConvergeReciprocalSpaceSum(std::vector<PolarSeg *> &target);
  EWD::triple<> CalculateForegroundCorrection(std::vector<PolarSeg *> &target);
  EWD::triple<> CalculateShapeCorrection(std::vector<PolarSeg *> &target);
  EWD::triple<> CalculateHigherRankCorrection(std::vector<PolarSeg *> &target) {
    return EWD::triple<>(0, 0, 0);
  }
  EWD::triple<> CalculateK0Correction(std::vector<PolarSeg *> &target) {
    return EWD::triple<>(0, 0, 0);
  }

  void Field_ConvergeRealSpaceSum();
  void Field_ConvergeReciprocalSpaceSum();
  void Field_CalculateForegroundCorrection();
  void Field_CalculateShapeCorrection();

  void Potential_ConvergeRealSpaceSum(std::vector<PolarSeg *> &target);
  void Potential_ConvergeReciprocalSpaceSum(std::vector<PolarSeg *> &target);
  void Potential_CalculateForegroundCorrection(std::vector<PolarSeg *> &target);
  void Potential_CalculateShapeCorrection(std::vector<PolarSeg *> &target);

  void ScanCutoff();

 private:
  const double int2eV =
      1 / (4 * EWD::Pi * 8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
  const double int2V_m =
      1 / (4 * EWD::Pi * 8.854187817e-12) * 1.602176487e-19 / 1.000e-18;
};

struct TinyNeighbour {
  TinyNeighbour(PolarSeg *seg, vec L) : _nb(seg), _L(L) { ; }
  PolarSeg *_nb;
  vec _L;
};

}  // namespace xtp
}  // namespace votca

#endif