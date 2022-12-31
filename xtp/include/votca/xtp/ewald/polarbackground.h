#ifndef VOTCA_XTP_POLARBACKGROUND_H
#define VOTCA_XTP_POLARBACKGROUND_H

#include "votca/tools/property.h"
#include "votca/xtp/backgroundregion.h"

#include <votca/xtp/ewald/ewaldactor.h>
#include <votca/xtp/ewald/threadforce.h>
#include <votca/xtp/ewald/xinteractor.h>
#include <votca/xtp/logger.h>

namespace votca {
namespace xtp {
namespace EWD {

class PolarBackground {
 public:
  // PolarBackground() : /*_top(NULL), _ptop(NULL),*/ _log(NULL),
  //_n_threads(1){};
  // PolarBackground(Topology *top, PolarTop *ptop,
  //                                tools::Property opt, Logger *log)
  //   : _top(top), _ptop(ptop), _log(log), _n_threads(1){};
  // PolarBackground(Topology *, BackgroundRegion &, tools::Property, Logger *);
  //: _top(top), _BGN(BGN), _log(log), _n_threads(1){};*/

  PolarBackground() : /*_top(NULL), _ptop(NULL),*/ _log(NULL), _n_threads(1){};
  PolarBackground(Topology *, PolarTop *, tools::Property, Logger *);
  ~PolarBackground();

  //~PolarBackground();

  // void Initialize();
  void Threaded(int n_threads) { _n_threads = n_threads; }
  void Polarize(int n_threads);
  void Checkpoint(int iter, bool converged);
  bool HasConverged() { return _converged; }

  void FX_RealSpace(std::string mode, bool do_setup_nbs);
  void FX_ReciprocalSpace(std::string S_mode, std::string F_mode,
                          bool gen_kvecs);
  void GenerateKVectors(std::vector<PolarSeg *> &ps1,
                        std::vector<PolarSeg *> &ps2);

  // TODO Would be nicer to have RThread / KThread as a local type in ...
  //      ... but you cannot use local types as template parameters
  //      (which is however possible by compiling with --std=c++0x)

  class RThread : public MultiModeThread<votca::tools::Thread,
                                         EWD::PolarBackground::RThread> {
   public:
    RThread(PolarBackground *master, bool do_setup_nbs) {
      _master = master;
      _do_setup_nbs = do_setup_nbs;
      _not_converged_count = 0;
      //_ewdactor = EwdInteractor(_master->_alpha, _master->_polar_aDamp);
      _ewdactor.Init(_master->_alpha, _master->_polar_aDamp);
      _actor = XInteractor(NULL, _master->_polar_aDamp);

      RegisterStart("FP_MODE", &RThread::FP_FieldCalc);
      RegisterStart("FU_MODE", &RThread::FU_FieldCalc);
      RegisterReset("FP_MODE", &RThread::FX_FieldReset);
      RegisterReset("FU_MODE", &RThread::FX_FieldReset);
      RegisterWload("FP_MODE", &RThread::FX_FieldWload);
      RegisterWload("FU_MODE", &RThread::FX_FieldWload);
    }
    ~RThread() {
      _full_bg_P.clear();
      _part_bg_P.clear();
    }

    RThread *Clone() { return new RThread(_master, _do_setup_nbs); }

    // Input
    void AddSharedInput(std::vector<PolarSeg *> &fullbg) {
      _full_bg_P = fullbg;
    }
    void AddAtomicInput(PolarSeg *pseg) { _part_bg_P.push_back(pseg); }
    // Mode targets
    void FU_FieldCalc();
    void FP_FieldCalc();
    void FX_FieldReset() {
      _not_converged_count = 0;
      _part_bg_P.clear();
    }
    double FX_FieldWload() {
      return 100. * double(_part_bg_P.size()) / double(_full_bg_P.size());
    }
    // Convergence & output related
    const int &NotConverged() { return _not_converged_count; }
    double AvgRco() { return _avg_R_co; }

   private:
    PolarBackground *_master;
    EwdInteractor _ewdactor;  // Long-range tensors (default)
    XInteractor _actor;       // Cut-off tensors (option)
    // Shared thread data
    std::vector<PolarSeg *> _full_bg_P;
    // Atomized thread data
    std::vector<PolarSeg *> _part_bg_P;
    // Convergence & output-related
    int _not_converged_count;
    double _avg_R_co;
    bool _do_setup_nbs;
    const double int2V_m =
        1 / (4 * Pi * 8.854187817e-12) * 1.602176487e-19 / 1.000e-18;
  };

  class KThread : public MultiModeThread<votca::tools::Thread,
                                         EWD::PolarBackground::KThread> {
   public:
    KThread(PolarBackground *master) {
      _master = master;
      //_ewdactor = EwdInteractor(_master->_alpha, _master->_polar_aDamp);
      _ewdactor.Init(_master->_alpha, _master->_polar_aDamp);
      RegisterStart("SP_MODE", &KThread::SP_SFactorCalc);
      RegisterStart("FP_MODE", &KThread::FP_KFieldCalc);
      RegisterStart("SU_MODE", &KThread::SU_SFactorCalc);
      RegisterStart("FU_MODE", &KThread::FU_KFieldCalc);
      RegisterReset("SP_MODE", &KThread::SFactorReset);
      RegisterReset("FP_MODE", &KThread::KFieldReset);
      RegisterReset("SU_MODE", &KThread::SFactorReset);
      RegisterReset("FU_MODE", &KThread::KFieldReset);
      RegisterWload("SP_MODE", &KThread::SFactorWload);
      RegisterWload("FP_MODE", &KThread::KFieldWload);
      RegisterWload("SU_MODE", &KThread::SFactorWload);
      RegisterWload("FU_MODE", &KThread::KFieldWload);
    }
    ~KThread() {
      _full_bg_P.clear();
      _full_kvecs.clear();
      _part_bg_P.clear();
      _part_kvecs.clear();
    }

    KThread *Clone() { return new KThread(_master); }

    void AddSharedInput(std::vector<EWD::KVector *> &full_kvecs) {
      _full_kvecs.clear();
      _full_kvecs = full_kvecs;
    }
    void AddSharedInput(std::vector<PolarSeg *> &full_bg_P) {
      _full_bg_P.clear();
      _full_bg_P = full_bg_P;
    }

    // MODE 1 : Compute structure factor of _part_kvecs using _full_bg_P
    void SP_SFactorCalc();
    void SU_SFactorCalc();
    void SFactorReset() {
      _part_kvecs.clear();
      _full_kvecs.clear();
    }
    void AddAtomicInput(EWD::KVector *k) { _part_kvecs.push_back(k); }
    double SFactorWload() {
      return 1. * double(_part_kvecs.size()) / double(_full_kvecs.size());
    }

    // MODE 2 : Increment fields of _part_bg_P using _full_kvecs
    void FP_KFieldCalc();
    void FU_KFieldCalc();
    void KFieldReset() { ; }
    void AddAtomicInput(PolarSeg *pseg) { _part_bg_P.push_back(pseg); }
    double KFieldWload() {
      return 1. * double(_part_bg_P.size()) / double(_full_bg_P.size());
    }

   private:
    PolarBackground *_master;
    EwdInteractor _ewdactor;
    // Shared thread data
    std::vector<PolarSeg *> _full_bg_P;
    std::vector<EWD::KVector *> _full_kvecs;
    // Atomized thread data
    std::vector<EWD::KVector *> _part_kvecs;
    std::vector<PolarSeg *> _part_bg_P;

   public:
    // Convergence info
    double _rms_sum_re;
    double _sum_im;
  };

 private:
  EwdInteractor _ewdactor;  // Long-range scheme (default)
  XInteractor _actor;       // Cut-off scheme (option)

  // CHECKPOINTING & RESTART OPTIONS
  bool _do_checkpointing;
  bool _do_restart;
  int _restart_from_iter;
  int _max_iter;
  bool _converged;

  // PERIODIC BOUNDARY
  Topology *_top;
  std::string _shape;
  bool _do_use_cutoff;

  // POLAR SEGMENTS
  // Part I - Ewald
  PolarTop *_ptop;
  // BackgroundRegion &_BGN;

  std::vector<PolarSeg *> _bg_P;  // Period. density = _bg_N v _fg_N
  // std::vector<PolarSegment> _bg_P;  // Period. density = _bg_N v _fg_N

  bool _do_compensate_net_dipole;
  std::string _dipole_compensation_type;       // "system" or "segment"
  std::string _dipole_compensation_direction;  // "xyz" or "z"

  // Logger and number of threads
  Logger *_log;
  int _n_threads;

  // CONVERGENCE
  // Part I - Ewald
  double _alpha;  // _a = 1/(sqrt(2)*sigma)
  double _kfactor;
  double _rfactor;
  double _K_co;     // k-space c/o
  double _R_co;     // r-space c/o
  double _crit_dE;  // Energy convergence criterion [eV]
  // bool   _converged_R;            // Did R-space sum converge?
  // bool   _converged_K;            // Did K-space sum converge?
  // bool   _field_converged_R;
  bool _field_converged_K;
  // Part II - Thole
  double _polar_aDamp;
  double _polar_wSOR_N;
  double _polar_cutoff;

  // LATTICE (REAL, RECIPROCAL)
  vec _a;
  vec _b;
  vec _c;                         // Real-space lattice vectors
  int _na_max, _nb_max, _nc_max;  // Max. cell indices to sum over (R)
  vec _A;
  vec _B;
  vec _C;                         // Reciprocal-space lattice vectors
  int _NA_max, _NB_max, _NC_max;  // Max. cell indices to sum over (K)
  double _LxLy;                   // |a^b|
  double _LxLyLz;                 // a*|b^c|

  EWD::VectorSort<EWD::KNorm, EWD::KVector> _kvecsort;
  std::vector<EWD::KVector *> _kvecs_2_0;  // K-vectors with two components =
                                           // zero
  std::vector<EWD::KVector *> _kvecs_1_0;  // K-vectors with one component  =
                                           // zero
  std::vector<EWD::KVector *> _kvecs_0_0;  // K-vectors with no  components =
                                           // zero
  double _kxyz_s1s2_norm;
  const double int2V_m =
      1 / (4 * Pi * 8.854187817e-12) * 1.602176487e-19 / 1.000e-18;
};

}  // namespace EWD
}  // namespace xtp
}  // namespace votca

#endif
