#include "votca/xtp/ewald/polarbackground.h"
#include "votca/tools/property.h"
#include <boost/format.hpp>
#include <votca/tools/globals.h>

namespace votca {
namespace xtp {
namespace EWD {

using boost::format;

PolarBackground::PolarBackground(Topology *top, PolarTop *ptop,
                                 tools::Property opt, Logger *log)
    : _top(top), _ptop(ptop), _log(log), _n_threads(1) {

  // EVALUATE OPTIONS
  std::string pfx = "";
  // Coulomb method
  std::string method = opt.get(pfx + ".coulombmethod.method").as<std::string>();
  if (method == "ewald")
    _do_use_cutoff = false;
  else if (method == "cutoff")
    _do_use_cutoff = true;
  else
    throw std::runtime_error(
        "Invalid parameter in options.ewdbgpol."
        "coulombmethod.method");
  // Dipole compensation options
  if (opt.exists(pfx + ".coulombmethod.dipole_corr"))
    _do_compensate_net_dipole =
        opt.get(pfx + ".coulombmethod.dipole_corr").as<bool>();
  else
    _do_compensate_net_dipole = false;
  if (opt.exists(pfx + ".coulombmethod.dipole_corr_type"))
    _dipole_compensation_type =
        opt.get(pfx + ".coulombmethod.dipole_corr_type").as<std::string>();
  else
    _dipole_compensation_type = "system";
  if (opt.exists(pfx + ".coulombmethod.dipole_corr_direction"))
    _dipole_compensation_direction =
        opt.get(pfx + ".coulombmethod.dipole_corr_direction").as<std::string>();
  else
    _dipole_compensation_direction = "xyz";
  // Ewald parameters
  _shape = opt.get(pfx + ".coulombmethod.shape").as<std::string>();
  _R_co = opt.get(pfx + ".coulombmethod.cutoff").as<double>();
  _crit_dE = opt.get(pfx + ".convergence.energy").as<double>();
  if (opt.exists(pfx + ".convergence.kfactor"))
    _kfactor = opt.get(pfx + ".convergence.kfactor").as<double>();
  else
    _kfactor = 100.;
  if (opt.exists(pfx + ".convergence.rfactor"))
    _rfactor = opt.get(pfx + ".convergence.rfactor").as<double>();
  else
    _rfactor = 6.;
  // Polar parameters
  if (opt.exists(pfx + ".polarmethod.cutoff"))
    _polar_cutoff = opt.get(pfx + ".polarmethod.cutoff").as<double>();
  else
    _polar_cutoff = 0.0;
  if (opt.exists(pfx + ".polarmethod.wSOR_N"))
    _polar_wSOR_N = opt.get(pfx + ".polarmethod.wSOR_N").as<double>();
  else
    _polar_wSOR_N = 0.35;
  if (opt.exists(pfx + ".polarmethod.aDamp"))
    _polar_aDamp = opt.get(pfx + ".polarmethod.aDamp").as<double>();
  else
    _polar_aDamp = 0.390;
  // Checkpointing
  if (opt.exists(pfx + ".control.checkpointing"))
    _do_checkpointing = opt.get(pfx + ".control.checkpointing").as<bool>();
  else
    _do_checkpointing = false;
  if (opt.exists(pfx + ".control.max_iter"))
    _max_iter = opt.get(pfx + ".control.max_iter").as<int>();
  else
    _max_iter = -1;

  // EWALD INTERACTION PARAMETERS (GUESS ONLY)
  _K_co = _kfactor / _R_co;
  _alpha = _rfactor / _R_co;
  //_ewdactor = EwdInteractor(_alpha, _polar_aDamp);
  _ewdactor.Init(_alpha, _polar_aDamp);
  _actor = XInteractor(NULL, _polar_aDamp);

  // SET-UP REAL & RECIPROCAL SPACE
  _a = votca::tools::conv::bohr2nm *
       _top->getBox().col(0);  // NEW: in a_0, needed here: nm
  _b = votca::tools::conv::bohr2nm * _top->getBox().col(1);
  _c = votca::tools::conv::bohr2nm * _top->getBox().col(2);
  _LxLyLz = _a.dot(_b.cross(_c));
  _LxLy = (_a.cross(_b)).norm();

  _A = 2 * M_PI / _LxLyLz * _b.cross(_c);
  _B = 2 * M_PI / _LxLyLz * _c.cross(_a);
  _C = 2 * M_PI / _LxLyLz * _a.cross(_b);

  _na_max = int(ceil(_R_co / _a.cwiseAbs().maxCoeff() - 0.5)) + 1;
  _nb_max = int(ceil(_R_co / _b.cwiseAbs().maxCoeff() - 0.5)) + 1;
  _nc_max = int(ceil(_R_co / _c.cwiseAbs().maxCoeff() - 0.5)) + 1;

  _NA_max = int(ceil(_K_co / _A.cwiseAbs().maxCoeff()));
  _NB_max = int(ceil(_K_co / _B.cwiseAbs().maxCoeff()));
  _NC_max = int(ceil(_K_co / _C.cwiseAbs().maxCoeff()));

  // SET-UP POLAR GROUNDS (FORE-, MID-, BACK-)
  _bg_P.clear();
  _bg_P = ptop->BGN();

  //_bg_P = BGN->segments_;
  // RESTART / CONVERGENCE OPTIONS
  _converged = false;
  _do_restart = false;
  _restart_from_iter = ptop->getPolarizationIter();
  if (_restart_from_iter > -1) {
    XTP_LOG(Log::info, *_log) << "Restarting from iteration "
                              << _restart_from_iter << std::flush << std::flush;
    _do_restart = true;
  }

  // CALCULATE COG POSITIONS, NET CHARGE
  // std::vector<PolarSegment>::iterator sit;
  // std::vector<PolarSite>::iterator pit;
  double Q_bg_P = 0.0;
  int estat_count = 0;
  int polar_count = 0;
  // for (sit = BGN->begin(); sit < BGN->end(); ++sit) {
  for (auto &pseg : _bg_P) {
    pseg->CalcPos();
    pseg->CalcIsCharged();
    pseg->CalcIsPolarizable();
    Q_bg_P += pseg->CalcTotQ();
    if (pseg->IsCharged()) estat_count += 1;
    if (pseg->IsPolarizable()) polar_count += 1;
  }

  XTP_LOG(Log::info, *_log)
      << (format("Net ground charge and size:")).str() << std::flush
      << (format("  o Q(BGP) = %1$+1.3fe |BGP| = %2$+5d") % Q_bg_P %
          _bg_P.size())
      << std::flush
      << (format("  o Activity <qdQ %1$d/%3$d> <P %2$d/%3$d>") % estat_count %
          polar_count % _bg_P.size())
      << std::flush;

  if (std::abs(Q_bg_P) > 1e-2) {
    std::cout << std::endl;
    std::cout << std::endl
              << format(
                     "***************************** ERROR "
                     "******************************");
    std::cout << std::endl
              << format(
                     "       Background charge |Q(BGP)| is larger than 0.01e.");
    std::cout << std::endl
              << format("       Be more precise: e.g. rounding error?");
    std::cout << std::endl
              << format(
                     "       Or think again: e.g. erroneous parametrization?");
    std::cout << std::endl
              << format(
                     "*********************************************************"
                     "*********");
    std::cout << std::endl;
  }

  // APPLY SYSTEM DIPOLE COMPENSATION
  if (_do_compensate_net_dipole) {
    XTP_LOG(Log::info, *_log)
        << (format("  o Dpl. compensation: type '%1$s', "
                   "direction '%2$s' ") %
            _dipole_compensation_type % _dipole_compensation_direction)
        << std::flush;
    vec system_dpl(0, 0, 0);
    int charged_count = 0;
    // for (sit = _BGN->begin(); sit < _BGN->end(); ++sit) {
    for (auto &pseg : _bg_P) {
      // PolarSegment pseg = *sit;
      // for (pit = pseg.begin(); pit < pseg.end(); ++pit) {
      if (!pseg->IsCharged()) continue;
      for (auto &psite : *pseg) {
        charged_count += 1;
        system_dpl += psite->getPos() * psite->getQ00();
        if (psite->getRank() > 0) {
          system_dpl += psite->getQ1();
        }
      }
    }
    XTP_LOG(Log::info, *_log)
        << "    - System Q1: " << -system_dpl << "  (apply to " << charged_count
        << " polar sites)" << std::flush;
    vec atomic_compensation_dpl_system = -system_dpl / charged_count;
    [[maybe_unused]] int compensation_count =
        0;  // otherwise trigger unused warning as only in assert
    for (auto &pseg : _bg_P) {
      if (!pseg->IsCharged()) continue;
      // Dipole compensation type
      vec atomic_compensation_dpl = vec(0, 0, 0);
      if (_dipole_compensation_type == "system")
        atomic_compensation_dpl = atomic_compensation_dpl_system;
      else if (_dipole_compensation_type == "segment") {
        vec pseg_dpl = pseg->CalcTotD();
        atomic_compensation_dpl = -pseg_dpl / pseg->size();
      } else
        assert(false);  // Compensation type not implemented
      // Dipole compensation direction
      if (_dipole_compensation_direction == "xyz")
        ;
      else if (_dipole_compensation_direction == "z") {
        atomic_compensation_dpl = vec(0., 0., atomic_compensation_dpl(2));
      } else
        assert(false);  // Compensation direction not implemented
                        // Apply dipolar compensation
      for (auto &psite : *pseg) {
        compensation_count += 1;
        if (psite->getRank() < 1) {
          psite->setQ1(atomic_compensation_dpl);
          psite->setRank(1);
        } else {
          psite->setQ1(psite->getQ1() + atomic_compensation_dpl);
        }
      }
    }
    assert(compensation_count == charged_count);
  }

  // CHARGE APPROPRIATELY & DEPOLARIZE
  if (_do_restart) {
    // Restarting from previous iteration, hence no depolarization
    ;
  } else {
    for (auto &pseg : _bg_P) {
      for (auto &psite : *pseg) {
        psite->Depolarize();
      }
    }
  }

  // CALCULATE NET DIPOLE OF BGP & FGC
  vec netdpl_bgP = vec(0, 0, 0);
  double qzz_bgP = 0.0;
  for (auto &pseg : _bg_P) {
    for (auto &psite : *pseg) {

      netdpl_bgP += psite->getPos() * psite->getQ00();
      if (psite->getRank() > 0) {
        netdpl_bgP += psite->getQ1();
      }
      qzz_bgP += psite->getQ00() * psite->getPos()(2) * psite->getPos()(2);
    }
  }

  XTP_LOG(Log::info, *_log)
      << (format("Net dipole moment of background density")).str() << std::flush
      << (format("  o D(BGP) [e*nm]           = %1$+1.3f %2$+1.3f %3$+1.3f  ") %
          netdpl_bgP(0) % netdpl_bgP(1) % netdpl_bgP(2))
             .str();
  XTP_LOG(Log::info, *_log)
      << std::flush
      << (format("  o Sigma q|z|**2 [e*nm**2] = %1$+1.7f   ") % qzz_bgP)
      << std::flush;

  return;
}

PolarBackground::~PolarBackground() {
  std::vector<EWD::KVector *>::iterator kit;
  for (kit = _kvecs_2_0.begin(); kit < _kvecs_2_0.end(); ++kit) delete *kit;
  for (kit = _kvecs_1_0.begin(); kit < _kvecs_1_0.end(); ++kit) delete *kit;
  for (kit = _kvecs_0_0.begin(); kit < _kvecs_0_0.end(); ++kit) delete *kit;
  _kvecs_2_0.clear();
  _kvecs_1_0.clear();
  _kvecs_0_0.clear();
}

void PolarBackground::Checkpoint(int iter, bool converged) {
  XTP_LOG(Log::debug, *_log)
      << "  o Checkpointing (iteration " << iter << ") ... ";
  _ptop->setPolarizationIter(iter, converged);
  _ptop->SaveToDrive("bgp_check.ptop");
  XTP_LOG(Log::debug, *_log) << "done." << std::flush;
  return;
}

void PolarBackground::Polarize(int n_threads = 1) {

  XTP_LOG(Log::debug, *_log) << std::flush;
  XTP_LOG(Log::debug, *_log) << "System & Ewald parameters" << std::flush;
  XTP_LOG(Log::debug, *_log) << "  o Real-space unit cell:      " << _a << " x "
                             << _b << " x " << _c << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o Real-space c/o (guess):    " << _R_co << " nm" << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o na(max), nb(max), nc(max): " << _na_max << ", " << _nb_max << ", "
      << _nc_max << std::flush;
  XTP_LOG(Log::debug, *_log) << "  o 1st Brillouin zone:        " << _A << " x "
                             << _B << " x " << _C << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o Reciprocal-space c/o:      " << _K_co << " 1/nm" << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o R-K switching param.       " << _alpha << " 1/nm" << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o Unit-cell volume:          " << _LxLyLz << " nm**3" << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o LxLy (for 3D2D EW):        " << _LxLy << " nm**2" << std::flush;
  XTP_LOG(Log::debug, *_log)
      << "  o kx(max), ky(max), kz(max): " << _NA_max << ", " << _NB_max << ", "
      << _NC_max << std::flush;

  // TLogLevel dbg = Log::debug;

  //[-Wunused-variable]
  // TLogLevel inf = logINFO;
  // TLogLevel err = logERROR;

  Logger &log = *_log;
  _n_threads = n_threads;

  std::vector<PolarSeg *>::iterator sit1;
  std::vector<APolarSite *>::iterator pit1;
  std::vector<APolarSite *>::iterator pit2;

  /*
  Verify neutrality & depolarize
  Generate permanent fields (FP)
    o Converge intermolecular real-space contribution, remember cut-off
    o Converge reciprocal-space contribution, remember K-vectors
    o Calculate shape fields
    o Apply MOLECULAR foreground correction
  Induce to 1st order
  Loop until 2nd-order fields converged
    | Reset 2nd-order fields
    | (Re-)generate induction fields (FU)
    | o Real-space INTRAmolecular contribution to 2nd-order fields
    | o Real-space INTERmolecular contribution to 2nd-order fields
    | o Reciprocal-space contribution, work off remembered K-vectors
    | o Calculate shape fields
    | o Apply ATOMIC foreground correction
    | Induce to 2nd order
    + Check convergence
  Extract (or serialize) induction state to hard-drive
  */

  double rms = 0.0;
  int rms_count = 0;

  if (!_do_restart) {
    // I GENERATE PERMANENT FIELDS (FP)
    XTP_LOG(Log::debug, log) << std::flush;
    XTP_LOG(Log::debug, log) << "Generate permanent fields (FP)" << std::flush;
    // I.A Intermolecular real-space contribution
    XTP_LOG(Log::debug, log) << "  o Real-space, intermolecular" << std::flush;
    this->FX_RealSpace("FP_MODE", true);
    if (!_do_use_cutoff) {
      // I.B Reciprocal-space contribution
      XTP_LOG(Log::debug, log) << "  o Reciprocal-space" << std::flush;
      this->FX_ReciprocalSpace("SP_MODE", "FP_MODE", true);
      // I.C Shape fields
      XTP_LOG(Log::debug, log)
          << "  o Shape fields ('" << _shape << "')" << std::flush;
      _ewdactor.FP12_ShapeField_At_By(_bg_P, _bg_P, _shape, _LxLyLz);
      // I.D Molecular ERF self-interaction correction
      XTP_LOG(Log::debug, log) << "  o Molecular SI correction" << std::flush;
      for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
          for (pit2 = (*sit1)->begin(); pit2 < (*sit1)->end(); ++pit2) {
            rms += _ewdactor.FP12_ERF_At_By(*(*pit1), *(*pit2));
            rms_count += 1;
          }
        }
      }
      rms = sqrt(rms / rms_count) * int2V_m;
    }

    // TEASER OUTPUT PERMANENT FIELDS
    XTP_LOG(Log::debug, *_log)
        << std::flush << "Foreground fields:" << std::flush;
    int fieldCount = 0;
    for (sit1 = _bg_P.begin() + 16; sit1 < _bg_P.end(); ++sit1) {
      PolarSeg *pseg = *sit1;
      Segment seg = _top->getSegment(pseg->getId());
      XTP_LOG(Log::debug, *_log) << "ID = " << pseg->getId() << " ("
                                 << seg.getId() << ") " << std::flush;
      for (pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
        vec fp = (*pit1)->getFieldP();
        XTP_LOG(Log::debug, *_log)
            << (format("FP = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m") %
                (fp(0) * int2V_m) % (fp(1) * int2V_m) % (fp(2) * int2V_m))
                   .str()
            << std::flush;
        fieldCount += 1;
        if (fieldCount > 10) {
          XTP_LOG(Log::debug, *_log) << "FP = ... ... ..." << std::flush;
          break;
        }
      }
      if (fieldCount > 10) break;
    }

    // II INDUCE TO 1ST ORDER
    XTP_LOG(Log::debug, log)
        << std::flush << "Induce to first order" << std::flush;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        (*pit1)->InduceDirect();
      }
    }

    // Checkpointing for potential restart
    if (_do_checkpointing) this->Checkpoint(0, false);
  } else {
    XTP_LOG(Log::debug, log) << std::flush;
    XTP_LOG(Log::debug, log)
        << "Restarting from checkpoint => "
        << "Permanent fields already generated." << std::flush;
  }

  // III CONVERGE INDUCTION FIELDS (FU)
  // NOTE Real-space neighbours should always be regenerated during first
  //      induction iteration
  // NOTE Reciprocal-space vectors need only be regenerated if not already
  //      done on permanent level
  int iter = (_do_restart) ? _restart_from_iter + 1 : 1;
  int setup_nbs_iter = (_do_restart) ? _restart_from_iter + 1 : 1;
  int generate_kvecs_iter = (_do_restart) ? _restart_from_iter + 1 : -1;

  XTP_LOG(Log::debug, log) << "  o Setup real-space neighbours at iteration "
                           << setup_nbs_iter << std::flush;
  XTP_LOG(Log::debug, log) << "  o Generate k-vectors at iteration "
                           << generate_kvecs_iter << std::flush;

  int max_iter = iter + _max_iter;
  double epstol = 1e-3;
  for (; iter != max_iter; ++iter) {
    XTP_LOG(Log::debug, log) << std::flush;
    XTP_LOG(Log::debug, log) << "Iter " << iter << " started" << std::flush;

    // III.A Reset 2nd-order fields (FU)
    XTP_LOG(Log::debug, log) << "  o Reset fields (FU)" << std::flush;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        (*pit1)->ResetFieldU();
      }
    }
    // III.B (Re-)generate induction fields FU
    XTP_LOG(Log::debug, log)
        << "  o (Re-)generate induction fields" << std::flush;
    // (1) Real-space intramolecular contribution
    XTP_LOG(Log::debug, log) << "  o Real-space, intramolecular" << std::flush;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = pit1 + 1; pit2 < (*sit1)->end(); ++pit2) {
          _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2));
          _ewdactor.FU12_ERFC_At_By(*(*pit2), *(*pit1));
          //_actor.BiasIndu(*(*pit1),*(*pit2));
          //_actor.FieldIndu(*(*pit1),*(*pit2));
        }
      }
    }
    // (2) Real-space intermolecular contribution
    bool do_setup_nbs = (iter == setup_nbs_iter) ? true : false;
    bool generate_kvecs = (iter == generate_kvecs_iter) ? true : false;
    this->FX_RealSpace("FU_MODE", do_setup_nbs);
    if (!_do_use_cutoff) {
      // (3) Reciprocal-space contribution
      XTP_LOG(Log::debug, log) << "  o Reciprocal-space" << std::flush;
      this->FX_ReciprocalSpace("SU_MODE", "FU_MODE", generate_kvecs);
      // (4) Calculate shape fields
      XTP_LOG(Log::debug, log)
          << "  o Shape fields ('" << _shape << "')" << std::flush;
      _ewdactor.FU12_ShapeField_At_By(_bg_P, _bg_P, _shape, _LxLyLz);
      // (5) Apply atomic ERF self-interaction correction
      XTP_LOG(Log::debug, log) << "  o Atomic SI correction" << std::flush;
      rms = 0.0;
      rms_count = 0;
      for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
          rms += _ewdactor.FU12_ERF_At_By(*(*pit1), *(*pit1));
          rms_count += 1;
        }
      }
      rms = sqrt(rms / rms_count) * int2V_m;
    }

    // TEASER OUTPUT INDUCTION FIELDS
    XTP_LOG(Log::debug, *_log)
        << std::flush << "Foreground fields:" << std::flush;
    int fieldCount = 0;
    for (sit1 = _bg_P.begin() + 288; sit1 < _bg_P.end(); ++sit1) {
      PolarSeg *pseg = *sit1;
      Segment seg = _top->getSegment(pseg->getId());
      XTP_LOG(Log::debug, *_log) << "ID = " << pseg->getId() << " ("
                                 << seg.getId() << ") " << std::flush;
      for (pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
        vec fu = (*pit1)->getFieldU();
        XTP_LOG(Log::debug, *_log)
            << (format("FU = (%1$+1.7e %2$+1.7e %3$+1.7e) V/m") %
                (fu(0) * int2V_m) % (fu(1) * int2V_m) % (fu(2) * int2V_m))
                   .str()
            << std::flush;
        fieldCount += 1;
        if (fieldCount > 10) {
          XTP_LOG(Log::debug, *_log)
              << "FU = ... ... ..." << std::flush << std::flush;
          break;
        }
      }
      if (fieldCount > 10) break;
    }

    // III.C Induce again
    XTP_LOG(Log::debug, log) << "  o Induce again" << std::flush;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        (*pit1)->Induce(_polar_wSOR_N);
      }
    }

    // III.D Check for convergence
    XTP_LOG(Log::debug, log) << "  o Convergence check" << std::flush;
    bool converged = true;
    double maxdU = -1;
    double avgdU = 0.0;
    int baseN = 0;
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        double dU = (*pit1)->HistdU();
        avgdU += dU;
        ++baseN;
        if (dU > maxdU) {
          maxdU = dU;
        }
        if (dU > epstol) {
          converged = false;
        }
      }
    }
    avgdU /= baseN;
    if (avgdU < epstol * 0.1) {
      converged = true;
    }
    if (_do_checkpointing) this->Checkpoint(iter, converged);
    if (converged) {
      _converged = true;
      XTP_LOG(Log::debug, log) << std::flush;
      XTP_LOG(Log::debug, log) << ":: Converged induction fields" << std::flush;
      break;
    } else if (iter == max_iter) {
      throw std::runtime_error("Not converged.");
      break;
    }
  }

  if (iter == max_iter) {
    XTP_LOG(Log::debug, log) << std::flush;
    XTP_LOG(Log::debug, log)
        << "Reached maximum number of iterations. Stop here." << std::flush;
  }

  if (Log::verbose()) {
    std::ofstream ofs;
    ofs.open("ewdbgpol.indu_state.tab", std::ofstream::out);
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      PolarSeg *pseg = *sit1;

      //[-Wunused-variable]
      // Segment *seg = _top->getSegment(pseg->getId());

      for (pit1 = pseg->begin(); pit1 < pseg->end(); ++pit1) {
        vec fp = (*pit1)->getFieldP();
        vec fu = (*pit1)->getFieldU();
        vec u1 = (*pit1)->getU1();
        vec pos = (*pit1)->getPos();
        ofs << (format("SEGID2 %1$4d   ") % (pseg->getId()));
        ofs << (format("XYZ6 %1$+1.7e %2$+1.7e %3$+1.7e    ") % (pos(0)) %
                (pos(1)) % (pos(2)))
                   .str();
        ofs << (format("FP10 %1$+1.7e %2$+1.7e %3$+1.7e    ") %
                (fp(0) * int2V_m) % (fp(1) * int2V_m) % (fp(2) * int2V_m))
                   .str();
        ofs << (format("FU14 %1$+1.7e %2$+1.7e %3$+1.7e    ") %
                (fu(0) * int2V_m) % (fu(1) * int2V_m) % (fu(2) * int2V_m))
                   .str();
        ofs << (format("U118 %1$+1.7e %2$+1.7e %3$+1.7e   ") % (u1(0)) %
                (u1(1)) % (u1(2)))
                   .str()
            << std::endl;
      }
    }
    ofs.close();
  }

  return;
}

// ========================================================================== //
// FP/FU REAL SPACE
// ========================================================================== //

void PolarBackground::RThread::FP_FieldCalc() {

  std::vector<PolarSeg *>::iterator sit1;
  std::vector<APolarSite *>::iterator pit1;
  std::vector<PolarSeg *>::iterator sit2;
  std::vector<APolarSite *>::iterator pit2;
  std::vector<PolarNb *>::iterator nit;
  _not_converged_count = 0;

  // CLEAR POLAR NEIGHBOR-LIST BEFORE SET-UP
  if (_do_setup_nbs) {
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "   - Clearing polar nb-list" << std::endl;
    }
    for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
      (*sit1)->ClearPolarNbs();
    }
  }

  double R_co_sum = 0.0;
  int R_co_sum_count = 0;

  for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
    PolarSeg *pseg1 = *sit1;
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "\rMST DBG     - Progress " << pseg1->getId() << "/"
          << _full_bg_P.size() << std::flush;
    }

    // GENERATE NEIGHBOUR SHELLS
    double dR_shell = 0.5;
    double R_co_max = 2 * _master->_R_co;
    int N_shells = int(R_co_max / dR_shell) + 1;
    std::vector<std::vector<PolarNb *> > shelled_nbs;
    shelled_nbs.resize(N_shells);
    [[maybe_unused]] unsigned int allocated_count = 0;
    [[maybe_unused]] unsigned int deleted_count = 0;

    for (sit2 = _full_bg_P.begin(); sit2 < _full_bg_P.end(); ++sit2) {
      PolarSeg *pseg2 = *sit2;
      // Active segment?
      if (!pseg2->IsCharged() && !pseg2->IsPolarizable()) continue;
      for (int na = -_master->_na_max; na < _master->_na_max + 1; ++na) {
        for (int nb = -_master->_nb_max; nb < _master->_nb_max + 1; ++nb) {
          for (int nc = -_master->_nc_max; nc < _master->_nc_max + 1; ++nc) {
            // Identical?
            if (na == 0 && nb == 0 && nc == 0 && pseg1 == pseg2) continue;
            if (na == 0 && nb == 0 && nc == 0 &&
                pseg1->getId() == pseg2->getId())
              assert(false);
            // Apply periodic-boundary correction, check c/o, shift
            vec dr12_pbc = _master->_top->PbShortestConnect(pseg1->getPos(),
                                                            pseg2->getPos());
            vec dr12_dir = pseg2->getPos() - pseg1->getPos();
            // Image box correction
            vec L = na * _master->_a + nb * _master->_b + nc * _master->_c;
            vec dr12_pbc_L = dr12_pbc + L;
            vec s22x_L = dr12_pbc_L - dr12_dir;
            double R = (dr12_pbc_L).norm();
            if (R > R_co_max) continue;
            // Add to shell
            int shell_idx = int(R / dR_shell);
            shelled_nbs[shell_idx].push_back(
                new PolarNb(pseg2, dr12_pbc_L, s22x_L));
            allocated_count += 1;
          }
        }
      }
    }

    int shell_idx = 0;
    double shell_R = 0.;
    // LONG-RANGE TREATMENT: REAL-SPACE SUM
    if (!_master->_do_use_cutoff) {
      // SUM OVER CONSECUTIVE SHELLS & STORE NBS FOR REUSE
      bool converged = false;
      int charged_nbs_count = 0;
      for (int sidx = 0; sidx < N_shells; ++sidx) {
        // Figure out shell parameters
        shell_idx = sidx;
        shell_R = (sidx + 1) * dR_shell;
        std::vector<PolarNb *> &nb_shell = shelled_nbs[sidx];
        if (nb_shell.size() < 1) continue;
        double shell_rms = 0.0;
        int shell_rms_count = 0;
        // Interact ...
        for (nit = nb_shell.begin(); nit < nb_shell.end(); ++nit) {
          PolarSeg *pseg2 = (*nit)->getNb();
          // Add neighbour for later use
          pseg1->AddPolarNb(*nit);
          if (!pseg2->IsCharged()) continue;
          charged_nbs_count += 1;
          if (((*nit)->getR()).norm() > R_co_max) assert(false);
          // Interact taking into account shift
          for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
            for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
              shell_rms +=
                  _ewdactor.FP12_ERFC_At_By(*(*pit1), *(*pit2), (*nit)->getS());
              shell_rms_count += 1;
            }
          }
        }
        // Determine convergence - measure is the energy of a dipole
        // of size 0.1*e*nm summed over the shell in an rms manner
        shell_rms = sqrt(shell_rms / shell_rms_count) * int2V_m;
        double e_measure = shell_rms * 1e-10 * shell_rms_count;
        if (shell_rms_count > 0 && e_measure <= _master->_crit_dE &&
            shell_R >= _master->_R_co) {
          converged = true;
          break;
        }
      }
      if (!converged && charged_nbs_count > 0) {
        _not_converged_count += 1;
      }
      if (charged_nbs_count > 0) {
        R_co_sum += shell_R;
        R_co_sum_count += 1;
      }
    }
    // CUTOFF TREATMENT: STANDARD REAL-SPACE SUM
    else {
      double epsilon = 1.;
      for (int sidx = 0; sidx < N_shells; ++sidx) {
        // Still in cutoff sphere?
        if ((sidx + 1) * dR_shell > _master->_R_co) break;
        // Figure out shell parameters
        shell_idx = sidx;
        shell_R = (sidx + 1) * dR_shell;
        std::vector<PolarNb *> &nb_shell = shelled_nbs[sidx];
        if (nb_shell.size() < 1) continue;
        for (nit = nb_shell.begin(); nit < nb_shell.end(); ++nit) {
          PolarSeg *pseg2 = (*nit)->getNb();
          // Add neighbour for later use
          pseg1->AddPolarNb(*nit);
          if (!pseg2->IsCharged()) continue;
          if (((*nit)->getR()).norm() > R_co_max) assert(false);
          // Interact taking into account shift
          for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
            for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
              _actor.BiasStat(*(*pit1), *(*pit2), (*nit)->getS());
              _actor.FieldPerm_At_By(*(*pit1), *(*pit2), epsilon);
            }
          }
        }
      }
    }

    // DELETE ALL NEIGHBOURS THAT WERE NOT NEEDED TO CONVERGE SUM
    for (int sidx = shell_idx + 1; sidx < N_shells; ++sidx) {
      std::vector<PolarNb *> &nb_shell = shelled_nbs[sidx];
      for (nit = nb_shell.begin(); nit < nb_shell.end(); ++nit) {
        delete *nit;
        deleted_count += 1;
      }
    }
    shelled_nbs.clear();
    assert(pseg1->PolarNbs().size() + deleted_count == allocated_count);
  }
  if (R_co_sum_count == 0)
    _avg_R_co = 0.0;
  else
    _avg_R_co = R_co_sum / R_co_sum_count;
  return;
}

void PolarBackground::RThread::FU_FieldCalc() {
  std::vector<PolarSeg *>::iterator sit1;
  std::vector<APolarSite *>::iterator pit1;
  std::vector<PolarSeg *>::iterator sit2;
  std::vector<APolarSite *>::iterator pit2;
  std::vector<PolarNb *>::iterator nit;

  // SET-UP NB CONTAINER
  if (_do_setup_nbs) {

    // CLEAR POLAR NEIGHBOR-LIST BEFORE SET-UP
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "   - Clearing polar nb-list" << std::endl;
    }
    for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
      (*sit1)->ClearPolarNbs();
    }

    double R_co_sum = 0.0;
    int R_co_sum_count = 0;

    for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
      PolarSeg *pseg1 = *sit1;
      if (Log::verbose()) {
        XTP_LOG(Log::debug, *(_master->_log))
            << "\rMST DBG     - Progress " << pseg1->getId() << "/"
            << _full_bg_P.size() << std::flush;
      }

      // GENERATE NEIGHBOUR SHELLS
      double dR_shell = 0.5;
      double R_co_max = 2 * _master->_R_co;
      int N_shells = int(R_co_max / dR_shell) + 1;
      std::vector<std::vector<PolarNb *> > shelled_nbs;
      shelled_nbs.resize(N_shells);
      [[maybe_unused]] unsigned int allocated_count = 0;
      [[maybe_unused]] unsigned int deleted_count = 0;

      for (sit2 = _full_bg_P.begin(); sit2 < _full_bg_P.end(); ++sit2) {
        PolarSeg *pseg2 = *sit2;
        // Active segment?
        if (!pseg2->IsCharged() && !pseg2->IsPolarizable()) continue;
        for (int na = -_master->_na_max; na < _master->_na_max + 1; ++na) {
          for (int nb = -_master->_nb_max; nb < _master->_nb_max + 1; ++nb) {
            for (int nc = -_master->_nc_max; nc < _master->_nc_max + 1; ++nc) {
              // Identical?
              if (na == 0 && nb == 0 && nc == 0 && pseg1 == pseg2) continue;
              if (na == 0 && nb == 0 && nc == 0 &&
                  pseg1->getId() == pseg2->getId())
                assert(false);
              // Apply periodic-boundary correction, check c/o, shift
              vec dr12_pbc = _master->_top->PbShortestConnect(pseg1->getPos(),
                                                              pseg2->getPos());
              vec dr12_dir = pseg2->getPos() - pseg1->getPos();
              // Image box correction
              vec L = na * _master->_a + nb * _master->_b + nc * _master->_c;
              vec dr12_pbc_L = dr12_pbc + L;
              vec s22x_L = dr12_pbc_L - dr12_dir;
              double R = (dr12_pbc_L).norm();
              if (R > R_co_max) continue;
              // Add to shell
              int shell_idx = int(R / dR_shell);
              shelled_nbs[shell_idx].push_back(
                  new PolarNb(pseg2, dr12_pbc_L, s22x_L));
              allocated_count += 1;
            }
          }
        }
      }

      int shell_idx = 0;
      double shell_R = 0.;
      // LONG-RANGE TREATMENT: REAL-SPACE SUM
      if (!_master->_do_use_cutoff) {
        // SUM OVER CONSECUTIVE SHELLS & STORE NBS FOR REUSE
        bool converged = false;
        int polarizable_nbs_count = 0;
        for (int sidx = 0; sidx < N_shells; ++sidx) {
          // Figure out shell parameters
          shell_idx = sidx;
          shell_R = (sidx + 1) * dR_shell;
          std::vector<PolarNb *> &nb_shell = shelled_nbs[sidx];
          if (nb_shell.size() < 1) continue;
          double shell_rms = 0.0;
          int shell_rms_count = 0;
          // Interact ...
          for (nit = nb_shell.begin(); nit < nb_shell.end(); ++nit) {
            PolarSeg *pseg2 = (*nit)->getNb();
            // Add neighbour for later use
            pseg1->AddPolarNb(*nit);
            if (!pseg2->IsPolarizable()) continue;
            polarizable_nbs_count += 1;
            if (((*nit)->getR()).norm() > R_co_max) assert(false);
            // Interact taking into account shift
            for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
              for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
                shell_rms += _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2),
                                                       (*nit)->getS());
                shell_rms_count += 1;
              }
            }
          }
          // Determine convergence - measure is the energy of a dipole
          // of size 0.1*e*nm summed over the shell in an rms manner
          shell_rms = sqrt(shell_rms / shell_rms_count) * int2V_m;
          double e_measure = shell_rms * 1e-10 * shell_rms_count;
          if (shell_rms_count > 0 && e_measure <= _master->_crit_dE &&
              shell_R >= _master->_R_co) {
            converged = true;
            break;
          }
        }
        if (!converged && polarizable_nbs_count > 0) {
          _not_converged_count += 1;
        }
        if (polarizable_nbs_count > 0) {
          R_co_sum += shell_R;
          R_co_sum_count += 1;
        }
      }
      // CUTOFF TREATMENT: STANDARD REAL-SPACE SUM
      else {
        double epsilon = 1.;
        for (int sidx = 0; sidx < N_shells; ++sidx) {
          // Still in cutoff sphere?
          if ((sidx + 1) * dR_shell > _master->_R_co) break;
          // Figure out shell parameters
          shell_idx = sidx;
          shell_R = (sidx + 1) * dR_shell;
          std::vector<PolarNb *> &nb_shell = shelled_nbs[sidx];
          if (nb_shell.size() < 1) continue;
          for (nit = nb_shell.begin(); nit < nb_shell.end(); ++nit) {
            PolarSeg *pseg2 = (*nit)->getNb();
            // Add neighbour for later use
            pseg1->AddPolarNb(*nit);
            if (!pseg2->IsPolarizable()) continue;
            if (((*nit)->getR()).norm() > R_co_max) assert(false);
            // Interact taking into account shift
            for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
              for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
                _actor.BiasIndu(*(*pit1), *(*pit2), (*nit)->getS());
                _actor.FieldIndu_At_By(*(*pit1), *(*pit2), epsilon);
              }
            }
          }
        }
      }

      // DELETE ALL NEIGHBOURS THAT WERE NOT NEEDED TO CONVERGE SUM
      for (int sidx = shell_idx + 1; sidx < N_shells; ++sidx) {
        std::vector<PolarNb *> &nb_shell = shelled_nbs[sidx];
        for (nit = nb_shell.begin(); nit < nb_shell.end(); ++nit) {
          delete *nit;
          deleted_count += 1;
        }
      }
      shelled_nbs.clear();
      assert(pseg1->PolarNbs().size() + deleted_count == allocated_count);
    }
    if (R_co_sum_count == 0)
      _avg_R_co = 0.0;
    else
      _avg_R_co = R_co_sum / R_co_sum_count;
  }
  // REUSE NB CONTAINER
  else {
    double rms = 0.0;
    int rms_count = 0;
    for (sit1 = _part_bg_P.begin(); sit1 < _part_bg_P.end(); ++sit1) {
      PolarSeg *pseg1 = *sit1;
      if (Log::verbose()) {
        XTP_LOG(Log::debug, *(_master->_log))
            << "\rMST DBG     - Progress " << pseg1->getId() << "/"
            << _full_bg_P.size() << std::flush;
      }
      // LONG-RANGE TREATMENT: REAL-SPACE SUM
      if (!_master->_do_use_cutoff) {
        for (nit = pseg1->PolarNbs().begin(); nit < pseg1->PolarNbs().end();
             ++nit) {
          PolarSeg *pseg2 = (*nit)->getNb();
          if (!pseg2->IsPolarizable()) continue;
          // Interact taking into account shift
          for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
            for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
              rms +=
                  _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2), (*nit)->getS());
              rms_count += 1;
            }
          }
        }
      }
      // CUTOFF TREATMENT: STANDARD REAL-SPACE SUM
      else {
        double epsilon = 1.;
        for (nit = pseg1->PolarNbs().begin(); nit < pseg1->PolarNbs().end();
             ++nit) {
          PolarSeg *pseg2 = (*nit)->getNb();
          if (!pseg2->IsPolarizable()) continue;
          // Interact taking into account shift
          for (pit1 = pseg1->begin(); pit1 < pseg1->end(); ++pit1) {
            for (pit2 = pseg2->begin(); pit2 < pseg2->end(); ++pit2) {
              _actor.BiasIndu(*(*pit1), *(*pit2), (*nit)->getS());
              _actor.FieldIndu_At_By(*(*pit1), *(*pit2), epsilon);
            }
          }
        }
      }
    }
    if (rms_count > 0) rms = sqrt(rms / rms_count) * int2V_m;
    _avg_R_co = 0.0;  // not available, since reusing NB container
  }
  return;
}

void PolarBackground::FX_RealSpace(std::string mode, bool do_setup_nbs) {

  RThread prototype(this, do_setup_nbs);
  ThreadForce<RThread, PrototypeCreator> tforce;
  ThreadForce<RThread, PrototypeCreator>::iterator tfit;
  tforce.setPrototype(&prototype);
  tforce.Initialize(_n_threads);

  tforce.AddSharedInput<std::vector<PolarSeg *> >(_bg_P);
  tforce.AddAtomicInput<PolarSeg *>(_bg_P);
  tforce.AssignMode<std::string>(mode);

  // Output workload
  XTP_LOG(Log::debug, *_log) << "    - Thread workload = [ ";
  for (tfit = tforce.begin(); tfit != tforce.end(); ++tfit) {
    XTP_LOG(Log::debug, *_log)
        << (format("%1$1.2f%% ") % (*tfit)->Workload(mode));
  }
  XTP_LOG(Log::debug, *_log) << "]" << std::flush;

  // Start & wait
  XTP_LOG(Log::debug, *_log)
      << "    - Start & wait until done" << std::flush << std::flush;
  _log->setPreface(Log::debug, "");
  tforce.StartAndWait();
  _log->setPreface(Log::debug, "\nMST DBG");

  // Assert convergence
  int not_converged_count = 0;
  for (tfit = tforce.begin(); tfit != tforce.end(); ++tfit)
    not_converged_count += (*tfit)->NotConverged();
  if (not_converged_count == 0)
    XTP_LOG(Log::debug, *_log) << "    - Converged" << std::flush;
  else
    XTP_LOG(Log::error, *_log) << "    - ERROR " << not_converged_count
                               << " items not converged." << std::flush;

  // Neighbor-list info: radius & neighbours/site
  double avg_R_co = 0;
  for (tfit = tforce.begin(); tfit != tforce.end(); ++tfit)
    avg_R_co += 0.01 * (*tfit)->Workload(mode) * (*tfit)->AvgRco();
  XTP_LOG(Log::debug, *_log)
      << "    - Real-space nb-list set: <R(c/o)> = " << avg_R_co << std::flush;
  int total_nbs_count = 0;
  std::vector<PolarSeg *>::iterator sit1;
  for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1)
    total_nbs_count += (*sit1)->PolarNbs().size();
  XTP_LOG(Log::debug, *_log)
      << "    - Real-space nb-list set: <nbs/seg> = "
      << (double)total_nbs_count / double(_bg_P.size()) << std::flush;

  return;
}

// ========================================================================== //
// FP & FU RECIPROCAL SPACE
// ========================================================================== //

void PolarBackground::KThread::SP_SFactorCalc() {

  // Calculate structure factors for each k and store with KVector
  int kvec_count = 0;
  for (std::vector<EWD::KVector *>::iterator kit = _part_kvecs.begin();
       kit < _part_kvecs.end(); ++kit) {
    kvec_count += 1;
    EWD::cmplx sfactor =
        _ewdactor.PStructureAmplitude(_full_bg_P, (*kit)->getK());
    (*kit)->setStructureFactor(sfactor);
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "\rMST DBG     - " << _current_mode << "(SP) Progress "
          << kvec_count << "/" << _part_kvecs.size() << std::flush;
    }
  }

  return;
}

void PolarBackground::KThread::FP_KFieldCalc() {

  _rms_sum_re = 0.0;
  _sum_im = 0.0;

  double rV = 1. / _master->_LxLyLz;

  // Increment fields within _part_bg_P for each k-vector
  int kvec_count = 0;
  for (std::vector<EWD::KVector *>::iterator kit = _full_kvecs.begin();
       kit < _full_kvecs.end(); ++kit) {
    kvec_count += 1;
    vec k = (*kit)->getK();
    EWD::cmplx S = (*kit)->getStructureFactor();
    EWD::cmplx f_rms = _ewdactor.FP12_At_ByS2(k, _part_bg_P, S, rV);
    _rms_sum_re += f_rms._re;
    _sum_im += f_rms._im;
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "\rMST DBG     - " << _current_mode << "(FP) Progress "
          << kvec_count << "/" << _full_kvecs.size() << std::flush;
    }
  }

  return;
}

void PolarBackground::KThread::SU_SFactorCalc() {

  // Calculate structure factors for each k and store with KVector
  int kvec_count = 0;
  for (std::vector<EWD::KVector *>::iterator kit = _part_kvecs.begin();
       kit < _part_kvecs.end(); ++kit) {
    kvec_count += 1;
    EWD::cmplx sfactor =
        _ewdactor.UStructureAmplitude(_full_bg_P, (*kit)->getK());
    (*kit)->setStructureFactor(sfactor);
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "\rMST DBG     - " << _current_mode << "(SU) Progress "
          << kvec_count << "/" << _part_kvecs.size() << std::flush;
    }
  }

  return;
}

void PolarBackground::KThread::FU_KFieldCalc() {

  _rms_sum_re = 0.0;
  _sum_im = 0.0;

  double rV = 1. / _master->_LxLyLz;

  // Increment fields within _part_bg_P for each k-vector
  int kvec_count = 0;
  for (std::vector<EWD::KVector *>::iterator kit = _full_kvecs.begin();
       kit < _full_kvecs.end(); ++kit) {
    kvec_count += 1;
    vec k = (*kit)->getK();
    EWD::cmplx S = (*kit)->getStructureFactor();
    EWD::cmplx f_rms = _ewdactor.FU12_At_ByS2(k, _part_bg_P, S, rV);
    _rms_sum_re += f_rms._re;
    _sum_im += f_rms._im;
    if (Log::verbose()) {
      XTP_LOG(Log::debug, *(_master->_log))
          << "\rMST DBG     - " << _current_mode << "(FU) Progress "
          << kvec_count << "/" << _full_kvecs.size() << std::flush;
    }
  }

  return;
}

void PolarBackground::FX_ReciprocalSpace(std::string mode1, std::string mode2,
                                         bool generate_kvecs) {

  double sum_re = 0.0;
  double sum_im = 0.0;
  _field_converged_K = false;

  if (mode1 == "SP_MODE")
    assert(mode2 == "FP_MODE");
  else if (mode1 == "SU_MODE")
    assert(mode2 == "FU_MODE");
  else
    assert(false);

  // GENERATE K-VECTORS
  if (generate_kvecs) {
    XTP_LOG(Log::debug, *_log) << "  o Generate k-vectors" << std::flush;
    _log->setPreface(Log::debug, "\nMST DBG     - ");
    GenerateKVectors(_bg_P, _bg_P);
    _log->setPreface(Log::debug, "\nMST DBG");
  }

  // THREAD FORCE & PROTOTYPE
  KThread prototype(this);
  ThreadForce<KThread, PrototypeCreator> threadforce;
  ThreadForce<KThread, PrototypeCreator>::iterator tfit;
  threadforce.setPrototype(&prototype);
  threadforce.Initialize(_n_threads);
  threadforce.AddSharedInput<std::vector<PolarSeg *> >(_bg_P);
  threadforce.AddAtomicInput<PolarSeg *>(_bg_P);

  // TWO COMPONENTS ZERO, ONE NON-ZERO
  XTP_LOG(Log::debug, *_log)
      << "  o Two components zero, one non-zero" << std::flush;
  // Assign k-vectors
  threadforce.AddSharedInput<std::vector<KVector *> >(_kvecs_2_0);
  threadforce.AddAtomicInput<EWD::KVector *>(_kvecs_2_0);
  // Compute structure factors
  XTP_LOG(Log::debug, *_log) << std::flush;
  threadforce.AssignMode<std::string>(mode1);
  _log->setPreface(Log::debug, "");
  threadforce.StartAndWait();
  _log->setPreface(Log::debug, "\nMST DBG");
  // Increment fields
  XTP_LOG(Log::debug, *_log) << std::flush;
  threadforce.AssignMode<std::string>(mode2);
  _log->setPreface(Log::debug, "");
  threadforce.StartAndWait();
  _log->setPreface(Log::debug, "\nMST DBG");
  // Collect r.m.s. information
  double rms_sum_re = 0.0;
  for (tfit = threadforce.begin(); tfit != threadforce.end(); ++tfit) {
    rms_sum_re += (*tfit)->Workload(mode2) * (*tfit)->_rms_sum_re;
    sum_im += (*tfit)->_sum_im;
  }
  double shell_rms = sqrt(rms_sum_re / double(_kvecs_2_0.size())) * int2V_m;
  double e_measure = shell_rms * 1e-10 * double(_kvecs_2_0.size());
  if (_kvecs_2_0.size() > 0) {
    XTP_LOG(Log::debug, *_log)
        << (format("    - M = %1$04d   G = %2$+1.3e   dF(rms) = %3$+1.3e V/m   "
                   "[1eA => %4$+1.3e eV]") %
            _kvecs_2_0.size() % 0.0 % shell_rms % e_measure)
               .str()
        << std::flush;
  }
  // Clear k-vector containers
  threadforce.Reset<std::string>(mode1);
  threadforce.Reset<std::string>(mode2);

  // ONE COMPONENT ZERO, TWO NON-ZERO
  XTP_LOG(Log::debug, *_log)
      << "  o One component zero, two non-zero" << std::flush;
  double crit_grade = 1. * _kxyz_s1s2_norm;
  bool converged10 = false;
  std::vector<EWD::KVector *>::iterator kit;
  kit = _kvecs_1_0.begin();
  while (!converged10 && kit < _kvecs_1_0.end()) {
    // Construct k-vector shell from critical grade
    std::vector<KVector *> shell_kvecs;
    while (kit < _kvecs_1_0.end()) {
      if ((*kit)->getGrade() < crit_grade) break;
      shell_kvecs.push_back(*kit);
      ++kit;
    }
    if (shell_kvecs.size() > 0) {
      // Assign k-vectors
      threadforce.AddSharedInput<std::vector<KVector *> >(shell_kvecs);
      threadforce.AddAtomicInput<EWD::KVector *>(shell_kvecs);
      // Compute structure factors
      XTP_LOG(Log::debug, *_log) << std::flush;
      threadforce.AssignMode<std::string>(mode1);
      _log->setPreface(Log::debug, "");
      threadforce.StartAndWait();
      _log->setPreface(Log::debug, "\nMST DBG");
      // Increment fields
      XTP_LOG(Log::debug, *_log) << std::flush;
      threadforce.AssignMode<std::string>(mode2);
      _log->setPreface(Log::debug, "");
      threadforce.StartAndWait();
      _log->setPreface(Log::debug, "\nMST DBG");
      // Collect r.m.s. information
      rms_sum_re = 0.0;
      for (tfit = threadforce.begin(); tfit != threadforce.end(); ++tfit) {
        rms_sum_re += (*tfit)->Workload(mode2) * (*tfit)->_rms_sum_re;
        sum_im += (*tfit)->_sum_im;
      }
      shell_rms = sqrt(rms_sum_re / double(shell_kvecs.size())) * int2V_m;
      e_measure = shell_rms * 1e-10 * double(shell_kvecs.size());
      // Log & assert convergence
      XTP_LOG(Log::debug, *_log)
          << (format("    - M = %1$04d   G = %2$+1.3e   dF(rms) = %3$+1.3e V/m "
                     "  [1eA => %4$+1.3e eV]") %
              shell_kvecs.size() % crit_grade % shell_rms % e_measure)
                 .str()
          << std::flush;
      if (shell_kvecs.size() > 10 && e_measure <= _crit_dE) {
        XTP_LOG(Log::debug, *_log)
            << (format("    :: RE %1$+1.7e IM %2$+1.7e") %
                (sqrt(sum_re) * int2V_m) % (sum_im * int2V_m))
                   .str()
            << std::flush;
        converged10 = true;
      }
      // Clear k-vector containers
      threadforce.Reset<std::string>(mode1);
      threadforce.Reset<std::string>(mode2);
      shell_kvecs.clear();
    }
    crit_grade *= 0.1;
  }

  // ZERO COMPONENTS ZERO, THREE NON-ZERO
  XTP_LOG(Log::debug, *_log)
      << "  o Zero components zero, three non-zero" << std::flush;
  crit_grade = 1. * _kxyz_s1s2_norm;
  bool converged00 = false;
  kit = _kvecs_0_0.begin();
  while (!converged00 && kit < _kvecs_0_0.end()) {
    // Construct k-vector shell from critical grade
    std::vector<KVector *> shell_kvecs;
    while (kit < _kvecs_0_0.end()) {
      if ((*kit)->getGrade() < crit_grade) break;
      shell_kvecs.push_back(*kit);
      ++kit;
    }
    if (shell_kvecs.size() > 0) {
      // Assign k-vectors
      threadforce.AddSharedInput<std::vector<KVector *> >(shell_kvecs);
      threadforce.AddAtomicInput<EWD::KVector *>(shell_kvecs);
      // Compute structure factors
      XTP_LOG(Log::debug, *_log) << std::flush;
      threadforce.AssignMode<std::string>(mode1);
      _log->setPreface(Log::debug, "");
      threadforce.StartAndWait();
      _log->setPreface(Log::debug, "\nMST DBG");
      // Increment fields
      XTP_LOG(Log::debug, *_log) << std::flush;
      threadforce.AssignMode<std::string>(mode2);
      _log->setPreface(Log::debug, "");
      threadforce.StartAndWait();
      _log->setPreface(Log::debug, "\nMST DBG");
      // Collect r.m.s. information
      rms_sum_re = 0.0;
      for (tfit = threadforce.begin(); tfit != threadforce.end(); ++tfit) {
        rms_sum_re += (*tfit)->Workload(mode2) * (*tfit)->_rms_sum_re;
        sum_im += (*tfit)->_sum_im;
      }
      shell_rms = sqrt(rms_sum_re / double(shell_kvecs.size())) * int2V_m;
      e_measure = shell_rms * 1e-10 * double(shell_kvecs.size());
      // Log & assert convergence
      XTP_LOG(Log::debug, *_log)
          << (format("    - M = %1$04d   G = %2$+1.3e   dF(rms) = %3$+1.3e V/m "
                     "  [1eA => %4$+1.3e eV]") %
              shell_kvecs.size() % crit_grade % shell_rms % e_measure)
                 .str()
          << std::flush;
      if (shell_kvecs.size() > 10 && e_measure <= _crit_dE) {
        XTP_LOG(Log::debug, *_log)
            << (format("    :: RE %1$+1.7e IM %2$+1.7e") %
                (sqrt(sum_re) * int2V_m) % (sum_im * int2V_m))
                   .str()
            << std::flush;
        converged00 = true;
      }
      // Clear k-vector containers
      threadforce.Reset<std::string>(mode1);
      threadforce.Reset<std::string>(mode2);
      shell_kvecs.clear();
    }
    crit_grade *= 0.1;
  }

  _field_converged_K = converged10 && converged00;

  if (_field_converged_K) {
    XTP_LOG(Log::debug, *_log)
        << (format("  o Converged to precision, {2-1}, {1-2}, {0-3}."))
        << std::flush;
  }

  return;
}

void PolarBackground::GenerateKVectors(std::vector<PolarSeg *> &ps1,
                                       std::vector<PolarSeg *> &ps2) {

  // Take care of norm for grading function
  // All three components non-zero
  //              S(kx)*S(ky)*S(kz)
  // G = A(k) * ---------------------
  //            (<S(kx)><S(ky)><S(kz)>)**(2/3)
  // Component i zero
  //                   S(kj)*S(kk)
  // G = A(k) * -------------------------
  //             (<S(kj)><S(kk)>)**(1/2)
  // Components i,j zero
  // => All S(k) calculated anyway, no need to grade
  // We can use the same grading function if we set
  //
  // S(ki=0) = <S(ki)>**(2/3) (<S(kj)><S(kk)>)**(1/6)

  std::vector<EWD::KVector *> kvecs_2_0;  // 2 components zero
  std::vector<EWD::KVector *> kvecs_1_0;  // 1 component zero
  std::vector<EWD::KVector *> kvecs_0_0;  // 0 components zero

  // CONTAINERS FOR GRADING K-VECTORS
  std::vector<double> kx_s1s2;
  kx_s1s2.push_back(1);
  std::vector<double> ky_s1s2;
  ky_s1s2.push_back(1);
  std::vector<double> kz_s1s2;
  kz_s1s2.push_back(1);
  double avg_kx_s1s2 = 0.0;
  double avg_ky_s1s2 = 0.0;
  double avg_kz_s1s2 = 0.0;

  // TWO COMPONENTS ZERO, ONE NON-ZERO
  XTP_LOG(Log::debug, *_log)
      << "Generating K-vectors: Exploring K resonances" << std::flush;
  for (int i = 1; i < _NA_max + 1; ++i) {
    vec k = +i * _A;
    EWD::triple<EWD::cmplx> ppuu_posk = _ewdactor.S1S2(k, ps1, ps2);
    kx_s1s2.push_back(0.5 * std::abs(ppuu_posk._pp._re));
    avg_kx_s1s2 += 0.5 * std::abs(ppuu_posk._pp._re);
    EWD::KVector *kvec_pos = new EWD::KVector(+1 * k, 0.);
    EWD::KVector *kvec_neg = new EWD::KVector(-1 * k, 0.);
    kvecs_2_0.push_back(kvec_pos);
    kvecs_2_0.push_back(kvec_neg);
  }
  avg_kx_s1s2 /= _NA_max;

  for (int i = 1; i < _NB_max + 1; ++i) {
    vec k = +i * _B;
    EWD::triple<EWD::cmplx> ppuu_posk = _ewdactor.S1S2(k, ps1, ps2);
    ky_s1s2.push_back(0.5 * std::abs(ppuu_posk._pp._re));
    avg_ky_s1s2 += 0.5 * std::abs(ppuu_posk._pp._re);
    EWD::KVector *kvec_pos = new EWD::KVector(+1 * k, 0);
    EWD::KVector *kvec_neg = new EWD::KVector(-1 * k, 0);
    kvecs_2_0.push_back(kvec_pos);
    kvecs_2_0.push_back(kvec_neg);
  }
  avg_ky_s1s2 /= _NB_max;

  for (int i = 1; i < _NC_max + 1; ++i) {
    vec k = +i * _C;
    EWD::triple<EWD::cmplx> ppuu_posk = _ewdactor.S1S2(k, ps1, ps2);
    kz_s1s2.push_back(0.5 * std::abs(ppuu_posk._pp._re));
    avg_kz_s1s2 += 0.5 * std::abs(ppuu_posk._pp._re);
    EWD::KVector *kvec_pos = new EWD::KVector(+1 * k, 0);
    EWD::KVector *kvec_neg = new EWD::KVector(-1 * k, 0);
    kvecs_2_0.push_back(kvec_pos);
    kvecs_2_0.push_back(kvec_neg);
  }
  avg_kz_s1s2 /= _NC_max;
  const double int2eV =
      1 / (4 * Pi * 8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
  double kxyz_s1s2_norm =
      1. / pow(avg_kx_s1s2 * avg_ky_s1s2 * avg_kz_s1s2, 2. / 3.) * int2eV /
      _LxLyLz;
  kx_s1s2[0] =
      pow(avg_ky_s1s2 * avg_kz_s1s2, 1. / 6.) * pow(avg_kx_s1s2, 2. / 3.);
  ky_s1s2[0] =
      pow(avg_kz_s1s2 * avg_kx_s1s2, 1. / 6.) * pow(avg_ky_s1s2, 2. / 3.);
  kz_s1s2[0] =
      pow(avg_kx_s1s2 * avg_ky_s1s2, 1. / 6.) * pow(avg_kz_s1s2, 2. / 3.);

  // ONE COMPONENT ZERO, TWO NON-ZERO
  XTP_LOG(Log::debug, *_log)
      << "K-planes through origin: Applying K resonances" << std::flush;

  int kx, ky, kz;
  kx = 0;
  for (ky = -_NB_max; ky < _NB_max + 1; ++ky) {
    if (ky == 0) continue;
    for (kz = -_NC_max; kz < _NC_max + 1; ++kz) {
      if (kz == 0) continue;
      vec k = kx * _A + ky * _B + kz * _C;
      double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] *
                     ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] *
                     kxyz_s1s2_norm;
      EWD::KVector *kvec = new EWD::KVector(k, grade);
      kvecs_1_0.push_back(kvec);
    }
  }
  ky = 0;
  for (kx = -_NA_max; kx < _NA_max + 1; ++kx) {
    if (kx == 0) continue;
    for (kz = -_NC_max; kz < _NC_max + 1; ++kz) {
      if (kz == 0) continue;
      vec k = kx * _A + ky * _B + kz * _C;
      double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] *
                     ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] *
                     kxyz_s1s2_norm;
      EWD::KVector *kvec = new EWD::KVector(k, grade);
      kvecs_1_0.push_back(kvec);
    }
  }
  kz = 0;
  for (kx = -_NA_max; kx < _NA_max + 1; ++kx) {
    if (kx == 0) continue;
    for (ky = -_NB_max; ky < _NB_max + 1; ++ky) {
      if (ky == 0) continue;
      vec k = kx * _A + ky * _B + kz * _C;
      double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] *
                     ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] *
                     kxyz_s1s2_norm;
      EWD::KVector *kvec = new EWD::KVector(k, grade);
      kvecs_1_0.push_back(kvec);
    }
  }
  _kvecsort._p = 1e-300;
  std::sort(kvecs_1_0.begin(), kvecs_1_0.end(), _kvecsort);

  // ZERO COMPONENTS ZERO, THREE NON-ZERO
  XTP_LOG(Log::debug, *_log)
      << "K-space (off-axis): Applying K resonances" << std::flush;

  for (kx = -_NA_max; kx < _NA_max + 1; ++kx) {
    if (kx == 0) continue;
    for (ky = -_NB_max; ky < _NB_max + 1; ++ky) {
      if (ky == 0) continue;
      for (kz = -_NC_max; kz < _NC_max + 1; ++kz) {
        if (kz == 0) continue;
        vec k = kx * _A + ky * _B + kz * _C;
        double grade = _ewdactor.Ark2Expk2(k) * kx_s1s2[std::abs(kx)] *
                       ky_s1s2[std::abs(ky)] * kz_s1s2[std::abs(kz)] *
                       kxyz_s1s2_norm;
        EWD::KVector *kvec = new EWD::KVector(k, grade);
        kvecs_0_0.push_back(kvec);
      }
    }
  }

  _kvecsort._p = 1e-300;
  std::sort(kvecs_0_0.begin(), kvecs_0_0.end(), _kvecsort);

  // STORE K-VECTORS
  std::vector<EWD::KVector *>::iterator kvit;
  for (kvit = _kvecs_2_0.begin(); kvit < _kvecs_2_0.end(); ++kvit) delete *kvit;
  for (kvit = _kvecs_1_0.begin(); kvit < _kvecs_1_0.end(); ++kvit) delete *kvit;
  for (kvit = _kvecs_0_0.begin(); kvit < _kvecs_0_0.end(); ++kvit) delete *kvit;
  _kvecs_2_0.clear();
  _kvecs_1_0.clear();
  _kvecs_0_0.clear();

  _kvecs_2_0 = kvecs_2_0;
  _kvecs_1_0 = kvecs_1_0;
  _kvecs_0_0 = kvecs_0_0;
  _kxyz_s1s2_norm = kxyz_s1s2_norm;

  return;
}

}  // namespace EWD
}  // namespace xtp
}  // namespace votca