#include "votca/xtp/ewald/polarbackground.h"
#include "votca/tools/property.h"
#include <boost/format.hpp>
#include <mutex>
#include <set>
#include <votca/tools/globals.h>

namespace votca {
namespace xtp {
namespace EWD {

using boost::format;

namespace {
// Debug/experimental (this session). Shared by both the per-pair
// static- and induced-field dumps below -- see their own gate's
// context (search this file for IsPerPairDumpTarget) for what this is
// for. Matches the same target sample used on the new (non-legacy)
// code's own equivalent extension, so the two can be compared
// directly, site by site.
bool IsPerPairDumpTarget(votca::Index segment_id) {
  static const std::set<votca::Index> sample_targets = {
      0, 100, 200, 300, 400, 500, 600, 700, 800, 900};
  return sample_targets.count(segment_id) > 0;
}
}  // namespace

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
  // Debug/experimental, see _debug_dump_first_iteration_and_stop's own
  // declaration.
  if (opt.exists(pfx + ".polarmethod.debug_dump_first_iteration_and_stop"))
    _debug_dump_first_iteration_and_stop =
        opt.get(pfx + ".polarmethod.debug_dump_first_iteration_and_stop")
            .as<bool>();
  else
    _debug_dump_first_iteration_and_stop = false;
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

    // Debug/experimental, no counterpart other than a matching flag
    // added to the "new" (non-legacy) Ewald code this session -- see
    // _debug_dump_first_iteration_and_stop's own declaration for what
    // this exists for. Placed HERE, right after InduceDirect() and
    // before the main iteration loop, not inside that loop's own first
    // pass (where an earlier version of this dump lived) -- confirmed
    // directly this session that the main loop's own first "Induce
    // again" call is NOT the genuine no-coupling-yet baseline it looks
    // like: InduceDirect() already set every site's own U1 to a real,
    // nonzero value beforehand, so by the time the main loop's own
    // first "(Re-)generate induction fields (FU)" step runs, FU is
    // already nonzero (a genuine dipole-dipole field from those already-
    // set U1 values, via FU12_ERFC_At_By), not the clean, coupling-free
    // ground state the "new" code's own JOR genuinely starts from
    // (x=0). InduceDirect() itself, by contrast, IS that clean
    // baseline: mu = -P*F_perm, no wSOR relaxation, no FU term at all
    // (FUx=FUy=FUz are still their own zero-initialized default here).
    if (_debug_dump_first_iteration_and_stop) {
      std::ofstream dump("legacy_iteration1_dump.csv");
      dump.precision(15);
      dump << "segment_id,site_index,pos_x_bohr,pos_y_bohr,pos_z_bohr,"
              "mu_x_bohr,mu_y_bohr,mu_z_bohr,mu_x_enm,mu_y_enm,mu_z_enm,"
              "Fperm_x_native,Fperm_y_native,Fperm_z_native,"
              "Pxx_native,Pyy_native,Pzz_native,"
              "Pxy_native,Pxz_native,Pyz_native\n";
      // votca::tools::conv (bohr2nm/nm2bohr) is not accessible from this
      // legacy code path (a different, older unit-constants convention
      // is used throughout this file already -- see e.g. int2V_m just
      // above); legacy's own internal positions/dipoles are already in
      // nm directly (see this class's own header comments and its own
      // "e*nm" convention, confirmed directly against APolarSite::
      // HistdU's own source this session), so the nm columns need no
      // conversion at all here, and the bohr columns are the ones
      // computed via nm2bohr -- the opposite direction from the new
      // code's own dump, which starts in bohr and converts to nm.
      //
      // Fperm/P are dumped in legacy's own NATIVE units, unconverted --
      // deliberately not run through nm2bohr (a linear length
      // conversion), since field and polarizability are DIFFERENT
      // powers of length than dipole moment (E ~ 1/length^2,
      // polarizability ~ length^3), and applying a linear dipole-moment
      // conversion factor to either would be wrong, not just imprecise.
      // The precise name of legacy's own native field/polarizability
      // unit convention was not independently confirmed this session
      // (unlike the "e*nm" dipole convention, confirmed directly
      // against APolarSite::HistdU's own source) -- "native" is an
      // honest label given that, not a claim about what unit it
      // actually is.
      const double nm2bohr = 18.897259886;
      for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        int site_index = 0;
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
          vec pos = (*pit1)->getPos();
          vec mu_enm = (*pit1)->getU1();
          vec f_perm = (*pit1)->getFieldP();
          vec p_diag = (*pit1)->getPDiag();
          vec p_offdiag = (*pit1)->getPOffDiag();
          dump << (*sit1)->getId() << "," << site_index << ","
               << pos.x() * nm2bohr << "," << pos.y() * nm2bohr << ","
               << pos.z() * nm2bohr << "," << mu_enm.x() * nm2bohr << ","
               << mu_enm.y() * nm2bohr << "," << mu_enm.z() * nm2bohr << ","
               << mu_enm.x() << "," << mu_enm.y() << "," << mu_enm.z()
               << "," << f_perm.x() << "," << f_perm.y() << ","
               << f_perm.z() << "," << p_diag.x() << "," << p_diag.y()
               << "," << p_diag.z() << "," << p_offdiag.x() << ","
               << p_offdiag.y() << "," << p_offdiag.z() << "\n";
          ++site_index;
        }
      }
      dump.close();
      // Deliberately does NOT throw/stop here anymore -- see the
      // second dump point (right after this main loop's own first
      // "III.C Induce again" call) for why this run continues one
      // more step before stopping.
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
    // Debug/experimental, staged FU capture -- see
    // _debug_dump_first_iteration_and_stop's own declaration. FU is
    // captured cumulatively after each of the five stages below
    // (intramolecular, real-space intermolecular, reciprocal-space,
    // shape, self-interaction), in site-visitation order, so that a
    // later cross-codebase comparison can isolate exactly which stage
    // a real discrepancy enters at, rather than only seeing it already
    // mixed into the final summed FU -- deliberately not attempting to
    // separate the two real-space contributions (intramolecular here,
    // intermolecular in FX_RealSpace) any further than this, since
    // FX_RealSpace's own internal structure is not something this
    // session inspected directly.
    std::vector<vec> fu_stage_a, fu_stage_b, fu_stage_c, fu_stage_d,
        fu_stage_e;
    // (1) Real-space intramolecular contribution
    XTP_LOG(Log::debug, log) << "  o Real-space, intramolecular" << std::flush;
    std::ofstream intrapair_dump;
    if (_debug_dump_first_iteration_and_stop) {
      intrapair_dump.open("legacy_intrapair_dump.csv");
      intrapair_dump.precision(15);
      intrapair_dump << "segment_id,site_i,site_j,r_native,tu3,ta1_tu3,l3,l5,"
                        "B0,B1,B2\n";
    }
    for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
      for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
        for (pit2 = pit1 + 1; pit2 < (*sit1)->end(); ++pit2) {
          _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2));
          // Debug/experimental -- see EwdInteractor::GetDebugPairState's
          // own declaration. Captured here, right after this specific
          // call, since FU12_ERFC_At_By's own ApplyBiasPolar +
          // UpdateAllBls calls are what set _ewdactor's own internal
          // r/tu3/l3/l5/B0/B1/B2 state for THIS (pit1, pit2) pair -- the
          // very next line's own call (pit2, pit1) overwrites that same
          // state with the pair's own reverse-direction values before
          // this dump would otherwise see it.
          if (_debug_dump_first_iteration_and_stop) {
            auto ds = _ewdactor.GetDebugPairState();
            intrapair_dump << (*sit1)->getId() << ","
                          << std::distance((*sit1)->begin(), pit1) << ","
                          << std::distance((*sit1)->begin(), pit2) << ","
                          << ds.r << "," << ds.tu3 << "," << ds.ta1_tu3 << ","
                          << ds.l3 << "," << ds.l5 << "," << ds.B0 << ","
                          << ds.B1 << "," << ds.B2 << "\n";
          }
          _ewdactor.FU12_ERFC_At_By(*(*pit2), *(*pit1));
          //_actor.BiasIndu(*(*pit1),*(*pit2));
          //_actor.FieldIndu(*(*pit1),*(*pit2));
        }
      }
    }
    if (_debug_dump_first_iteration_and_stop) {
      intrapair_dump.close();
      for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1)
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1)
          fu_stage_a.push_back((*pit1)->getFieldU());
    }
    // (2) Real-space intermolecular contribution
    bool do_setup_nbs = (iter == setup_nbs_iter) ? true : false;
    bool generate_kvecs = (iter == generate_kvecs_iter) ? true : false;
    this->FX_RealSpace("FU_MODE", do_setup_nbs);
    if (_debug_dump_first_iteration_and_stop) {
      for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1)
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1)
          fu_stage_b.push_back((*pit1)->getFieldU());
    }
    if (!_do_use_cutoff) {
      // (3) Reciprocal-space contribution
      XTP_LOG(Log::debug, log) << "  o Reciprocal-space" << std::flush;
      this->FX_ReciprocalSpace("SU_MODE", "FU_MODE", generate_kvecs);
      if (_debug_dump_first_iteration_and_stop) {
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1)
          for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1)
            fu_stage_c.push_back((*pit1)->getFieldU());
      }
      // (4) Calculate shape fields
      XTP_LOG(Log::debug, log)
          << "  o Shape fields ('" << _shape << "')" << std::flush;
      _ewdactor.FU12_ShapeField_At_By(_bg_P, _bg_P, _shape, _LxLyLz);
      if (_debug_dump_first_iteration_and_stop) {
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1)
          for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1)
            fu_stage_d.push_back((*pit1)->getFieldU());
      }
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
      if (_debug_dump_first_iteration_and_stop) {
        for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1)
          for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1)
            fu_stage_e.push_back((*pit1)->getFieldU());
      }
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

    // Debug/experimental, second of two dump points -- see the first
    // (right after InduceDirect(), before this loop even starts) for
    // the base rationale. This one exists specifically to test where
    // induced-induced coupling first enters: mu_1 (this site's own
    // pre-Induce value, from InduceDirect() and already dumped
    // separately to legacy_iteration1_dump.csv) is confirmed to match
    // the "new" code's own mu_1 (once compared correctly, up to the
    // already-understood wSOR-vs-unrelaxed scaling); mu_2 here (this
    // call's own result) and FU (the field that produced it, via
    // getFieldU(), computed from every OTHER site's own mu_1 through
    // FU12_ERFC_At_By and the intermolecular/reciprocal/shape
    // contributions above) are the genuinely new information. A first,
    // aggregate comparison this session found a real (scattered, not a
    // clean scalar) mismatch here -- fu_stage_a..e (captured above,
    // cumulatively, right after each of the five physical contributions
    // to FU) exist so a follow-up comparison can isolate exactly which
    // stage that mismatch enters at, rather than only seeing it already
    // mixed into the final summed FU.
    if (_debug_dump_first_iteration_and_stop) {
      if (_do_use_cutoff) {
        throw std::runtime_error(
            "PolarBackground: debug_dump_first_iteration_and_stop's own "
            "staged FU capture (fu_stage_c/d/e) assumes the non-cutoff "
            "Ewald path (reciprocal-space/shape/self-interaction all "
            "run) -- this run has _do_use_cutoff=true, where those "
            "stages never execute, so those columns would be empty.");
      }
      std::ofstream dump("legacy_iteration2_dump.csv");
      dump.precision(15);
      dump << "segment_id,site_index,"
              "mu2_x_enm,mu2_y_enm,mu2_z_enm,"
              "FU_x_native,FU_y_native,FU_z_native,"
              "FUa_x_native,FUa_y_native,FUa_z_native,"
              "FUb_x_native,FUb_y_native,FUb_z_native,"
              "FUc_x_native,FUc_y_native,FUc_z_native,"
              "FUd_x_native,FUd_y_native,FUd_z_native,"
              "FUe_x_native,FUe_y_native,FUe_z_native\n";
      std::size_t n = 0;
      for (sit1 = _bg_P.begin(); sit1 < _bg_P.end(); ++sit1) {
        int site_index = 0;
        for (pit1 = (*sit1)->begin(); pit1 < (*sit1)->end(); ++pit1) {
          vec mu2_enm = (*pit1)->getU1();
          vec fu = (*pit1)->getFieldU();
          dump << (*sit1)->getId() << "," << site_index << ","
               << mu2_enm.x() << "," << mu2_enm.y() << "," << mu2_enm.z()
               << "," << fu.x() << "," << fu.y() << "," << fu.z() << ","
               << fu_stage_a[n].x() << "," << fu_stage_a[n].y() << ","
               << fu_stage_a[n].z() << "," << fu_stage_b[n].x() << ","
               << fu_stage_b[n].y() << "," << fu_stage_b[n].z() << ","
               << fu_stage_c[n].x() << "," << fu_stage_c[n].y() << ","
               << fu_stage_c[n].z() << "," << fu_stage_d[n].x() << ","
               << fu_stage_d[n].y() << "," << fu_stage_d[n].z() << ","
               << fu_stage_e[n].x() << "," << fu_stage_e[n].y() << ","
               << fu_stage_e[n].z() << "\n";
          ++site_index;
          ++n;
        }
      }
      dump.close();
      throw std::runtime_error(
          "PolarBackground: stopped after dumping iteration-2 induced "
          "dipoles/coupling field to legacy_iteration2_dump.csv "
          "(debug_dump_first_iteration_and_stop was set) -- this is not "
          "a real failure, just the requested early stop.");
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
              // Debug/experimental -- per-pair STATIC field dump,
              // mirroring FU_FieldCalc's own equivalent
              // legacy_perpairfield_dump.csv addition (see that one's
              // own comment for why), but for FP12_ERFC_At_By (this
              // pair's own contribution to the PERMANENT field, FPx/
              // FPy/FPz) instead of FU12_ERFC_At_By's own induced-field
              // one. Added specifically to let a direct per-pair
              // comparison isolate whether ApplyStaticField's own
              // new-code counterpart (already checked against
              // ApplyInducedField's own equivalent, which came back
              // clean) is where a real, still-unexplained discrepancy
              // in the AGGREGATE F_perm/FUb comparison actually lives.
              vec fp_before = (*pit1)->getFieldP();
              shell_rms +=
                  _ewdactor.FP12_ERFC_At_By(*(*pit1), *(*pit2), (*nit)->getS());
              shell_rms_count += 1;
              if (_master->_debug_dump_first_iteration_and_stop &&
                  IsPerPairDumpTarget(pseg1->getId())) {
                vec fp_after = (*pit1)->getFieldP();
                vec pair_field = fp_after - fp_before;
                static std::once_flag pps_dump_header_once;
                static std::mutex pps_dump_mutex;
                static std::ofstream pps_dump;
                std::call_once(pps_dump_header_once, []() {
                  pps_dump.open("legacy_perpairstaticfield_dump.csv",
                              std::ios::trunc);
                  pps_dump.precision(15);
                  pps_dump << "target_segment_id,target_site_index,"
                              "source_segment_id,source_site_index,"
                              "field_x_native,field_y_native,"
                              "field_z_native\n";
                });
                std::lock_guard<std::mutex> lock(pps_dump_mutex);
                pps_dump << pseg1->getId() << ","
                       << (pit1 - pseg1->begin()) << ","
                       << pseg2->getId() << ","
                       << (pit2 - pseg2->begin()) << ","
                       << pair_field.x() << "," << pair_field.y()
                       << "," << pair_field.z() << "\n";
              }
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
              // Debug/experimental -- see this session's own added
              // EwaldRealSpaceSum::DumpNeighborList (new-code side) for
              // what this is for. Dumps every pseg1 (not just one target
              // segment -- a first version of this restricted to a
              // single segment risked hiding a genuine PBC-scheme
              // discrepancy that only shows up for segments in a
              // different geometric relationship to the box, e.g. near
              // a boundary rather than deep in the interior), so
              // target_segment_id is written as its own column. Every
              // pseg1 hits its own "first write for this pseg1" moment
              // (allocated_count==1), so a std::call_once guards the
              // one-time file truncation + header (never more than
              // once, regardless of which pseg1/thread reaches it
              // first) and a mutex guards every individual row write
              // against concurrent multi-thread corruption (this loop
              // runs inside a genuinely multi-threaded worker -- see
              // FX_RealSpace's own ThreadForce setup).
              // CRITICAL: the stream itself is opened exactly once (a
              // static, function-local std::ofstream, left open for the
              // rest of the process's own lifetime) and reused for
              // every row -- an earlier version of this opened a fresh
              // std::ofstream on EVERY row, under the same mutex,
              // serializing every worker thread onto a full file
              // open+seek-to-end+close for each individual neighbor
              // pair; for a real system this made the whole run many
              // orders of magnitude slower (observed: unfinished after
              // 24 hours) than the underlying physics warrants.
              if (_master->_debug_dump_first_iteration_and_stop) {
                static std::once_flag nb_dump_header_once;
                static std::mutex nb_dump_mutex;
                static std::ofstream nb_dump;
                std::call_once(nb_dump_header_once, []() {
                  nb_dump.open("legacy_neighborlist_dump.csv",
                              std::ios::trunc);
                  nb_dump.precision(15);
                  nb_dump << "target_segment_id,source_segment_id,"
                             "r_nm,t_x_nm,t_y_nm,t_z_nm\n";
                });
                std::lock_guard<std::mutex> lock(nb_dump_mutex);
                nb_dump << pseg1->getId() << "," << pseg2->getId() << ","
                       << R << "," << dr12_pbc_L.x() << ","
                       << dr12_pbc_L.y() << "," << dr12_pbc_L.z() << "\n";
                static std::size_t nb_dump_row_count = 0;
                if (++nb_dump_row_count % 10000 == 0) {
                  nb_dump.flush();
                }
              }
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
                // Debug/experimental -- per-pair field dump, mirroring
                // this session's own EwaldRealSpaceSum::
                // DumpPerPairFieldAppend (new-code side). Isolates each
                // individual pair's own field contribution by
                // differencing pit1's own FUx/FUy/FUz before/after this
                // one FU12_ERFC_At_By call -- the SAME call already
                // confirmed (via the !_do_use_cutoff branch check) to
                // be the one this run's own config actually executes,
                // after an earlier turn this session mistakenly looked
                // at the OTHER (_do_use_cutoff-only, dead-for-this-run)
                // branch's own _actor/XInteractor calls instead.
                // Restricted to a single target segment (id 0) to keep
                // file size manageable at this per-pair granularity --
                // unlike the earlier per-segment neighbor-list dump,
                // which dumped every target since segment-level
                // granularity kept it small.
                vec fu_before = (*pit1)->getFieldU();
                shell_rms += _ewdactor.FU12_ERFC_At_By(*(*pit1), *(*pit2),
                                                       (*nit)->getS());
                shell_rms_count += 1;
                if (_master->_debug_dump_first_iteration_and_stop &&
                    IsPerPairDumpTarget(pseg1->getId())) {
                  vec fu_after = (*pit1)->getFieldU();
                  vec pair_field = fu_after - fu_before;
                  static std::once_flag pp_dump_header_once;
                  static std::mutex pp_dump_mutex;
                  static std::ofstream pp_dump;
                  std::call_once(pp_dump_header_once, []() {
                    pp_dump.open("legacy_perpairfield_dump.csv",
                                std::ios::trunc);
                    pp_dump.precision(15);
                    pp_dump << "target_segment_id,target_site_index,"
                               "source_segment_id,source_site_index,"
                               "field_x_native,field_y_native,"
                               "field_z_native\n";
                  });
                  std::lock_guard<std::mutex> lock(pp_dump_mutex);
                  pp_dump << pseg1->getId() << ","
                         << (pit1 - pseg1->begin()) << ","
                         << pseg2->getId() << ","
                         << (pit2 - pseg2->begin()) << ","
                         << pair_field.x() << "," << pair_field.y()
                         << "," << pair_field.z() << "\n";
                }
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

  // Debug/experimental -- see this session's own log-line addition
  // right after _field_converged_K's own assignment for what this is
  // for. Captured here, before kit gets reassigned for the 0-0 group
  // right below, since it's the SAME iterator/cursor reused for both
  // groups.
  const Index kvecs_1_0_used = Index(std::distance(_kvecs_1_0.begin(), kit));

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

  // Debug/experimental (this session): reports the ACTUAL number of
  // k-vectors this run consumed before declaring convergence -- not
  // K_co (an outer safety bound this adaptive, grade-sorted scheme
  // rarely if ever fully exhausts) and not the total candidate set
  // GenerateKVectors built (which also just reflects K_co, not actual
  // usage). Added to let a direct, apples-to-apples comparison against
  // the new (non-legacy) code's own NumKVectors() -- which DOES always
  // fully evaluate its own k_max sphere, no adaptive stopping -- be
  // possible, after a raw K_co-vs-k_max comparison was correctly
  // pointed out as not meaningful on its own, precisely because of
  // this adaptive stopping behavior.
  const Index kvecs_0_0_used = Index(std::distance(_kvecs_0_0.begin(), kit));
  const Index kvecs_total_used =
      Index(_kvecs_2_0.size()) + kvecs_1_0_used + kvecs_0_0_used;
  XTP_LOG(Log::debug, *_log)
      << (format("  o K-vectors actually used (this run): 2-0=%1$d (all, "
                 "no convergence check) + 1-0=%2$d (of %3$d candidates) + "
                 "0-0=%4$d (of %5$d candidates) = %6$d total") %
          _kvecs_2_0.size() % kvecs_1_0_used % _kvecs_1_0.size() %
          kvecs_0_0_used % _kvecs_0_0.size() % kvecs_total_used)
             .str()
      << std::flush;

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
