/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#include <cmath>
#include <vector>

// Local includes
#include "votca/xtp/kspace.h"

namespace votca {
namespace xtp {

void KSpace::computeStaticField() {
  for (Index i = 0; i < Index(ewaldSegments_.size()); ++i) {
    EwdSegment& seg = ewaldSegments_[i];
    for (EwdSite& site : seg) {
      for (const KVector& kvec : kvector_list_) {
        site.addToStaticField(
            fourPiVolume * kvec.getVector() * kvec.getAk() *
            (ii * std::exp(ii * kvec.getVector().dot(site.getPos())) *
             std::conj(kvec.getSk()))
                .real());
      }
    }
  }
}

/**
 * This functions computes the TOTAL kspace field.
 * Note that the kvectors are computed and setup at initialization of this class
 * If the segments/sites contain induced dipoles (which they should if loaded
 * from a state file) then the kvectors Sk values etc will all have been
 * computed with the
 */
void KSpace::computeTotalField(PolarSegment& seg) {
  EwdSegment& currentSeg = ewaldSegments_[seg.getId()];
  for (EwdSite& site : currentSeg) {
    for (const KVector& kvec : kvector_list_) {
      site.addToStaticField(
          fourPiVolume * kvec.getVector() * kvec.getAk() *
          (ii * std::exp(ii * kvec.getVector().dot(site.getPos())) *
           std::conj(kvec.getSk()))
              .real());
    }
  }
}

void KSpace::addInducedDipoleInteractionTo(Eigen::MatrixXd& result) {
  for (const KVector& kvec : kvector_list_) {
    Eigen::VectorXcd SkInteraction =
        getSkInteractionVector(kvec.getVector()).conjugate();
    for (Index i = 0; i < Index(ewaldSegments_.size()); ++i) {
      EwdSegment& seg = ewaldSegments_[i];
      Index startRow = segmentOffSet[i];
      for (EwdSite& site : seg) {
        double kPosDot = kvec.getVector().dot(site.getPos());
        Eigen::VectorXd realPart =
            (std::exp(ii * kPosDot) * SkInteraction).imag();
        result.row(startRow + 0) +=
            -fourPiVolume * kvec.getVector()[0] * kvec.getAk() * realPart;
        result.row(startRow + 1) +=
            -fourPiVolume * kvec.getVector()[1] * kvec.getAk() * realPart;
        result.row(startRow + 2) +=
            -fourPiVolume * kvec.getVector()[2] * kvec.getAk() * realPart;
        startRow += 3;
      }
    }
  }
}

Eigen::VectorXcd KSpace::getSkInteractionVector(
    const Eigen::Vector3d& kvector) {
  Eigen::VectorXcd result(systemSize);
  result.fill(0);

  for (Index segId = 0; segId < Index(ewaldSegments_.size()); ++segId) {
    EwdSegment& currentSeg = ewaldSegments_[segId];
    Index offset = segmentOffSet[segId];
    Index siteCounter = 0;
    for (const EwdSite& site : currentSeg) {
      std::complex<double> expKR = std::exp(ii * kvector.dot(site.getPos()));
      result.segment<3>(offset + siteCounter) += expKR * ii * kvector;
      siteCounter += 3;
    }
  }
  return result;
}

void KSpace::addShapeCorrectionTo(Eigen::MatrixXd& result) {
  for (Index i = 0; i < Index(ewaldSegments_.size()); ++i) {
    EwdSegment& seg = ewaldSegments_[i];
    Index startRow = segmentOffSet[i];
    for (EwdSite& site : seg) {
      switch (options_.shape) {
        case Shape::xyslab:
          for (Index j = 2; j < systemSize; j = j + 3) {
            result(startRow + 2, j) += fourPiVolume;
          }
          startRow += 3;
          break;
        case Shape::cube:
        case Shape::sphere:
          for (Index j = 0; j < systemSize; j = j + 3) {
            result(startRow + 0, j + 0) += (fourPiVolume / 3.0);
            result(startRow + 1, j + 1) += (fourPiVolume / 3.0);
            result(startRow + 2, j + 2) += (fourPiVolume / 3.0);
          }
          startRow += 3;
          break;
        default:
          throw std::runtime_error("Shape not implemented.");
      }
    }
  }
}

void KSpace::addSICorrectionTo(Eigen::MatrixXd& result) {
  // Intramolecular correction
  for (Index segId = 0; segId < Index(ewaldSegments_.size()); ++segId) {
    EwdSegment& currentSeg = ewaldSegments_[segId];
    Index startRow = segmentOffSet[segId];
    Index startCol = startRow;
    for (Index site_ind1 = 0; site_ind1 < currentSeg.size(); ++site_ind1) {
      result.block<3, 3>(startRow + 3 * site_ind1, startCol + 3 * site_ind1) -=
          inducedDipoleInteractionAtBy(currentSeg[site_ind1],
                                       currentSeg[site_ind1]);
    }
  }
}

Eigen::Vector3d KSpace::computeShapeField() {
  Eigen::Vector3d dip_tot = Eigen::Vector3d::Zero();
  Eigen::Vector3d shapeField = Eigen::Vector3d::Zero();

  // Compute total dipole
  for (const EwdSegment& seg : ewaldSegments_) {
    for (const EwdSite& sit : seg) {
      dip_tot += sit.getCharge() * sit.getPos();
      dip_tot += sit.getTotalDipole();
    }
  }

  switch (options_.shape) {
    case Shape::xyslab:
      shapeField[2] = -fourPiVolume * dip_tot[2];
      break;
    case Shape::cube:
    case Shape::sphere:
      shapeField = -(fourPiVolume / 3.0) * dip_tot;
      break;
    default:
      throw std::runtime_error("Shape not implemented.");
  }
  return shapeField;
}

void KSpace::applyAPeriodicCorrection(PolarSegment& seg,
                                      std::vector<SegId> pCloud_indices) {
  EwdSegment currentSeg = ewaldSegments_[seg.getId()];
  for (SegId& id : pCloud_indices) {
    if (id.Id() != seg.getId()) {
      EwdSegment& nbSeg = ewaldSegments_[id.Id()];
      for (EwdSite& site : currentSeg) {
        for (EwdSite& site2 : nbSeg) {
        }
      }
    }
  }
}

Eigen::Vector3d KSpace::kspaceCorrectionFieldAtBy(EwdSite& site,
                                                  const EwdSite& nbSite) {
  computeDistanceVariables(site.getPos() - nbSite.getPos());
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  // charge
  field += -nbSite.getCharge() * dr * rR3s;
  if (rank > 0) {  // dipole
    field += nbSite.getTotalDipole() * rR3s;
    field += -rR5s * dr * dr.dot(nbSite.getStaticDipole());
    if (rank > 1) {  // quadrupole
      // Using that the trace of a quadrupole contributes nothing, we can skip
      // that part
      field += rR5s * 2 * nbSite.getQuadrupole() * dr;
      field += -rR7s * dr * dr.dot(nbSite.getQuadrupole() * dr);
    }
  }
  return field;
}

void KSpace::applyStaticShapeField() {
  Eigen::Vector3d shapeField = computeShapeField();

  // Apply the field to the sites
  for (EwdSegment& seg : ewaldSegments_) {
    for (EwdSite& sit : seg) {
      sit.addToStaticField(-shapeField);
    }
  }
}

void KSpace::applyTotalShapeField(PolarSegment& seg) {
  Eigen::Vector3d shapeField = computeShapeField();
  for (PolarSite& site : seg) {
    site.addToBackgroundField(-shapeField);
  }
}

void KSpace::computeIntraMolecularCorrection() {
  for (EwdSegment& seg : ewaldSegments_) {
    for (EwdSite& sit1 : seg) {
      for (EwdSite& sit2 : seg) {
        sit1.addToStaticField(-staticFieldAtBy(sit1, sit2));
      }
    }
  }
}

void KSpace::applySICorrection(PolarSegment& seg) {
  EwdSegment& currentSeg = ewaldSegments_[seg.getId()];
  for (Index i = 0; i < currentSeg.size(); i++) {
    for (Index j = 0; j < currentSeg.size(); j++) {
      seg[i].addToBackgroundField(
          -totalFieldAtBy(currentSeg[i], currentSeg[j]));
    }
  }
}

/**************************************************
 * PRIVATE FUNCTIONS                              *
 **************************************************/

void KSpace::computeScreenedInteraction() {
  double rSqrtPiExp = rSqrtPi * std::exp(-a2 * R2);
  // Note KSpace screening is with erf
  rR1s = std::erf(a1 * R1) * rR1;
  rR3s = rR2 * (rR1s - 2.0 * a1 * rSqrtPiExp);
  rR5s = rR2 * (3.0 * rR3s - (4.0 * a3) * rSqrtPiExp);
  rR7s = rR2 * (5.0 * rR5s - (8.0 * a5) * rSqrtPiExp);
}

void KSpace::computeDistanceVariables(Eigen::Vector3d distVec) {
  dr = distVec;
  R1 = dr.norm();
  R2 = R1 * R1;
  rR1 = 1.0 / R1;
  rR2 = rR1 * rR1;
}

void KSpace::computeTholeVariables(const Eigen::Matrix3d& pol1,
                                   const Eigen::Matrix3d& pol2) {
  thole_u3 =
      (R1 * R2) / std::sqrt((1.0 / 3.0) * (pol1.array() * pol2.array()).sum());

  if (thole * thole_u3 < 40) {
    double thole_exp = std::exp(-thole * thole_u3);
    double thole_u6 = thole_u3 * thole_u3;
    l3 = 1 - thole_exp;
    l5 = 1 - (1 + thole * thole_u3) * thole_exp;
    l7 = 1 - (1 + thole * thole_u3 + (3. / 5.) * thole2 * thole_u6);
    l9 = 1 - (1 + thole * thole_u3 + (18. / 35.) * thole2 * thole_u6 +
              (9. / 35.) * thole3 * thole_u6 * thole_u3) *
                 thole_exp;
  } else {
    l3 = l5 = l7 = l9 = 1.0;
  }
}

Eigen::Matrix3d KSpace::inducedDipoleInteractionAtBy(EwdSite& site,
                                                     const EwdSite& nbSite) {
  computeDistanceVariables(site.getPos() - nbSite.getPos());
  computeScreenedInteraction();
  computeTholeVariables(site.getPolarizationMatrix(),
                        nbSite.getPolarizationMatrix());
  Eigen::Matrix3d interaction = Eigen::Matrix3d::Zero();
  if (R1 < 1e-4) {  // if the same atom, add a correction term
    interaction.diagonal().array() += 4. / 3. * a3 * rSqrtPi;
  } else {  // if not the same atom, add the actual interaction
    interaction.diagonal().array() += l3 * rR3s;
    interaction -= dr * dr.transpose() * l5 * rR5s;
  }
  return interaction;
}

Eigen::Vector3d KSpace::staticFieldAtBy(const EwdSite& site,
                                        const EwdSite& nbSite) {
  computeDistanceVariables(site.getPos() - nbSite.getPos());
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  if (R1 < 1e-2) {
    field += 4.0 / 3.0 * a3 * rSqrtPi * nbSite.getStaticDipole();
  } else {
    // charge
    field += -nbSite.getCharge() * dr * rR3s;
    if (rank > 0) {  // dipole
      field += nbSite.getStaticDipole() * rR3s;
      field += -rR5s * dr * dr.dot(nbSite.getStaticDipole());
      if (rank > 1) {  // quadrupole
        // Using that the trace of a quadrupole contributes nothing, we can skip
        // that part
        field += rR5s * 2 * nbSite.getQuadrupole() * dr;
        field += -rR7s * dr * dr.dot(nbSite.getQuadrupole() * dr);
      }
    }
  }
  return field;
}

Eigen::Vector3d KSpace::totalFieldAtBy(const EwdSite& site,
                                       const EwdSite& nbSite) {
  computeDistanceVariables(site.getPos() - nbSite.getPos());
  computeScreenedInteraction();

  Eigen::Vector3d field = Eigen::Vector3d::Zero();
  Index rank = nbSite.getRank();
  if (R1 < 1e-2) {
    field += 4.0 / 3.0 * a3 * rSqrtPi * nbSite.getStaticDipole();
  } else {
    // charge
    field += -nbSite.getCharge() * dr * rR3s;
    if (rank > 0) {  // dipole
      field += nbSite.getTotalDipole() * rR3s;
      field += -rR5s * dr * dr.dot(nbSite.getTotalDipole());
      if (rank > 1) {  // quadrupole
        // Using that the trace of a quadrupole contributes nothing, we can skip
        // that part
        field += rR5s * 2 * nbSite.getQuadrupole() * dr;
        field += -rR7s * dr * dr.dot(nbSite.getQuadrupole() * dr);
      }
    }
  }
  return field;
}

std::complex<double> KSpace::computeSk(const Eigen::Vector3d& kvector) const {
  std::complex<double> sk(0.0, 0.0);
  for (const EwdSegment& seg : ewaldSegments_) {
    for (const EwdSite& site : seg) {
      std::complex<double> generalizedCharge =
          (site.getCharge() + ii * kvector.dot(site.getTotalDipole()) -
           kvector.dot(site.getQuadrupole() * kvector));
      std::complex<double> expKR = std::exp(ii * kvector.dot(site.getPos()));
      sk += generalizedCharge * expKR;
    }
  }
  return sk;
}

double KSpace::computeAk(const Eigen::Vector3d& kvector) const {
  /* Compute the A_k factor, i.e. k^(-2) exp(-k^2/(4\alpha^2)) */
  double k_squared = kvector.squaredNorm();
  return std::exp(-k_squared / (4 * a2)) / k_squared;
}

void KSpace::computeKVectors() {
  for (Index ix = -max_K[0]; ix <= max_K[0]; ++ix) {
    for (Index iy = -max_K[1]; iy <= max_K[1]; ++iy) {
      for (Index iz = -max_K[2]; iz <= max_K[2]; ++iz) {
        if (ix == 0 && iy == 0 && iz == 0) {
          continue;
        }
        Eigen::Vector3d kvector = unit_cell_.getKVector(ix, iy, iz);
        kvector_list_.push_back(KVector(kvector));
      }
    }
  }
  std::sort(kvector_list_.begin(), kvector_list_.end());
  for (Index i = 0; i < static_cast<Index>(kvector_list_.size()); ++i) {
    KVector& kvec = kvector_list_[i];
    kvec.setSk(computeSk(kvec.getVector()));
    kvec.setAk(computeAk(kvec.getVector()));
  }
}

}  // namespace xtp
}  // namespace votca