/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include <locale>

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/gnode.h"
#include "votca/xtp/kmccalculator.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/rate_engine.h"
#include "votca/xtp/topology.h"

using namespace std;

namespace votca {
namespace xtp {

void KMCCalculator::ParseCommonOptions(const tools::Property& options) {
  seed_ = options.get(".seed").as<Index>();

  numberofcarriers_ = options.get(".numberofcarriers").as<Index>();
  injection_name_ = options.get(".injectionpattern").as<std::string>();
  maxrealtime_ = options.get(".maxrealtime").as<double>();
  trajectoryfile_ = options.get(".trajectoryfile").as<std::string>();
  temperature_ = options.get(".temperature").as<double>();

  temperature_ *= (tools::conv::kB * tools::conv::ev2hrt);
  occfile_ = options.get(".occfile").as<std::string>();
  ratefile_ = options.get(".ratefile").as<std::string>();

  injectionmethod_ = options.get(".injectionmethod").as<std::string>();
}

void KMCCalculator::LoadGraph(Topology& top) {

  std::vector<Segment>& segs = top.Segments();

  if (segs.size() < 1) {
    throw std::runtime_error("Your state file contains no segments!");
  }
  nodes_.reserve(segs.size());
  for (Segment& seg : segs) {
    bool injectable = false;
    if (tools::wildcmp(injection_name_, seg.getType())) {
      injectable = true;
    }
    nodes_.push_back(GNode(seg, carriertype_, injectable));
  }

  QMNBList& nblist = top.NBList();
  if (nblist.size() < 1) {
    throw std::runtime_error("neighborlist contains no pairs!");
  }
  if (temperature_ <= 0) {
    throw std::runtime_error(
        "Your Temperature is negative or zero, please specify the temperature "
        "in Kelvin.");
  }

  Rate_Engine rate_engine(temperature_, field_);
  XTP_LOG(Log::error, log_) << "\nCalculating initial rates." << std::flush;
  XTP_LOG(Log::error, log_) << rate_engine << std::flush;
  XTP_LOG(Log::error, log_)
      << "    carriertype: " << carriertype_.ToLongString() << std::flush;

  for (const QMPair* pair : nblist) {
    Rate_Engine::PairRates rates = rate_engine.Rate(*pair, carriertype_);
    nodes_[pair->Seg1()->getId()].AddEventfromQmPair(*pair, nodes_,
                                                     rates.rate12);
    nodes_[pair->Seg2()->getId()].AddEventfromQmPair(*pair, nodes_,
                                                     rates.rate21);
  }
  RandomVariable_.setMaxInt(Index(nodes_.size()));
  XTP_LOG(Log::error, log_) << "    Rates for " << nodes_.size()
                            << " sites are computed." << std::flush;
  WriteRatestoFile(ratefile_, nblist);

  Index events = 0;
  Index max = std::numeric_limits<Index>::min();
  Index min = std::numeric_limits<Index>::max();
  double minlength = std::numeric_limits<double>::max();
  double maxlength = 0;
  for (const auto& node : nodes_) {

    Index size = Index(node.Events().size());
    for (const auto& event : node.Events()) {
      if (event.isDecayEvent()) {
        continue;
      }
      double dist = event.getDeltaR().norm();
      if (dist > maxlength) {
        maxlength = dist;
      } else if (dist < minlength) {
        minlength = dist;
      }
    }

    events += size;
    if (size == 0) {
      XTP_LOG(Log::info, log_)
          << "Node " << node.getId() << " has 0 jumps" << std::flush;
    } else if (size < min) {
      min = size;
    } else if (size > max) {
      max = size;
    }
  }
  double avg = double(events) / double(nodes_.size());
  double deviation = 0.0;
  for (const auto& node : nodes_) {
    double size = double(node.Events().size());
    deviation += (size - avg) * (size - avg);
  }
  deviation = std::sqrt(deviation / double(nodes_.size()));

  XTP_LOG(Log::error, log_)
      << "Nblist has " << nblist.size() << " pairs. Nodes contain " << events
      << " jump events\n"
      << "with avg=" << avg << " std=" << deviation << " max=" << max
      << " min=" << min << " jumps per site\n"
      << "Minimum jumpdistance =" << minlength * tools::conv::bohr2nm
      << " nm Maximum distance =" << maxlength * tools::conv::bohr2nm << " nm\n"
      << std::flush;
  double conv = std::pow(tools::conv::bohr2nm, 3);
  XTP_LOG(Log::error, log_)
      << "spatial carrier density: "
      << double(numberofcarriers_) / (top.BoxVolume() * conv) << " nm^-3"
      << std::flush;

  for (auto& node : nodes_) {
    node.InitEscapeRate();
    node.MakeHuffTree();
  }
  return;
}

void KMCCalculator::ResetForbiddenlist(
    std::vector<GNode*>& forbiddenlist) const {
  forbiddenlist.clear();
  return;
}

void KMCCalculator::AddtoForbiddenlist(
    GNode& node, std::vector<GNode*>& forbiddenlist) const {
  forbiddenlist.push_back(&node);
  return;
}

bool KMCCalculator::CheckForbidden(
    const GNode& node, const std::vector<GNode*>& forbiddenlist) const {
  bool forbidden = false;
  for (const GNode* fnode : forbiddenlist) {
    if (&node == fnode) {
      forbidden = true;
      break;
    }
  }
  return forbidden;
}

bool KMCCalculator::CheckSurrounded(
    const GNode& node, const std::vector<GNode*>& forbiddendests) const {
  bool surrounded = true;
  for (const auto& event : node.Events()) {
    bool thisevent_possible = true;
    for (const GNode* fnode : forbiddendests) {
      if (event.getDestination() == fnode) {
        thisevent_possible = false;
        break;
      }
    }
    if (thisevent_possible == true) {
      surrounded = false;
      break;
    }
  }
  return surrounded;
}

void KMCCalculator::RandomlyCreateCharges() {

  XTP_LOG(Log::error, log_) << "looking for injectable nodes..." << std::flush;
  for (Index i = 0; i < numberofcarriers_; i++) {
    Chargecarrier newCharge(i);
    RandomlyAssignCarriertoSite(newCharge);

    XTP_LOG(Log::error, log_)
        << "starting position for charge " << i << ": segment "
        << newCharge.getCurrentNodeId() << std::flush;
    carriers_.push_back(newCharge);
  }
  return;
}

void KMCCalculator::RandomlyAssignCarriertoSite(Chargecarrier& Charge) {
  Index nodeId_guess = -1;
  do {
    nodeId_guess = RandomVariable_.rand_uniform_int();
  } while (nodes_[nodeId_guess].isOccupied() ||
           nodes_[nodeId_guess].isInjectable() ==
               false);  // maybe already occupied? or maybe not injectable?
  if (Charge.hasNode()) {
    Charge.ReleaseNode();
  }
  Charge.settoNote(&nodes_[nodeId_guess]);

  return;
}

double KMCCalculator::Promotetime(double cumulated_rate) {
  double dt = 0;
  double rand_u = 1 - RandomVariable_.rand_uniform();
  dt = -1 / cumulated_rate * std::log(rand_u);
  return dt;
}

const GLink& KMCCalculator::ChooseHoppingDest(const GNode& node) {
  double u = 1 - RandomVariable_.rand_uniform();
  return *(node.findHoppingDestination(u));
}

Chargecarrier* KMCCalculator::ChooseAffectedCarrier(double cumulated_rate) {
  if (carriers_.size() == 1) {
    return &carriers_[0];
  }
  Chargecarrier* carrier = nullptr;
  double u = 1 - RandomVariable_.rand_uniform();
  for (Index i = 0; i < numberofcarriers_; i++) {
    u -= carriers_[i].getCurrentEscapeRate() / cumulated_rate;
    if (u <= 0 || i == numberofcarriers_ - 1) {
      carrier = &carriers_[i];
      break;
    }
  }
  return carrier;
}
void KMCCalculator::WriteRatestoFile(std::string filename,
                                     const QMNBList& nblist) {
  XTP_LOG(Log::error, log_)
      << "\nRates are written to " << filename << std::flush;
  fstream ratefs;
  ratefs.open(filename, fstream::out);
  ratefs << "#PairID,SiteID1,SiteID2, ,rate12[1/s],rate21[1/s] at "
         << temperature_ * tools::conv::hrt2ev / tools::conv::kB
         << "K for carrier:" << carriertype_.ToString() << endl;

  Rate_Engine rate_engine(temperature_, field_);
  for (const QMPair* pair : nblist) {
    Rate_Engine::PairRates rates = rate_engine.Rate(*pair, carriertype_);
    ratefs << pair->getId() << " " << pair->Seg1()->getId() << " "
           << pair->Seg2()->getId() << " " << rates.rate12 << " "
           << rates.rate21 << "\n";
  }
  ratefs << std::flush;
  ratefs.close();
}

void KMCCalculator::WriteOccupationtoFile(double simtime,
                                          std::string filename) {
  XTP_LOG(Log::error, log_)
      << "\nOccupations are written to " << filename << std::flush;
  fstream probs;
  probs.open(filename, fstream::out);
  probs << "#SiteID, Occupation prob at "
        << temperature_ * tools::conv::hrt2ev / tools::conv::kB
        << "K for carrier:" << carriertype_.ToString() << endl;
  for (const GNode& node : nodes_) {
    double occupationprobability = node.OccupationTime() / simtime;
    probs << node.getId() << "\t" << occupationprobability << endl;
  }
  probs.close();
}

}  // namespace xtp
}  // namespace votca
