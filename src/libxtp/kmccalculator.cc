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
  _seed = options.get(".seed").as<Index>();

  _numberofcarriers = options.get(".numberofcarriers").as<Index>();
  _injection_name = options.get(".injectionpattern").as<std::string>();
  _maxrealtime = options.get(".maxrealtime").as<double>();
  _trajectoryfile = options.get(".trajectoryfile").as<std::string>();
  _temperature = options.get(".temperature").as<double>();

  _temperature *= (tools::conv::kB * tools::conv::ev2hrt);
  _occfile = options.get(".occfile").as<std::string>();
  _ratefile = options.get(".ratefile").as<std::string>();

  _injectionmethod = options.get(".injectionmethod").as<std::string>();
}

void KMCCalculator::LoadGraph(Topology& top) {

  std::vector<Segment>& segs = top.Segments();

  if (segs.size() < 1) {
    throw std::runtime_error("Your state file contains no segments!");
  }
  _nodes.reserve(segs.size());
  for (Segment& seg : segs) {
    bool injectable = false;
    if (tools::wildcmp(_injection_name, seg.getType())) {
      injectable = true;
    }
    _nodes.push_back(GNode(seg, _carriertype, injectable));
  }

  QMNBList& nblist = top.NBList();
  if (nblist.size() < 1) {
    throw std::runtime_error("neighborlist contains no pairs!");
  }
  if (_temperature <= 0) {
    throw std::runtime_error(
        "Your Temperature is negative or zero, please specify the temperature "
        "in Kelvin.");
  }

  Rate_Engine rate_engine(_temperature, _field);
  XTP_LOG(Log::error, _log) << "\nCalculating initial rates." << std::flush;
  XTP_LOG(Log::error, _log) << rate_engine << std::flush;
  XTP_LOG(Log::error, _log)
      << "    carriertype: " << _carriertype.ToLongString() << std::flush;

  for (const QMPair* pair : nblist) {
    Rate_Engine::PairRates rates = rate_engine.Rate(*pair, _carriertype);
    _nodes[pair->Seg1()->getId()].AddEventfromQmPair(*pair, _nodes,
                                                     rates.rate12);
    _nodes[pair->Seg2()->getId()].AddEventfromQmPair(*pair, _nodes,
                                                     rates.rate21);
  }
  _RandomVariable.setMaxInt(Index(_nodes.size()));
  XTP_LOG(Log::error, _log) << "    Rates for " << _nodes.size()
                            << " sites are computed." << std::flush;
  WriteRatestoFile(_ratefile, nblist);

  Index events = 0;
  Index max = std::numeric_limits<Index>::min();
  Index min = std::numeric_limits<Index>::max();
  double minlength = std::numeric_limits<double>::max();
  double maxlength = 0;
  for (const auto& node : _nodes) {

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
      XTP_LOG(Log::info, _log)
          << "Node " << node.getId() << " has 0 jumps" << std::flush;
    } else if (size < min) {
      min = size;
    } else if (size > max) {
      max = size;
    }
  }
  double avg = double(events) / double(_nodes.size());
  double deviation = 0.0;
  for (const auto& node : _nodes) {
    double size = double(node.Events().size());
    deviation += (size - avg) * (size - avg);
  }
  deviation = std::sqrt(deviation / double(_nodes.size()));

  XTP_LOG(Log::error, _log)
      << "Nblist has " << nblist.size() << " pairs. Nodes contain " << events
      << " jump events\n"
      << "with avg=" << avg << " std=" << deviation << " max=" << max
      << " min=" << min << " jumps per site\n"
      << "Minimum jumpdistance =" << minlength * tools::conv::bohr2nm
      << " nm Maximum distance =" << maxlength * tools::conv::bohr2nm << " nm\n"
      << std::flush;
  double conv = std::pow(tools::conv::bohr2nm, 3);
  XTP_LOG(Log::error, _log)
      << "spatial carrier density: "
      << double(_numberofcarriers) / (top.BoxVolume() * conv) << " nm^-3"
      << std::flush;

  for (auto& node : _nodes) {
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

  XTP_LOG(Log::error, _log) << "looking for injectable nodes..." << std::flush;
  for (Index i = 0; i < _numberofcarriers; i++) {
    Chargecarrier newCharge(i);
    RandomlyAssignCarriertoSite(newCharge);

    XTP_LOG(Log::error, _log)
        << "starting position for charge " << i << ": segment "
        << newCharge.getCurrentNodeId() << std::flush;
    _carriers.push_back(newCharge);
  }
  return;
}

void KMCCalculator::RandomlyAssignCarriertoSite(Chargecarrier& Charge) {
  Index nodeId_guess = -1;
  do {
    nodeId_guess = _RandomVariable.rand_uniform_int();
  } while (_nodes[nodeId_guess].isOccupied() ||
           _nodes[nodeId_guess].isInjectable() ==
               false);  // maybe already occupied? or maybe not injectable?
  if (Charge.hasNode()) {
    Charge.ReleaseNode();
  }
  Charge.settoNote(&_nodes[nodeId_guess]);

  return;
}

double KMCCalculator::Promotetime(double cumulated_rate) {
  double dt = 0;
  double rand_u = 1 - _RandomVariable.rand_uniform();
  dt = -1 / cumulated_rate * std::log(rand_u);
  return dt;
}

const GLink& KMCCalculator::ChooseHoppingDest(const GNode& node) {
  double u = 1 - _RandomVariable.rand_uniform();
  return *(node.findHoppingDestination(u));
}

Chargecarrier* KMCCalculator::ChooseAffectedCarrier(double cumulated_rate) {
  if (_carriers.size() == 1) {
    return &_carriers[0];
  }
  Chargecarrier* carrier = nullptr;
  double u = 1 - _RandomVariable.rand_uniform();
  for (Index i = 0; i < _numberofcarriers; i++) {
    u -= _carriers[i].getCurrentEscapeRate() / cumulated_rate;
    if (u <= 0 || i == _numberofcarriers - 1) {
      carrier = &_carriers[i];
      break;
    }
  }
  return carrier;
}
void KMCCalculator::WriteRatestoFile(std::string filename,
                                     const QMNBList& nblist) {
  XTP_LOG(Log::error, _log)
      << "\nRates are written to " << filename << std::flush;
  fstream ratefs;
  ratefs.open(filename, fstream::out);
  ratefs << "#PairID,SiteID1,SiteID2, ,rate12[1/s],rate21[1/s] at "
         << _temperature * tools::conv::hrt2ev / tools::conv::kB
         << "K for carrier:" << _carriertype.ToString() << endl;

  Rate_Engine rate_engine(_temperature, _field);
  for (const QMPair* pair : nblist) {
    Rate_Engine::PairRates rates = rate_engine.Rate(*pair, _carriertype);
    ratefs << pair->getId() << " " << pair->Seg1()->getId() << " "
           << pair->Seg2()->getId() << " " << rates.rate12 << " "
           << rates.rate21 << "\n";
  }
  ratefs << std::flush;
  ratefs.close();
}

void KMCCalculator::WriteOccupationtoFile(double simtime,
                                          std::string filename) {
  XTP_LOG(Log::error, _log)
      << "\nOccupations are written to " << filename << std::flush;
  fstream probs;
  probs.open(filename, fstream::out);
  probs << "#SiteID, Occupation prob at "
        << _temperature * tools::conv::hrt2ev / tools::conv::kB
        << "K for carrier:" << _carriertype.ToString() << endl;
  for (const GNode& node : _nodes) {
    double occupationprobability = node.OccupationTime() / simtime;
    probs << node.getId() << "\t" << occupationprobability << endl;
  }
  probs.close();
}

}  // namespace xtp
}  // namespace votca
