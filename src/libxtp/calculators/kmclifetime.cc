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
#include <exception>
#include <fstream>

// Third party includes
#include <boost/format.hpp>
#include <boost/progress.hpp>

// VOTCA includes
#include <votca/tools/constants.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/gnode.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/topology.h"

// Local private VOTCA includes
#include "kmclifetime.h"

namespace votca {
namespace xtp {

void KMCLifetime::ParseSpecificOptions(const tools::Property& options) {

  insertions_ = options.get(".numberofinsertions").as<unsigned long>();
  lifetimefile_ = options.get(".lifetimefile").as<std::string>();

  probfile_ = options.ifExistsReturnElseReturnDefault<std::string>(
      ".decayprobfile", "");

  const tools::Property& carrier_options = options.get(".carrierenergy");
  do_carrierenergy_ = carrier_options.get(".run").as<bool>();
  energy_outputfile_ = carrier_options.get(".outputfile").as<std::string>();
  alpha_ = carrier_options.get(".alpha").as<double>();
  outputsteps_ = carrier_options.get(".outputsteps").as<unsigned long>();

  log_.setReportLevel(Log::current_level);
  log_.setCommonPreface("\n ...");

  carriertype_ = QMStateType(QMStateType::Singlet);
  XTP_LOG(Log::error, log_)
      << "carrier type:" << carriertype_.ToLongString() << std::flush;
  field_ = Eigen::Vector3d::Zero();
}

void KMCLifetime::WriteDecayProbability(std::string filename) {

  Eigen::VectorXd outrates = Eigen::VectorXd::Zero(nodes_.size());
  Eigen::VectorXd inrates = Eigen::VectorXd::Zero(nodes_.size());
  Eigen::VectorXd decayrates = Eigen::VectorXd::Ones(nodes_.size());

  for (unsigned i = 0; i < nodes_.size(); i++) {
    GNode& node = nodes_[i];
    if (node.canDecay()) {
      for (const GLink& event : node.Events()) {
        if (event.isDecayEvent()) {
          decayrates[i] = event.getRate();
        } else {
          inrates[event.getDestination()->getId()] += event.getRate();
          outrates[i] += event.getRate();
        }
      }
    }
  }
  outrates = decayrates.cwiseQuotient(outrates + decayrates);
  inrates = decayrates.cwiseQuotient(inrates + decayrates);

  std::fstream probs;
  probs.open(filename, std::fstream::out);
  if (!probs.is_open()) {
    std::string error_msg = "Unable to write to file " + filename;
    throw std::runtime_error(error_msg);
  }
  probs << "#SiteID, Relative Prob outgoing, Relative Prob ingoing\n";
  for (unsigned i = 0; i < nodes_.size(); i++) {
    probs << nodes_[i].getId() << " " << outrates[i] << " " << inrates[i]
          << std::endl;
  }
  probs.close();
  return;
}

void KMCLifetime::ReadLifetimeFile(std::string filename) {
  tools::Property xml;
  xml.LoadFromXML(filename);
  std::vector<tools::Property*> jobProps = xml.Select("lifetimes.site");
  if (jobProps.size() != nodes_.size()) {
    throw std::runtime_error(
        (boost::format("The number of sites in the topology: %i does not match "
                       "the number in the lifetimefile: %i") %
         nodes_.size() % jobProps.size())
            .str());
  }

  for (tools::Property* prop : jobProps) {
    Index site_id = prop->getAttribute<Index>("id");
    double lifetime = boost::lexical_cast<double>(prop->value());
    bool check = false;
    for (auto& node : nodes_) {
      if (node.getId() == site_id && !(node.canDecay())) {
        node.AddDecayEvent(1.0 / lifetime);
        check = true;
        break;
      } else if (node.getId() == site_id && node.canDecay()) {
        throw std::runtime_error(
            (boost::format("Node %i appears twice in your list") % site_id)
                .str());
      }
    }
    if (!check) {
      throw std::runtime_error(
          (boost::format("Site from file with id: %i not found in state file") %
           site_id)
              .str());
    }
  }
  for (auto& node : nodes_) {
    node.InitEscapeRate();
    node.MakeHuffTree();
  }
  return;
}

void KMCLifetime::WriteToTraj(std::fstream& traj, unsigned long insertioncount,
                              double simtime,
                              const Chargecarrier& affectedcarrier) const {
  const Eigen::Vector3d& dr_travelled = affectedcarrier.get_dRtravelled();
  traj << simtime << "\t" << insertioncount << "\t" << affectedcarrier.getId()
       << "\t" << affectedcarrier.getLifetime() << "\t"
       << affectedcarrier.getSteps() << "\t"
       << affectedcarrier.getCurrentNodeId() + 1 << "\t"
       << dr_travelled.x() * tools::conv::bohr2nm << "\t"
       << dr_travelled.y() * tools::conv::bohr2nm << "\t"
       << dr_travelled.z() * tools::conv::bohr2nm << std::endl;
}

void KMCLifetime::RunVSSM() {

  std::chrono::time_point<std::chrono::system_clock> realtime_start =
      std::chrono::system_clock::now();
  XTP_LOG(Log::error, log_)
      << "\nAlgorithm: VSSM for Multiple Charges with finite Lifetime\n"
         "number of charges: "
      << numberofcarriers_ << "\nnumber of nodes: " << nodes_.size()
      << std::flush;

  if (numberofcarriers_ > Index(nodes_.size())) {
    throw std::runtime_error(
        "ERROR in kmclifetime: specified number of charges is greater than the "
        "number of nodes. This conflicts with single occupation.");
  }

  std::fstream traj;
  std::fstream energyfile;

  XTP_LOG(Log::error, log_)
      << "Writing trajectory to " << trajectoryfile_ << "." << std::flush;

  traj.open(trajectoryfile_, std::fstream::out);
  if (!traj.is_open()) {
    std::string error_msg = "Unable to write to file " + trajectoryfile_;
    throw std::runtime_error(error_msg);
  }

  traj << "#Simtime [s]\t Insertion\t Carrier ID\t Lifetime[s]\tSteps\t Last "
          "Segment\t x_travelled[nm]\t y_travelled[nm]\t z_travelled[nm]\n";

  if (do_carrierenergy_) {

    XTP_LOG(Log::error, log_)
        << "Tracking the energy of one charge carrier and exponential average "
           "with alpha="
        << alpha_ << " to " << energy_outputfile_ << std::flush;
    energyfile.open(energy_outputfile_, std::fstream::out);
    if (!energyfile.is_open()) {
      std::string error_msg = "Unable to write to file " + energy_outputfile_;
      throw std::runtime_error(error_msg);
    }

    energyfile << "Simtime [s]\tSteps\tCarrier ID\tEnergy_a=" << alpha_
               << "[eV]\n";
  }

  // Injection
  XTP_LOG(Log::error, log_)
      << "\ninjection method: " << injectionmethod_ << std::flush;

  RandomlyCreateCharges();

  unsigned long insertioncount = 0;
  unsigned long step = 0;
  double simtime = 0.0;

  std::vector<GNode*> forbiddennodes;
  std::vector<GNode*> forbiddendests;

  time_t now = time(nullptr);
  tm* localtm = localtime(&now);
  XTP_LOG(Log::error, log_)
      << "Run started at " << asctime(localtm) << std::flush;

  double avlifetime = 0.0;
  double meanfreepath = 0.0;
  Eigen::Vector3d difflength_squared = Eigen::Vector3d::Zero();

  double avgenergy = carriers_[0].getCurrentEnergy();
  Index carrieridold = carriers_[0].getId();

  while (insertioncount < insertions_) {
    std::chrono::duration<double> elapsed_time =
        std::chrono::system_clock::now() - realtime_start;
    if (elapsed_time.count() > (maxrealtime_ * 60. * 60.)) {
      XTP_LOG(Log::error, log_)
          << "\nReal time limit of " << maxrealtime_ << " hours ("
          << Index(maxrealtime_ * 60 * 60 + 0.5)
          << " seconds) has been reached. Stopping here.\n"
          << std::flush;
      break;
    }

    double cumulated_rate = 0;

    for (const auto& carrier : carriers_) {
      cumulated_rate += carrier.getCurrentEscapeRate();
    }
    if (cumulated_rate == 0) {  // this should not happen: no possible jumps
                                // defined for a node
      throw std::runtime_error(
          "ERROR in kmclifetime: Incorrect rates in the database file. All the "
          "escape rates for the current setting are 0.");
    }
    // go forward in time
    double dt = Promotetime(cumulated_rate);

    if (do_carrierenergy_) {
      bool print = false;
      if (carriers_[0].getId() > carrieridold) {
        avgenergy = carriers_[0].getCurrentEnergy();
        print = true;
        carrieridold = carriers_[0].getId();
      } else if (step % outputsteps_ == 0) {
        avgenergy =
            alpha_ * carriers_[0].getCurrentEnergy() + (1 - alpha_) * avgenergy;
        print = true;
      }
      if (print) {
        energyfile << simtime << "\t" << step << "\t" << carriers_[0].getId()
                   << "\t" << avgenergy * tools::conv::hrt2ev << std::endl;
      }
    }

    simtime += dt;
    step++;
    for (auto& carrier : carriers_) {
      carrier.updateLifetime(dt);
      carrier.updateSteps(1);
      carrier.updateOccupationtime(dt);
    }

    ResetForbiddenlist(forbiddennodes);
    bool secondlevel = true;
    while (secondlevel) {

      // determine which carrier will escape
      GNode* newnode = nullptr;
      Chargecarrier* affectedcarrier = ChooseAffectedCarrier(cumulated_rate);

      if (CheckForbidden(affectedcarrier->getCurrentNode(), forbiddennodes)) {
        continue;
      }

      // determine where it will jump to
      ResetForbiddenlist(forbiddendests);
      boost::progress_display progress(insertions_);

      while (true) {
        // LEVEL 2

        newnode = nullptr;
        const GLink& event =
            ChooseHoppingDest(affectedcarrier->getCurrentNode());

        if (event.isDecayEvent()) {
          const Eigen::Vector3d& dr_travelled =
              affectedcarrier->get_dRtravelled();
          avlifetime += affectedcarrier->getLifetime();
          meanfreepath += dr_travelled.norm();
          difflength_squared += dr_travelled.cwiseAbs2();
          WriteToTraj(traj, insertioncount, simtime, *affectedcarrier);
          ++progress;
          RandomlyAssignCarriertoSite(*affectedcarrier);
          affectedcarrier->resetCarrier();
          insertioncount++;
          affectedcarrier->setId(numberofcarriers_ - 1 + insertioncount);
          secondlevel = false;
          break;
        } else {
          newnode = event.getDestination();
        }

        // check after the event if this was allowed
        if (CheckForbidden(*newnode, forbiddendests)) {
          continue;
        }

        // if the new segment is unoccupied: jump; if not: add to forbidden list
        // and choose new hopping destination
        if (newnode->isOccupied()) {
          if (CheckSurrounded(affectedcarrier->getCurrentNode(),
                              forbiddendests)) {
            AddtoForbiddenlist(affectedcarrier->getCurrentNode(),
                               forbiddennodes);
            break;  // select new escape node (ends level 2 but without setting
                    // level1step to 1)
          }
          AddtoForbiddenlist(*newnode, forbiddendests);
          continue;  // select new destination
        } else {
          affectedcarrier->jumpAccordingEvent(event);
          secondlevel = false;

          break;  // this ends LEVEL 2 , so that the time is updated and the
                  // next MC step started
        }

        // END LEVEL 2
      }
      // END LEVEL 1
    }
  }

  XTP_LOG(Log::error, log_)
      << "\nTotal runtime:\t\t\t\t\t" << simtime
      << " s\n"
         "Total KMC steps:\t\t\t\t"
      << step << "\nAverage lifetime:\t\t\t\t"
      << avlifetime / double(insertioncount) << " s\n"
      << "Mean freepath\t l=<|r_x-r_o|> :\t\t"
      << (meanfreepath * tools::conv::bohr2nm / double(insertioncount))
      << " nm\n"
      << "Average diffusionlength\t d=sqrt(<(r_x-r_o)^2>)\t"
      << std::sqrt(difflength_squared.norm() / double(insertioncount)) *
             tools::conv::bohr2nm
      << " nm\n"
      << std::flush;

  WriteOccupationtoFile(simtime, occfile_);
  traj.close();
  if (do_carrierenergy_) {
    energyfile.close();
  }
  return;
}

bool KMCLifetime::Evaluate(Topology& top) {
  XTP_LOG(Log::error, log_) << "\n-----------------------------------"
                               "\n      KMCLIFETIME started"
                               "\n-----------------------------------\n"
                            << std::flush;

  XTP_LOG(Log::info, log_) << "\nInitialising random number generator"
                           << std::flush;

  RandomVariable_.init(seed_);
  LoadGraph(top);
  ReadLifetimeFile(lifetimefile_);

  if (!probfile_.empty()) {
    WriteDecayProbability(probfile_);
  }
  RunVSSM();

  time_t now = time(nullptr);
  tm* localtm = localtime(&now);

  XTP_LOG(Log::info, log_) << " KMCLIFETIME finished at:" << asctime(localtm)
                           << std::flush;
  std::cout << log_;
  return true;
}

}  // namespace xtp
}  // namespace votca
