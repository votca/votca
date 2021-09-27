/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// Standard includes
#include <fstream>
#include <iostream>
#include <stdexcept>

// VOTCA includes
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/application.h>
#include <votca/tools/filesystem.h>
#include <votca/tools/globals.h>

// Local VOTCA includes
#include "votca/xtp/md2qmengine.h"
#include "votca/xtp/statesaver.h"
#include "votca/xtp/topology.h"
#include "votca/xtp/version.h"

using namespace std;

namespace CSG = votca::csg;
namespace XTP = votca::xtp;
namespace TOOLS = votca::tools;

class XtpMap : public TOOLS::Application {

 public:
  string ProgramName() override { return "xtp_map"; }
  void HelpText(ostream& out) override {
    out << "Generates QM|MD topology" << endl;
  }
  void ShowHelpText(std::ostream& out) override;

  void Initialize() override;
  bool EvaluateOptions() override;
  void Run() override;

 protected:
};

namespace propt = boost::program_options;

void XtpMap::Initialize() {

  CSG::TrajectoryWriter::RegisterPlugins();
  CSG::TrajectoryReader::RegisterPlugins();
  CSG::TopologyReader::RegisterPlugins();

  AddProgramOptions()("topology,t", propt::value<string>(), "  topology");
  AddProgramOptions()("coordinates,c", propt::value<string>(),
                      "  coordinates or trajectory");
  AddProgramOptions()("segments,s", propt::value<string>(),
                      "  definition of segments and fragments");
  AddProgramOptions()("makesegments,m", "  write out a skeleton segments file");
  AddProgramOptions()("file,f", propt::value<string>(), "  state file");
  AddProgramOptions()("first-frame,i",
                      propt::value<votca::Index>()->default_value(0),
                      "  start from this frame");
  AddProgramOptions()("begin,b", propt::value<double>()->default_value(0.0),
                      "  start time in simulation");
  AddProgramOptions()("nframes,n",
                      propt::value<votca::Index>()->default_value(1),
                      "  number of frames to process");
}

bool XtpMap::EvaluateOptions() {

  CheckRequired("topology", "Missing topology file");
  CheckRequired("segments", "Missing segment definition file");
  CheckRequired("coordinates", "Missing trajectory input");
  if (!(OptionsMap().count("makesegments"))) {
    CheckRequired("file", "Missing state file");
  }
  return 1;
}

void XtpMap::Run() {

  std::string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  XTP::HelpTextHeader(name);

  // ++++++++++++++++++++++++++++ //
  // Create MD topology from file //
  // ++++++++++++++++++++++++++++ //

  // Create topology reader
  string topfile = OptionsMap()["topology"].as<string>();
  std::unique_ptr<CSG::TopologyReader> topread =
      CSG::TopReaderFactory().Create(topfile);

  if (topread == nullptr) {
    throw runtime_error(string("Input format not supported: ") +
                        OptionsMap()["topology"].as<string>());
  }
  CSG::Topology mdtopol;
  topread->ReadTopology(topfile, mdtopol);
  if (votca::Log::verbose()) {
    cout << "Read MD topology from " << topfile << ": Found "
         << mdtopol.BeadCount() << " atoms in " << mdtopol.MoleculeCount()
         << " molecules. " << endl;
  }

  // ++++++++++++++++++++++++++++++ //
  // Create MD trajectory from file //
  // ++++++++++++++++++++++++++++++ //

  // Create trajectory reader and initialize
  string trjfile = OptionsMap()["coordinates"].as<string>();
  std::unique_ptr<CSG::TrajectoryReader> trjread =
      CSG::TrjReaderFactory().Create(trjfile);

  if (trjread == nullptr) {
    throw runtime_error(string("Input format not supported: ") +
                        OptionsMap()["coordinates"].as<string>());
  }
  trjread->Open(trjfile);
  trjread->FirstFrame(mdtopol);

  string mapfile = OptionsMap()["segments"].as<string>();
  if (OptionsMap().count("makesegments")) {
    if (TOOLS::filesystem::FileExists(mapfile)) {
      cout << endl
           << "xtp_map : map file '" << mapfile
           << "' already in use. Delete the current mapfile or specify a "
              "different name."
           << endl;
      return;
    }

    cout << " Writing template mapfile to " << mapfile << std::endl;

    TOOLS::Property mapfile_prop("topology", "", "");
    TOOLS::Property& molecules = mapfile_prop.add("molecules", "");

    std::map<std::string, const CSG::Molecule*> firstmolecule;

    std::map<std::string, votca::Index> molecule_names;
    for (const CSG::Molecule& mol : mdtopol.Molecules()) {
      if (!molecule_names.count(mol.getName())) {
        firstmolecule[mol.getName()] = &mol;
      }
      molecule_names[mol.getName()]++;
    }
    for (const auto& mol : molecule_names) {
      std::cout << "Found " << mol.second << " with name " << mol.first
                << std::endl;
    }
    for (const auto& mol : molecule_names) {
      TOOLS::Property& molecule = molecules.add("molecule", "");
      molecule.add("mdname", mol.first);
      TOOLS::Property& segments = molecule.add("segments", "");
      TOOLS::Property& segment = segments.add("segment", "");
      segment.add("name", "UPTOYOU_BUTUNIQUE");
      segment.add("qmcoords_n", "XYZFILE_GROUNDSTATE");
      segment.add("multipoles_n", "MPSFILE_GROUNDSTATE");
      segment.add("map2md", "WANTTOMAPTOMDGEOMETRY");
      segment.add("U_xX_nN_h", "REORG1_hole");
      segment.add("U_nX_nN_h", "REORG2_hole");
      segment.add("U_xN_xX_h", "REORG3_hole");
      TOOLS::Property& fragments = segment.add("fragments", "");
      TOOLS::Property& fragment = fragments.add("fragment", "");
      std::string atomnames = "";
      const CSG::Molecule* csgmol = firstmolecule[mol.first];
      std::vector<const CSG::Bead*> sortedbeads;
      sortedbeads.reserve(csgmol->BeadCount());
      for (const CSG::Bead* bead : csgmol->Beads()) {
        sortedbeads.push_back(bead);
      }
      std::sort(sortedbeads.begin(), sortedbeads.end(),
                [&](const CSG::Bead* b1, const CSG::Bead* b2) {
                  return b1->getId() < b2->getId();
                });

      for (const CSG::Bead* bead : sortedbeads) {
        atomnames += " " + std::to_string(bead->getResnr()) + ":" +
                     bead->getName() + ":" + std::to_string(bead->getId());
      }
      fragment.add("name", "UPTOYOU_BUTUNIQUE");
      fragment.add("mdatoms", atomnames);
      fragment.add("qmatoms", "IDS of QMATOMS i.e 0:C 1:H 2:C");
      fragment.add("mpoles", "IDS of MPOLES i.e 0:C 1:H 2:C");
      fragment.add("weights",
                   "weights for mapping(often atomic mass) i.e. 12  1 12");
      fragment.add("localframe", "IDs of up to 3 qmatoms or mpoles i.e. 0 1 2");
      std::ofstream template_mapfile(mapfile);
      template_mapfile << mapfile_prop << std::flush;
      template_mapfile.close();

      std::cout << "MOLECULETYPE " << csgmol->getName() << std::endl;
      std::cout << "SAMPLECOORDINATES" << std::endl;
      std::cout << "ID NAME COORDINATES[Angstroem] " << std::endl;
      for (const CSG::Bead* bead : sortedbeads) {
        Eigen::Vector3d pos = bead->getPos() * votca::tools::conv::nm2ang;
        std::string output =
            (boost::format("%1$i %2$s %3$+1.4f %4$+1.4f %5$+1.4f\n") %
             bead->getId() % bead->getName() % pos[0] % pos[1] % pos[2])
                .str();
        std::cout << output;
      }
    }
    std::cout << std::flush;
    return;
  }

  if (!TOOLS::filesystem::FileExists(mapfile)) {
    cout << endl
         << "xtp_map : map file '" << mapfile << "' could not be found."
         << endl;
    return;
  }
  XTP::Md2QmEngine md2qm(mapfile);

  votca::Index firstFrame = OptionsMap()["first-frame"].as<votca::Index>();
  votca::Index nFrames = OptionsMap()["nframes"].as<votca::Index>();
  bool beginAt = false;
  double time = OptionsMap()["begin"].as<double>();
  double startTime = mdtopol.getTime();
  if (time > 0.0) {
    beginAt = true;
    startTime = time;
  }

  // Extract first frame specified
  bool hasFrame;
  votca::Index frames_found = 0;
  votca::Index firstframecounter = firstFrame;
  for (hasFrame = true; hasFrame == true;
       hasFrame = trjread->NextFrame(mdtopol)) {
    frames_found++;
    if (((mdtopol.getTime() < startTime) && beginAt) || firstframecounter > 0) {
      firstframecounter--;
      continue;
    }
    break;
  }
  if (!hasFrame) {
    trjread->Close();

    throw runtime_error("Time or frame number exceeds trajectory length");
  }
  if (votca::Log::verbose()) {
    cout << "Read MD trajectory from " << trjfile << ": found " << frames_found
         << " frames, starting from frame " << firstFrame << endl;
  }
  // +++++++++++++++++++++++++ //
  // Convert MD to QM Topology //
  // +++++++++++++++++++++++++ //

  string statefile = OptionsMap()["file"].as<string>();
  if (TOOLS::filesystem::FileExists(statefile)) {
    cout << endl
         << "xtp_map : state file '" << statefile
         << "' already in use. Delete the current statefile or specify a "
            "different name."
         << endl;
    return;
  }

  XTP::StateSaver statsav(statefile);
  votca::Index laststep =
      -1;  // for some formats no step is given out so we check if the step
  for (votca::Index saved = 0; hasFrame && saved < nFrames;
       hasFrame = trjread->NextFrame(mdtopol), saved++) {
    if (mdtopol.getStep() == laststep) {
      mdtopol.setStep(laststep + 1);
    }
    laststep = mdtopol.getStep();
    XTP::Topology qmtopol = md2qm.map(mdtopol);
    statsav.WriteFrame(qmtopol);
  }
}

void XtpMap::ShowHelpText(std::ostream& out) {
  string name = ProgramName();
  if (VersionString() != "") {
    name = name + ", version " + VersionString();
  }
  XTP::HelpTextHeader(name);
  HelpText(out);
  out << "\n\n" << VisibleOptions() << endl;
}

int main(int argc, char** argv) {
  XtpMap xtpmap;
  return xtpmap.Exec(argc, argv);
}
