/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "votca/tools/application.h"
#include "votca/xtp/statesaver.h"
#include <votca/csg/topologyreader.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/tools/filesystem.h>
#include <votca/tools/globals.h>
#include <votca/xtp/md2qmengine.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/version.h>

using namespace std;

namespace CSG = votca::csg;
namespace XTP = votca::xtp;
namespace TOOLS = votca::tools;

class XtpMap : public TOOLS::Application {

 public:
  string ProgramName() { return "xtp_map"; }
  void HelpText(ostream& out) { out << "Generates QM|MD topology" << endl; }
  void ShowHelpText(std::ostream& out);

  void Initialize();
  bool EvaluateOptions();
  void Run();

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
}

bool XtpMap::EvaluateOptions() {

  CheckRequired("topology", "Missing topology file");
  CheckRequired("segments", "Missing segment definition file");
  CheckRequired("coordinates", "Missing trajectory input");
  if (!(_op_vm.count("makesegments"))) {
    CheckRequired("file", "Missing state file");
  }
  return 1;
}

void XtpMap::Run() {

  std::string name = ProgramName();
  if (VersionString() != "") name = name + ", version " + VersionString();
  XTP::HelpTextHeader(name);

  // ++++++++++++++++++++++++++++ //
  // Create MD topology from file //
  // ++++++++++++++++++++++++++++ //

  // Create topology reader
  string topfile = _op_vm["topology"].as<string>();
  std::unique_ptr<CSG::TopologyReader> topread =
      std::unique_ptr<CSG::TopologyReader>(
          CSG::TopReaderFactory().Create(topfile));

  if (topread == nullptr) {
    throw runtime_error(string("Input format not supported: ") +
                        _op_vm["topology"].as<string>());
  }
  CSG::Topology mdtopol;
  topread->ReadTopology(topfile, mdtopol);
  if (TOOLS::globals::verbose) {
    cout << "Read MD topology from " << topfile << ": Found "
         << mdtopol.BeadCount() << " atoms in " << mdtopol.MoleculeCount()
         << " molecules. " << endl;
  }

  string mapfile = _op_vm["segments"].as<string>();
  if (_op_vm.count("makesegments")) {
    if (TOOLS::filesystem::FileExists(mapfile)) {
      cout << endl
           << "xtp_map : map file '" << mapfile
           << "' already in use. Delete the current mapfile or specify a "
              "different name."
           << endl;
      return;
    }

    TOOLS::Property mapfile_prop("topology", "", "");
    TOOLS::Property& molecules = mapfile_prop.add("molecules", "");

    std::map<std::string, int> firstmoleculeid;

    std::map<std::string, int> molecule_names;
    for (const CSG::Molecule* mol : mdtopol.Molecules()) {
      if (!molecule_names.count(mol->getName())) {
        firstmoleculeid[mol->getName()] = mol->getId();
      }
      molecule_names[mol->getName()]++;
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
      TOOLS::Property& fragments = segment.add("fragments", "");
      TOOLS::Property& fragment = fragments.add("fragment", "");
      std::string atomnames = "";
      const CSG::Molecule* csgmol =
          mdtopol.getMolecule(firstmoleculeid[mol.first]);
      for (const CSG::Bead* bead : csgmol->Beads()) {
        atomnames +=
            "\t" + std::to_string(bead->getId()) + ":" + bead->getName();
      }
      TOOLS::Property& mdatoms = fragment.add("mdatoms", atomnames);
    }

    std::ofstream template_mapfile(mapfile);
    template_mapfile << mapfile_prop << std::flush;
    template_mapfile.close();
    return;
  }

  if (!TOOLS::filesystem::FileExists(mapfile)) {
    cout << endl
         << "xtp_map : map file '" << mapfile << "' could not be found."
         << endl;
    return;
  }
  XTP::Md2QmEngine md2qm(mapfile);

  // ++++++++++++++++++++++++++++++ //
  // Create MD trajectory from file //
  // ++++++++++++++++++++++++++++++ //

  // Create trajectory reader and initialize
  string trjfile = _op_vm["coordinates"].as<string>();
  std::unique_ptr<CSG::TrajectoryReader> trjread =
      std::unique_ptr<CSG::TrajectoryReader>(
          CSG::TrjReaderFactory().Create(trjfile));

  if (trjread == nullptr) {
    throw runtime_error(string("Input format not supported: ") +
                        _op_vm["coordinates"].as<string>());
  }
  trjread->Open(trjfile);
  trjread->FirstFrame(mdtopol);

  int firstFrame = 0;
  int nFrames = 1;
  bool beginAt = 0;
  double startTime = mdtopol.getTime();

  if (_op_vm.count("nframes")) {
    nFrames = _op_vm["nframes"].as<int>();
  }
  if (_op_vm.count("first-frame")) {
    firstFrame = _op_vm["first-frame"].as<int>();
  }
  if (_op_vm.count("begin")) {
    beginAt = true;
    startTime = _op_vm["begin"].as<double>();
  }

  // Extract first frame specified
  bool hasFrame;

  for (hasFrame = true; hasFrame == true;
       hasFrame = trjread->NextFrame(mdtopol)) {
    if (((mdtopol.getTime() < startTime) && beginAt) || firstFrame > 0) {
      firstFrame--;
      continue;
    }
    break;
  }
  if (!hasFrame) {
    trjread->Close();

    throw runtime_error("Time or frame number exceeds trajectory length");
  }
  if (TOOLS::globals::verbose) {
    cout << "Read MD trajectory from " << trjfile << ": found " << nFrames
         << " frames, starting from frame " << firstFrame << endl;
  }
  // +++++++++++++++++++++++++ //
  // Convert MD to QM Topology //
  // +++++++++++++++++++++++++ //
  int laststep =
      -1;  // for some formats no step is given out so we check if the step

  string statefile = _op_vm["file"].as<string>();
  if (TOOLS::filesystem::FileExists(statefile)) {
    cout << endl
         << "xtp_map : state file '" << statefile
         << "' already in use. Delete the current statefile or specify a "
            "different name."
         << endl;
    return;
  }

  XTP::StateSaver statsav(statefile);
  for (int saved = 0; hasFrame && saved < nFrames;
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
  if (VersionString() != "") name = name + ", version " + VersionString();
  XTP::HelpTextHeader(name);
  HelpText(out);
  out << "\n\n" << VisibleOptions() << endl;
}

int main(int argc, char** argv) {
  XtpMap xtpmap;
  return xtpmap.Exec(argc, argv);
}
