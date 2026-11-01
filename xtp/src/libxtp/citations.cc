/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#include <algorithm>
#include <iterator>
#include <map>
#include <mutex>
#include <sstream>

// Local VOTCA includes
#include "votca/xtp/citations.h"

namespace votca {
namespace xtp {
namespace citations {

namespace {

std::mutex& Lock() {
  static std::mutex lock;
  return lock;
}

std::vector<std::string>& UsedKeys() {
  static std::vector<std::string> used;
  return used;
}

const std::vector<Optional>& GWBSEOptions() {
  static const std::vector<Optional> options = {
      {"Wehner:2018", "GW-BSE"},
      {"Tirimbo:2020a", "GW-BSE"},
      {"vanSchilfgaarde:2006", "quasiparticle self-consistent GW (gw.do_qsgw)"},
      {"Marie:2023", "SRG regularization of QSGW"},
  };
  return options;
}

const std::vector<Optional>& DFTOptions() {
  static const std::vector<Optional> options = {
      {"Wu:2005", "constrained DFT (cdft)"},
      {"Helmich-Paris:2021", "second-order (TRAH) SCF steps"},
      {"Manby:2012", "projection-based DFT-in-DFT embedding"},
      {"Sundaram:2024", "projection-based DFT-in-DFT embedding"},
  };
  return options;
}

std::vector<Optional> Join(std::vector<Optional> a,
                           const std::vector<Optional>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

const std::map<std::string, CalculatorEntry>& Table() {
  const std::vector<std::string> structure = {"Ruehle:2011"};
  const std::vector<Optional> electronic = Join(DFTOptions(), GWBSEOptions());
  const std::vector<Optional> embedded =
      Join(electronic, {{"Baumeier:2014", "GW-BSE region"},
                        {"Tirimbo:2020b", "GW-BSE region"},
                        {"Tirimbo:2024a", "GW-BSE region"},
                        {"Thole:1981", "polar region"},
                        {"Poelking:2016", "Ewald region"}});
  const std::vector<Optional> couplings =
      Join(electronic, {{"Baumeier:2010", "DFT couplings (dftcoupling)"},
                        {"Wehner:2017", "BSE couplings (bsecoupling)"}});
  const std::vector<Optional> grw = {
      {"Chen:2025a", "graph random walk (solver exact or grw)"}};
  const std::vector<Optional> network = {
      {"Chen:2025b", "traps (trap finder or trap templates)"}};

  static const std::map<std::string, CalculatorEntry> table = {
      // calculators (xtp_run)
      {"neighborlist", {structure, {}}},
      {"mapchecker", {structure, {}}},
      {"ianalyze", {structure, {}}},
      {"eanalyze", {structure, {}}},
      {"einternal", {structure, {}}},
      {"vaverage", {structure, {}}},
      {"kmcmultiple", {structure, grw}},
      {"kmctof", {structure, grw}},
      {"kmclifetime", {{"Ruehle:2011", "Wehner:2020"}, {}}},
      {"grw", {{"Chen:2025a", "Stenzel:2014"}, {}}},
      {"ewaldbackground", {{"Poelking:2016", "Poelking:2015"}, {}}},
      // job calculators (xtp_parallel)
      {"iqm", {structure, couplings}},
      {"ipodcoupling", {{"Ghan:2020", "Baumeier:2010"}, DFTOptions()}},
      {"eqm", {{}, electronic}},
      {"iexcitoncl", {{"Thole:1981"}, {}}},
      {"qmmm", {structure, embedded}},
      // tools (xtp_tools)
      {"dftgwbse", {{}, electronic}},
      {"spectrum", {{}, GWBSEOptions()}},
      {"excitoncoupling", {{"Wehner:2017"}, electronic}},
      {"coupling", {{"Baumeier:2010"}, {}}},
      {"podcoupling", {{"Ghan:2020", "Baumeier:2010"}, {}}},
      {"molpol", {{"Thole:1981"}, {}}},
      {"partialcharges", {{"Breneman:1990"}, {}}},
      {"diabatization", {{"Tirimbo:2024b"}, {}}},
      {"latticegen", {{"Chen:2025a"}, {}}},
      {"networkfit", {{"Baumeier:2012", "Kordt:2014"}, network}},
      {"networkgen", {{"Baumeier:2012", "Kordt:2014"}, network}},
  };
  return table;
}

}  // namespace

const std::vector<Reference>& References() {
  static const std::vector<Reference> references = {
      // the VOTCA paper
      {"Baumeier:2024",
       "Baumeier, Wehner, Renaud, Zapata Ruiz, Halver, Madhikar, Gerritsen, "
       "Tirimbo, Sijen, Rosenberger, Brown, Sundaram, Krajniak, Bernhardt, and "
       "Junghans, VOTCA: multiscale frameworks for quantum and classical "
       "simulations in soft matter, J. Open Source Softw. 9, 6864 (2024)",
       "10.21105/joss.06864"},
      // methods developed with VOTCA-XTP
      {"Ruehle:2011",
       "Ruehle, Lukyanov, May, Schrader, Vehoff, Kirkpatrick, Baumeier, and "
       "Andrienko, Microscopic simulations of charge transport in disordered "
       "organic semiconductors, J. Chem. Theory Comput. 7, 3335 (2011)",
       "10.1021/ct200388s"},
      {"Baumeier:2010",
       "Baumeier, Kirkpatrick, and Andrienko, Density-functional based "
       "determination of intermolecular charge transfer properties for "
       "large-scale morphologies, Phys. Chem. Chem. Phys. 12, 11103 (2010)",
       "10.1039/C002337J"},
      {"Baumeier:2012",
       "Baumeier, Stenzel, Poelking, Andrienko, and Schmidt, Stochastic "
       "modeling of molecular charge transport networks, Phys. Rev. B 86, "
       "184202 (2012)",
       "10.1103/PhysRevB.86.184202"},
      {"Baumeier:2014",
       "Baumeier, Rohlfing, and Andrienko, Electronic excitations in push-pull "
       "oligomers and their complexes with fullerene from many-body Green's "
       "functions theory with polarizable embedding, J. Chem. Theory Comput. "
       "10, 3104 (2014)",
       "10.1021/ct500479f"},
      {"Kordt:2014",
       "Kordt, Stenzel, Baumeier, Schmidt, and Andrienko, Parametrization of "
       "extended Gaussian disorder models from microscopic charge transport "
       "simulations, J. Chem. Theory Comput. 10, 2508 (2014)",
       "10.1021/ct500269r"},
      {"Stenzel:2014",
       "Stenzel et al., A general framework for consistent estimation of "
       "charge transport properties via random walks in random environments, "
       "Multiscale Model. Simul. 12, 1108 (2014)",
       "10.1137/130942504"},
      {"Poelking:2015",
       "Poelking, Tietze, Elschner, Olthof, Hertel, Baumeier, Wuerthner, "
       "Meerholz, Leo, and Andrienko, Impact of mesoscale order on "
       "open-circuit voltage in organic solar cells, Nat. Mater. 14, 434 "
       "(2015)",
       "10.1038/nmat4167"},
      {"Poelking:2016",
       "Poelking and Andrienko, Long-range embedding of molecular ions and "
       "excitations in a polarizable molecular environment, J. Chem. Theory "
       "Comput. 12, 4516 (2016)",
       "10.1021/acs.jctc.6b00599"},
      {"Wehner:2017",
       "Wehner and Baumeier, Intermolecular singlet and triplet exciton "
       "transfer integrals from many-body Green's functions theory, J. Chem. "
       "Theory Comput. 13, 1584 (2017)",
       "10.1021/acs.jctc.6b00935"},
      {"Wehner:2018",
       "Wehner, Brombacher, Brown, Junghans, Caylak, Khalak, Madhikar, "
       "Tirimbo, and Baumeier, Electronic excitations in complex molecular "
       "environments: many-body Green's functions theory in VOTCA-XTP, J. "
       "Chem. Theory Comput. 14, 6253 (2018)",
       "10.1021/acs.jctc.8b00617"},
      {"Wehner:2020",
       "Wehner and Baumeier, Multiscale simulations of singlet and triplet "
       "exciton dynamics in energetically disordered molecular systems based "
       "on many-body Green's functions theory, New J. Phys. 22, 033033 (2020)",
       "10.1088/1367-2630/ab7a04"},
      {"Tirimbo:2020a",
       "Tirimbo, Sundaram, Caylak, Scharpach, Sijen, Junghans, Brown, Zapata "
       "Ruiz, Renaud, Wehner, and Baumeier, Excited-state electronic structure "
       "of molecules using many-body Green's functions: quasiparticles and "
       "electron-hole excitations with VOTCA-XTP, J. Chem. Phys. 152, 114103 "
       "(2020)",
       "10.1063/1.5144277"},
      {"Tirimbo:2020b",
       "Tirimbo, de Vries, Weijtens, Bobbert, Neumann, Coehoorn, and Baumeier, "
       "Quantitative predictions of photoelectron spectra in amorphous "
       "molecular solids from multiscale quasiparticle embedding, Phys. Rev. B "
       "101, 035402 (2020)",
       "10.1103/PhysRevB.101.035402"},
      {"Tirimbo:2024a",
       "Tirimbo, Sundaram, and Baumeier, Embedded many-body Green's function "
       "methods for complex molecular systems, WIREs Comput. Mol. Sci. 14, "
       "e1734 (2024)",
       "10.1002/wcms.1734"},
      {"Tirimbo:2024b",
       "Tirimbo and Baumeier, Electronic couplings and conversion dynamics "
       "between localized and charge transfer excitations from many-body "
       "Green's functions theory, J. Chem. Theory Comput. 20, 4605 (2024)",
       "10.1021/acs.jctc.4c00142"},
      {"Sundaram:2024",
       "Sundaram and Baumeier, Quantum-quantum and quantum-quantum-classical "
       "schemes for near-gap excitations with projection-based-embedded "
       "GW-Bethe-Salpeter equation, J. Chem. Theory Comput. 20, 5451 (2024)",
       "10.1021/acs.jctc.4c00163"},
      {"Chen:2025a",
       "Chen, van der Hoorn, and Baumeier, A graph random walk method for "
       "calculating time-of-flight charge mobility in organic semiconductors "
       "from multiscale simulations, Multiscale Model. Simul. 23, 1439 (2025)",
       "10.1137/24M1658085"},
      {"Chen:2025b",
       "Chen, van der Hoorn, and Baumeier, Automatic identification of traps "
       "in molecular charge transport networks of organic semiconductors, "
       "Phys. Rev. B 111, 224209 (2025)",
       "10.1103/shyl-kcfp"},
      // methods by others
      {"Thole:1981",
       "Thole, Molecular polarizabilities calculated with a modified dipole "
       "interaction, Chem. Phys. 59, 341 (1981)",
       "10.1016/0301-0104(81)85176-2"},
      {"Breneman:1990",
       "Breneman and Wiberg, Determining atom-centered monopoles from "
       "molecular electrostatic potentials, J. Comput. Chem. 11, 361 (1990)",
       "10.1002/jcc.540110311"},
      {"Wu:2005",
       "Wu and Van Voorhis, Direct optimization method to study constrained "
       "systems within density-functional theory, Phys. Rev. A 72, 024502 "
       "(2005)",
       "10.1103/PhysRevA.72.024502"},
      {"vanSchilfgaarde:2006",
       "van Schilfgaarde, Kotani, and Faleev, Quasiparticle self-consistent GW "
       "theory, Phys. Rev. Lett. 96, 226402 (2006)",
       "10.1103/PhysRevLett.96.226402"},
      {"Manby:2012",
       "Manby, Stella, Goodpaster, and Miller, A simple, exact "
       "density-functional-theory embedding scheme, J. Chem. Theory Comput. 8, "
       "2564 (2012)",
       "10.1021/ct300544e"},
      {"Ghan:2020",
       "Ghan, Kunkel, Reuter, and Oberhofer, Improved projection-operator "
       "diabatization schemes for the calculation of electronic coupling "
       "values, J. Chem. Theory Comput. 16, 7431 (2020)",
       "10.1021/acs.jctc.0c00887"},
      {"Helmich-Paris:2021",
       "Helmich-Paris, A trust-region augmented Hessian implementation for "
       "restricted and unrestricted Hartree-Fock and Kohn-Sham methods, J. "
       "Chem. Phys. 154, 164104 (2021)",
       "10.1063/5.0040798"},
      {"Marie:2023",
       "Marie and Loos, A similarity renormalization group approach to Green's "
       "function methods, J. Chem. Theory Comput. 19, 3943 (2023)",
       "10.1021/acs.jctc.3c00281"},
  };
  return references;
}

bool Exists(const std::string& key) {
  const std::vector<Reference>& refs = References();
  return std::any_of(refs.begin(), refs.end(),
                     [&key](const Reference& r) { return r.key == key; });
}

std::string Format(const std::string& key) {
  const std::vector<Reference>& refs = References();
  auto it = std::find_if(refs.begin(), refs.end(),
                         [&key](const Reference& r) { return r.key == key; });
  if (it == refs.end()) {
    return key;
  }
  return it->doi.empty() ? it->text : it->text + ", doi:" + it->doi;
}

void Cite(const std::string& key) {
  std::lock_guard<std::mutex> guard(Lock());
  std::vector<std::string>& used = UsedKeys();
  if (std::find(used.begin(), used.end(), key) == used.end()) {
    used.push_back(key);
  }
}

std::vector<std::string> Used() {
  std::lock_guard<std::mutex> guard(Lock());
  return UsedKeys();
}

void Reset() {
  std::lock_guard<std::mutex> guard(Lock());
  UsedKeys().clear();
}

namespace {
// Greedy word wrap to at most width characters per line. The first line
// starts with first, continuation lines with indent. Words (e.g. DOIs) are
// never split, so a single word longer than the line may exceed width.
std::string Wrap(const std::string& text, const std::string& first,
                 const std::string& indent, std::size_t width = 78) {
  std::istringstream words(text);
  std::string out = first;
  std::size_t line_length = first.size();
  bool line_empty = true;
  std::string word;
  while (words >> word) {
    if (!line_empty && line_length + 1 + word.size() > width) {
      out += "\n" + indent;
      line_length = indent.size();
      line_empty = true;
    }
    if (!line_empty) {
      out += " ";
      line_length++;
    }
    out += word;
    line_length += word.size();
    line_empty = false;
  }
  return out + "\n";
}
}  // namespace

std::string Report() {
  const std::vector<std::string> used = Used();
  if (used.empty()) {
    return "";
  }
  const std::string rule(78, '=');
  std::ostringstream out;
  out << rule << "\n"
      << "Please cite the methods used in this run:\n";
  std::size_t n = 1;
  for (const std::string& key : used) {
    const std::string label = "  [" + std::to_string(n++) + "] ";
    out << Wrap(Format(key), label, std::string(label.size(), ' '));
  }
  out << rule << "\n";
  return out.str();
}

CalculatorEntry ForCalculator(const std::string& name) {
  const auto& table = Table();
  auto it = table.find(name);
  return it == table.end() ? CalculatorEntry{} : it->second;
}

std::vector<std::string> Calculators() {
  const auto& table = Table();
  std::vector<std::string> names;
  names.reserve(table.size());
  std::transform(table.begin(), table.end(), std::back_inserter(names),
                 [](const auto& entry) { return entry.first; });
  return names;
}

std::string CalculatorReport(const std::string& name) {
  const CalculatorEntry entry = ForCalculator(name);
  std::ostringstream out;
  out << "References for " << name << ":\n";
  out << Wrap(Format("Baumeier:2024"), "  ", "    ");
  for (const std::string& key : entry.keys) {
    out << Wrap(Format(key), "  ", "    ");
  }
  if (!entry.optional.empty()) {
    out << "Depending on the options:\n";
    for (const Optional& opt : entry.optional) {
      out << Wrap(opt.when + ": " + Format(opt.key), "  ", "    ");
    }
  }
  return out.str();
}

}  // namespace citations
}  // namespace xtp
}  // namespace votca
