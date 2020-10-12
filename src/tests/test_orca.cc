/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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
#include <libint2/initialize.h>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE orca_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmpackagefactory.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(orca_test)

BOOST_AUTO_TEST_CASE(polar_test) {

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> orca =
      std::unique_ptr<QMPackage>(QMPackages().Create("orca"));
  Logger log;
  orca->setLog(&log);
  orca->setRunDir(std::string(XTP_TEST_DATA_FOLDER) + "/orca");
  orca->setLogFileName("polar_orca.log");
  Eigen::Matrix3d polar_mat = orca->GetPolarizability();

  Eigen::Matrix3d polar_ref = Eigen::Matrix3d::Zero();
  polar_ref << 11.40196, -0.00423, 0.00097, -0.00423, 11.42894, 0.01163,
      0.00097, 0.01163, 11.41930;
  bool polar_check = polar_ref.isApprox(polar_mat, 1e-5);
  if (!polar_check) {
    std::cout << "res" << std::endl;
    std::cout << polar_mat << std::endl;
    std::cout << "ref " << std::endl;
    std::cout << polar_ref << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(ext_charges_test) {
  libint2::initialize();
  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> orca =
      std::unique_ptr<QMPackage>(QMPackages().Create("orca"));
  Logger log;
  orca->setLog(&log);
  orca->setRunDir(std::string(XTP_TEST_DATA_FOLDER) + "/orca");
  orca->setLogFileName("orca_ext_charges.log");
  Orbitals orb;
  orca->ParseLogFile(orb);
  BOOST_CHECK_CLOSE(orb.getScaHFX(), 0.25, 1e-5);
  double ref_tot = -40.244504856356;  // HF - Self energy ext charges
  BOOST_CHECK_CLOSE(ref_tot, orb.getDFTTotalEnergy(), 1e-5);
  const QMMolecule& seg = orb.QMAtoms();
  double ang2bohr = votca::tools::conv::ang2bohr;
  QMMolecule ref("ref", 0);
  Eigen::Vector3d pos1 = {0.528800, 0.161000, 0.935900};
  QMAtom s1(0, "H", pos1 * ang2bohr);
  pos1 = {0.000000, 0.000000, 0.000000};
  QMAtom s2(1, "C", pos1 * ang2bohr);
  pos1 = {0.205100, 0.824000, -0.678600};
  QMAtom s3(2, "H", pos1 * ang2bohr);
  pos1 = {0.334500, -0.931400, -0.449600};
  QMAtom s4(3, "H", pos1 * ang2bohr);
  pos1 = {-1.068500, -0.053700, 0.192100};
  QMAtom s5(4, "H", pos1 * ang2bohr);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);
  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (votca::Index i = 0; i < votca::Index(seg.size()); i++) {
    BOOST_CHECK_EQUAL(ref[i].getPos().isApprox(seg[i].getPos(), 1e-5), true);
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }

  orb.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                      "/orca/3-21G_small.xml");
  orca->setMOsFileName("orca_ext_mos.gbw");

  orca->ParseMOsFile(orb);
  Eigen::VectorXd MOs_energy_ref = Eigen::VectorXd::Zero(17);
  MOs_energy_ref << -10.1441, -0.712523, -0.405481, -0.403251, -0.392086,
      0.157504, 0.196706, 0.202667, 0.240538, 0.718397, 0.727724, 0.728746,
      1.0577, 1.08125, 1.08607, 1.11627, 1.71596;
  bool check_eng = MOs_energy_ref.isApprox(orb.MOs().eigenvalues(), 1e-5);
  BOOST_CHECK_EQUAL(check_eng, true);
  if (!check_eng) {
    std::cout << "result eng" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
    std::cout << "ref eng" << std::endl;
    std::cout << MOs_energy_ref << std::endl;
  }
  Eigen::MatrixXd MOs_coeff_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/orca/MOs_coeff_ref.mm");

  // clang-format off
    std::array<Index, 49> votcaOrder_old = {
        0,                             // s
        0, -1, 1,                      // p
        0, -1, 1, -2, 2,               // d
        0, -1, 1, -2, 2, -3, 3,        // f
        0, -1, 1, -2, 2, -3, 3, -4, 4,  // g
        0, -1, 1, -2, 2, -3, 3, -4, 4,-5,5,  // h
        0, -1, 1, -2, 2, -3, 3, -4, 4,-5,5,-6,6  // i
    };
  // clang-format on

  std::array<votca::Index, 49> multiplier;
  multiplier.fill(1);
  OrbReorder ord(votcaOrder_old, multiplier);
  AOBasis aobasis = orb.SetupDftBasis();
  ord.reorderOrbitals(MOs_coeff_ref, aobasis);
  bool check_coeff = MOs_coeff_ref.isApprox(orb.MOs().eigenvectors(), 1e-5);
  BOOST_CHECK_EQUAL(check_coeff, true);
  if (!check_coeff) {
    std::cout << "result coeff" << std::endl;
    std::cout << orb.MOs().eigenvectors() << std::endl;
    std::cout << "ref coeff" << std::endl;
    std::cout << MOs_coeff_ref << std::endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(charges_test) {
  libint2::initialize();
  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> orca =
      std::unique_ptr<QMPackage>(QMPackages().Create("orca"));
  Logger log;
  orca->setLog(&log);
  orca->setRunDir(std::string(XTP_TEST_DATA_FOLDER) + "/orca");
  orca->setLogFileName("orca_charges.log");
  StaticSegment seg = orca->GetCharges();

  double ang2bohr = votca::tools::conv::ang2bohr;
  StaticSegment ref("ref", 0);
  Eigen::Vector3d pos1 = {0.528800, 0.161000, 0.935900};
  StaticSite s1(0, "H", pos1 * ang2bohr);
  s1.setCharge(0.178516);
  pos1 = {0.000000, 0.000000, 0.000000};
  StaticSite s2(1, "C", pos1 * ang2bohr);
  s2.setCharge(-0.709443);
  pos1 = {0.205100, 0.824000, -0.678600};
  StaticSite s3(2, "H", pos1 * ang2bohr);
  s3.setCharge(0.176227);
  pos1 = {0.334500, -0.931400, -0.449600};
  StaticSite s4(3, "H", pos1 * ang2bohr);
  s4.setCharge(0.177109);
  pos1 = {-1.068500, -0.053700, 0.192100};
  StaticSite s5(4, "H", pos1 * ang2bohr);
  s5.setCharge(0.177591);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);

  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (votca::Index i = 0; i < votca::Index(seg.size()); i++) {
    BOOST_CHECK_EQUAL(ref[i].Q().isApprox(seg[i].Q(), 1e-5), true);
    BOOST_CHECK_EQUAL(ref[i].getPos().isApprox(seg[i].getPos(), 1e-5), true);
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(opt_test) {

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> orca =
      std::unique_ptr<QMPackage>(QMPackages().Create("orca"));
  Logger log;
  orca->setLog(&log);
  orca->setRunDir(std::string(XTP_TEST_DATA_FOLDER) + "/orca");
  orca->setLogFileName("orca_opt.log");
  Orbitals orb;
  orca->ParseLogFile(orb);

  double ref_tot = -40.240899574987;  // HF - Self energy ext charges
  BOOST_CHECK_CLOSE(ref_tot, orb.getDFTTotalEnergy(), 1e-5);
  const QMMolecule& seg = orb.QMAtoms();
  double ang2bohr = votca::tools::conv::ang2bohr;
  QMMolecule ref("ref", 0);
  Eigen::Vector3d pos1 = {0.457282, 0.159495, 0.949475};
  QMAtom s1(0, "H", pos1 * ang2bohr);
  pos1 = {-0.060043, 0.000032, 0.000064};
  QMAtom s2(1, "C", pos1 * ang2bohr);
  pos1 = {0.161203, 0.828181, -0.679453};
  QMAtom s3(2, "H", pos1 * ang2bohr);
  pos1 = {0.279541, -0.938419, -0.446827};
  QMAtom s4(3, "H", pos1 * ang2bohr);
  pos1 = {-1.138083, -0.049389, 0.176541};
  QMAtom s5(4, "H", pos1 * ang2bohr);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);
  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (votca::Index i = 0; i < votca::Index(seg.size()); i++) {
    BOOST_CHECK_EQUAL(ref[i].getPos().isApprox(seg[i].getPos(), 1e-5), true);
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
}

BOOST_AUTO_TEST_CASE(input_generation_version_4_0_1) {
  unsetenv("VOTCASHARE");
  std::ofstream defaults("user_input.xml");

  defaults << "<package>\n"
           << "<name>orca</name>\n"
           << "<charge>0</charge>\n"
           << "<spin>1</spin>\n"
           << "<executable>some/path/orca</executable>\n"
           << "<basisset>" << std::string(XTP_TEST_DATA_FOLDER)
           << "/orca/3-21G.xml</basisset>\n"
           << "<functional>pbe0</functional>\n"
           << "<read_guess>false</read_guess>\n"
           << "<write_charges>false</write_charges>\n"
           << "<scratch>/tmp/qmpackage</scratch>\n"
           << "<optimize>false</optimize>\n"
           << "<convergence_tightness>tight</convergence_tightness>\n"
           << "<orca>\n"
           << "<method></method>\n"
           << "<scf>GUESS PMODEL</scf>\n"
           << "<maxcore>3000</maxcore>\n"
           << "</orca>\n"
           << "</package>";
  defaults.close();

  votca::tools::Property prop;
  prop.LoadFromXML("user_input.xml");

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> orca =
      std::unique_ptr<QMPackage>(QMPackages().Create("orca"));
  Logger log;
  orca->setLog(&log);
  orca->setRunDir(".");
  orca->Initialize(prop);

  Orbitals orb;
  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/orca/co.xyz");
  orca->WriteInputFile(orb);

  std::string input_file = "system.inp";
  if (!tools::filesystem::FileExists(input_file)) {
    throw std::runtime_error("Orca input file does not exists!");
  }
  std::ifstream file_input(input_file);
  std::stringstream buffer;
  buffer << file_input.rdbuf();
  std::string inp = buffer.str();

  // check basis section
  auto index1 = inp.find("%basis");
  auto index2 = inp.find("end", index1);
  BOOST_CHECK_EQUAL(inp.substr(index1, index2 - index1),
                    "%basis\nGTOName =\"system.bas\";\n");

  // check scf section multiline
  index1 = inp.find("%scf");
  index2 = inp.find("end", index1);
  BOOST_CHECK_EQUAL(inp.substr(index1, index2 - index1),
                    "%scf\nGUESS PMODEL\n");

  // Check method
  index1 = inp.find("!");
  BOOST_CHECK_EQUAL(inp.substr(index1), "! DFT pbe0   \n");

  // Check singleline orca kewords
  index1 = inp.find("%maxcore");
  index2 = inp.find("\n", index1);
  std::cout << "\ninp: " << inp.substr(index1) << "\n";
  BOOST_CHECK_EQUAL(inp.substr(index1, index2 - index1), "%maxcore 3000");
}

BOOST_AUTO_TEST_SUITE_END()
