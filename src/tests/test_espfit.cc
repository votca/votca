/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE espfit_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/espfit.h>
#include <votca/xtp/logger.h>

#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(espfit_test)

BOOST_AUTO_TEST_CASE(esp_charges) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"H\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"5.447178e+00\">" << endl;
  basisfile << "        <contractions factor=\"1.562850e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"8.245470e-01\">" << endl;
  basisfile << "        <contractions factor=\"9.046910e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.831920e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  orbitals.setDFTbasisName("3-21G.xml");
  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(5);

  Eigen::MatrixXd MOs = Eigen::MatrixXd::Zero(17, 17);
  MOs << 0.9859859225723715, -0.2120110387683121, 6.027624018483543E-10,
      -6.033794875416886E-10, 8.465530119117907E-11, -0.18523741323314466,
      3.527736165571619E-8, 1.037690086043051E-7, 5.255631038084852E-8,
      -1.0944961312553891E-10, -2.9968082771885056E-9, 6.517433004000378E-9,
      -0.12272340526993321, 2.5237950468317497E-8, -2.7690325768547364E-8,
      -9.537078236009389E-9, 0.01606292231800701, 0.10844483908732148,
      0.20199082507264948, 8.78232426985557E-10, -2.3585285979179067E-9,
      4.437708312076492E-9, 0.08346521435359845, -2.1218394715753976E-8,
      -8.447017801137761E-8, -4.76075692528739E-8, 3.5462765433826436E-9,
      2.5512503603366535E-8, -4.875181879833227E-8, 0.15999743825256385,
      -5.9456554468735114E-8, -5.8929148376675974E-8, 2.7792950690956708E-8,
      1.9963860644471587, 2.7984965565630298E-12, 1.6663617397340769E-9,
      0.28960842945281356, -0.12094737412196328, -0.2180463931995849,
      -8.006235552049359E-8, 0.23893645403592437, -0.14723231835289058,
      -0.1659942697418058, 0.31420115849340813, -0.4946468228163605,
      -0.5346808701707676, -4.212370478680589E-9, 0.5810636384110652,
      0.5158725310043325, -0.1701943326462476, 2.3690303256800467E-8,
      -8.118951902984378E-12, -1.1166528494836455E-9, -0.012908988356424136,
      -0.34101497628399874, 0.1720107068369689, -8.23557032081976E-8,
      -0.18602104453276286, 4.40456679313588E-5, -0.2678026340056062,
      0.6317918131605849, -0.10471336978523607, 0.4681406931693265,
      -1.3836928630150902E-7, -0.46113470132617324, 0.33669466077451127,
      -0.5538203777308582, 2.1067861115890944E-8, 1.025289351265045E-13,
      1.280958226338368E-9, -0.24900974741963455, -0.12298803993368872,
      -0.26251405095846053, 1.2360343702195112E-7, 0.12094466479006177,
      0.2909376780082672, -0.08396274712436566, 0.3624889003114479,
      0.6112618917841202, -0.352480031589774, 2.3260350032740282E-7,
      0.28713297598075155, -0.5032273526491667, -0.5450159777411522,
      -2.7576203215902084E-8, -0.08004075171448947, 0.626198040671726,
      -8.379601154033952E-9, 1.4866114094229967E-8, -2.0533499101327408E-8,
      2.4366759747742277, -4.472660729002343E-7, -1.2816020267729011E-6,
      -6.451306005576548E-7, -3.1148991171487757E-9, -1.1707876154394198E-8,
      5.190575831483262E-9, 0.14241427494562903, 1.8706472713769244E-8,
      1.6040202728739854E-7, -1.6861047266992964E-8, -3.927614713439364,
      -1.30937555026551E-11, 6.182375824284834E-10, 0.2357672052355552,
      -0.09846199439133324, -0.17750929212528355, -3.467507759150073E-7,
      0.9840551841492139, -0.6063733795767939, -0.6836441171383165,
      -0.5262603717174537, 0.8284915752163928, 0.8955452067833676,
      2.4126725749333145E-9, -0.7191255578424186, -0.638444944952469,
      0.2106328663760439, -4.80625783310161E-8, 2.901707163371214E-11,
      -2.719796974935695E-10, -0.010509073496188591, -0.27761673874171866,
      0.140032121374806, -4.933391684717672E-7, -0.7661241260882522,
      1.8138513129375454E-4, -1.1029398700931667, -1.0581978284977842,
      0.1753860286546577, -0.7840960270125937, 1.6072582247045665E-7,
      0.5707012525134849, -0.4166940266506984, 0.6854092823991711,
      -5.1655410468833874E-8, -1.2734160613584285E-11, -1.1299301498433548E-9,
      -0.2027162433873264, -0.10012328223717212, -0.21370995188562458,
      7.109427235509553E-7, 0.49810827240258226, 1.1982210997982332,
      -0.34579889210230996, -0.6071382432137781, -1.0238119589683745,
      0.5903742090269714, -2.6402693882536973E-7, -0.35535636188530373,
      0.6227952360496355, 0.6745129389585675, 6.67333967733626E-8,
      -8.919247986643991E-4, 0.12117898570381394, 0.011449767264516308,
      -0.24187865843875817, -0.1275861971359611, -0.018943332480253657,
      -0.03960381098454161, -0.03274486477091314, 0.11794119072117701,
      0.5721297515087533, 0.005203991830835487, -0.18321530943426048,
      -0.6669015194699461, -0.28482656314471205, -0.2444379586548266,
      0.8879577896757126, 0.1497578658344628, 0.015461061593680181,
      0.018817629250991876, 0.010559070324500083, -0.22306240933331595,
      -0.11766098930445101, -0.9233263670640812, -0.527125423742969,
      -0.43583266941959903, 1.569793327053212, 0.08389194927879248,
      7.630601953325533E-4, -0.02686503624957739, 0.3024500079579037,
      0.4444391365794962, 0.3814170455764181, -1.385555816698486,
      0.7783794082761369, -8.919248045470733E-4, 0.12117898526589616,
      0.2280580591302763, 0.1418544881277567, -0.05273941403141442,
      -0.018943312072462044, -0.06925149615804158, 0.09982139594444281,
      -0.042317052904871144, -0.29736272105570133, -0.43776964616602054,
      -0.28435985208662384, -0.6669013868736686, -0.528329281418086,
      -0.4774882899758197, -0.6497833559357065, 0.1497578680753049,
      0.01546106162661308, 0.018817628977754658, 0.21031693955350184,
      0.1308193190261492, -0.048636699597722884, -0.9233260583427542,
      -0.9217350564778704, 1.328619612816386, -0.5632383147818923,
      -0.04360259051493696, -0.06419059744422875, -0.041695964112123236,
      0.30244981905593243, 0.8243970270915912, 0.7450652139549778,
      1.0139119766431184, 0.7783794196473948, -8.919248004546498E-4,
      0.12117898487524767, -0.2173822522234407, 0.14016684568632184,
      -0.08951442158699169, -0.01894339086438439, -0.015496472210454727,
      -0.0998012975505603, -0.07968920412760487, -0.25513539246973976,
      0.5293408395926812, -0.12502631071778403, -0.6669017240363495,
      -0.11699479694446009, 0.9486675033734427, -0.12524787912673363,
      0.14975792844168354, 0.015461061616917033, 0.0188176310519791,
      -0.20047162218236553, 0.12926296669473397, -0.08255089767776494,
      -0.9233276076225175, -0.20625745937855017, -1.328351560654868,
      -1.0606603214033368, -0.037410756977543196, 0.07761777552480653,
      -0.01833273029903981, 0.3024503076758718, 0.18255701220434814,
      -1.4802864794669026, 0.19543482055132402, 0.7783792539657504,
      -8.919247912428705E-4, 0.1211789852140613, -0.02212557455741454,
      -0.04014267667681461, 0.26984003637600523, -0.018943263725554098,
      0.12435181014743463, 0.032724829534187684, 0.004065085775822225,
      -0.019631635569613906, -0.0967751917821364, 0.5926014707241773,
      -0.6669011958339979, 0.9301511743837797, -0.2267419818007438,
      -0.11292675586543484, 0.14975791017932275, 0.015461061591032913,
      0.018817623087706944, -0.020404377511759473, -0.037019896520445024,
      0.24884861803588862, -0.923325904897998, 1.6551186156196893,
      0.4355665955914489, 0.05410631833795338, -0.0028786036975301174,
      -0.014190247201704715, 0.0868937812871138, 0.30244953232979943,
      -1.4513934569079558, 0.35380449165354727, 0.1762091332478184,
      0.7783793751383817;

  orbitals.MOCoefficients() = MOs;
  QMState gs = QMState("n");
  Logger log;
  Espfit esp = Espfit(log);
  esp.setUseSVD(1e-8);
  esp.Fit2Density(orbitals, gs, "medium");
  Eigen::VectorXd pcharges = Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  int index = 0;
  for (const auto& site : orbitals.Multipoles()) {
    pcharges(index) = site.getCharge();
    index++;
  }
  Eigen::VectorXd p_ref = Eigen::VectorXd::Zero(5);
  p_ref << -0.513812, 0.128201, 0.128537, 0.128537, 0.128537;

  bool check_esp_num = p_ref.isApprox(pcharges, 0.01);
  if (!check_esp_num) {
    cout << "ref" << endl;
    cout << p_ref << endl;
    cout << "calc" << endl;
    cout << pcharges << endl;
  }
  BOOST_CHECK_EQUAL(check_esp_num, 1);

  std::vector<std::pair<int, int> > pairconstraint;
  std::pair<int, int> p1;
  p1.first = 1;
  p1.second = 2;
  pairconstraint.push_back(p1);
  std::pair<int, int> p2;
  p2.first = 3;
  p2.second = 4;
  pairconstraint.push_back(p2);
  Espfit esp2 = Espfit(log);
  esp2.setUseSVD(1e-8);
  esp2.setPairConstraint(pairconstraint);
  esp2.Fit2Density(orbitals, gs, "medium");
  Eigen::VectorXd pcharges_equal =
      Eigen::VectorXd::Zero(orbitals.Multipoles().size());
  index = 0;
  for (const auto& site : orbitals.Multipoles()) {
    pcharges_equal(index) = site.getCharge();
    index++;
  }

  bool check_p1 = (std::abs(pcharges_equal(1) - pcharges_equal(2)) < 1e-6);
  bool check_p2 = (std::abs(pcharges_equal(3) - pcharges_equal(4)) < 1e-6);
  BOOST_CHECK_EQUAL(check_p1 && check_p2, 1);

  std::vector<QMFragment<double> > regionconstraint;

  std::string indeces = "1...3";
  QMFragment<double> reg = QMFragment<double>("constraint", 0, indeces);
  reg.value() = 1.0;
  regionconstraint.push_back(reg);
  Espfit esp3 = Espfit(log);
  esp3.setRegionConstraint(regionconstraint);
  esp3.setUseSVD(1e-8);
  esp3.Fit2Density(orbitals, gs, "medium");
  Eigen::VectorXd pcharges_reg =
      Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  index = 0;

  for (const auto& site : orbitals.Multipoles()) {
    pcharges_reg(index) = site.getCharge();
    index++;
  }

  bool check_reg = (std::abs(pcharges_reg.segment(1, 3).sum() - 1.0) < 1e-6);
  if (!check_reg) {
    std::cout << "All charges " << pcharges_reg << std::endl;
    std::cout << "Sum of charges 1,2,3 should equal 1:"
              << pcharges_reg.segment(1, 3).sum() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_reg, 1);
}

BOOST_AUTO_TEST_CASE(analytic_vs_numeric) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 1" << endl;
  xyzfile << " carbon" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile.close();

  ofstream basisfile("3-21G.xml");
  basisfile << "<basis name=\"3-21G\">" << endl;
  basisfile << "  <element name=\"C\">" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"S\">" << endl;
  basisfile << "      <constant decay=\"1.722560e+02\">" << endl;
  basisfile << "        <contractions factor=\"6.176690e-02\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"2.591090e+01\">" << endl;
  basisfile << "        <contractions factor=\"3.587940e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"5.533350e+00\">" << endl;
  basisfile << "        <contractions factor=\"7.007130e-01\" type=\"S\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"3.664980e+00\">" << endl;
  basisfile << "        <contractions factor=\"-3.958970e-01\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"2.364600e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "      <constant decay=\"7.705450e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.215840e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"8.606190e-01\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "    <shell scale=\"1.0\" type=\"SP\">" << endl;
  basisfile << "      <constant decay=\"1.958570e-01\">" << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"S\"/>"
            << endl;
  basisfile << "        <contractions factor=\"1.000000e+00\" type=\"P\"/>"
            << endl;
  basisfile << "      </constant>" << endl;
  basisfile << "    </shell>" << endl;
  basisfile << "  </element>" << endl;
  basisfile << "</basis>" << endl;
  basisfile.close();

  Orbitals orbitals;
  orbitals.setDFTbasisName("3-21G.xml");
  orbitals.QMAtoms().LoadFromFile("molecule.xyz");
  QMState gs = QMState("n");
  orbitals.MOCoefficients() = Eigen::MatrixXd::Identity(9, 9);
  Logger log;

  Espfit esp = Espfit(log);
  esp.setUseSVD(1e-8);

  esp.Fit2Density_analytic(orbitals, gs);
  Eigen::VectorXd pcharges_anal =
      Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  int index = 0;
  for (const StaticSite& atom : orbitals.Multipoles()) {
    pcharges_anal(index) = atom.getCharge();
    index++;
  }

  esp.Fit2Density(orbitals, gs, "medium");
  Eigen::VectorXd pcharges = Eigen::VectorXd::Zero(orbitals.QMAtoms().size());
  index = 0;
  for (const StaticSite& atom : orbitals.Multipoles()) {
    pcharges(index) = atom.getCharge();
    index++;
  }

  bool check_esp_ana = pcharges.isApprox(pcharges_anal, 0.01);
  if (!check_esp_ana) {
    cout << "numeric" << endl;
    cout << pcharges << endl;
    cout << "analytic" << endl;
    cout << pcharges_anal << endl;
  }
  BOOST_CHECK_EQUAL(check_esp_ana, 1);
}

BOOST_AUTO_TEST_SUITE_END()
