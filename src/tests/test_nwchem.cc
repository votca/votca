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

#define BOOST_TEST_MODULE nwchem_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmpackagefactory.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(nwchem_test)

BOOST_AUTO_TEST_CASE(polar_test) {
  std::ofstream polar("polar_nwchem.log");
  polar << " argument  1 = system.nw" << std::endl;
  polar << "                                         " << std::endl;
  polar << "                                         " << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar
      << "              Northwest Computational Chemistry Package (NWChem) 6.8"
      << std::endl;
  polar
      << "              ------------------------------------------------------"
      << std::endl;
  polar << "                             -------------------------"
        << std::endl;
  polar << "" << std::endl;
  polar << " Output coordinates in angstroms (scale by  1.889725989 to convert "
           "to a.u.)"
        << std::endl;
  polar << "" << std::endl;
  polar << "  No.       Tag          Charge          X              Y          "
           "    Z"
        << std::endl;
  polar << " ---- ---------------- ---------- -------------- -------------- "
           "--------------"
        << std::endl;
  polar << "    1 H                    1.0000     0.52881000     0.16101000    "
           " 0.93592000"
        << std::endl;
  polar << "    2 C                    6.0000     0.00001000     0.00001000    "
           " 0.00002000"
        << std::endl;
  polar << "    3 H                    1.0000     0.20511000     0.82401000    "
           "-0.67858000"
        << std::endl;
  polar << "    4 H                    1.0000     0.33451000    -0.93139000    "
           "-0.44958000"
        << std::endl;
  polar << "    5 H                    1.0000    -1.06849000    -0.05369000    "
           " 0.19212000"
        << std::endl;
  polar << "" << std::endl;
  polar << "      Atomic Mass " << std::endl;
  polar << "      ----------- " << std::endl;
  polar << "" << std::endl;
  polar << "      H                  1.007825" << std::endl;
  polar << "      C                 12.000000" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " Effective nuclear repulsion energy (a.u.)      13.4728439460"
        << std::endl;
  polar << "" << std::endl;
  polar << "            Nuclear Dipole moment (a.u.) " << std::endl;
  polar << "            ----------------------------" << std::endl;
  polar << "        X                 Y               Z" << std::endl;
  polar << " ---------------- ---------------- ----------------" << std::endl;
  polar << "     0.0000000000    -0.0000000000     0.0000000000" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "            XYZ format geometry" << std::endl;
  polar << "            -------------------" << std::endl;
  polar << "     5" << std::endl;
  polar << " geometry" << std::endl;
  polar << " H                     0.52881000     0.16101000     0.93592000"
        << std::endl;
  polar << " C                     0.00001000     0.00001000     0.00002000"
        << std::endl;
  polar << " H                     0.20511000     0.82401000    -0.67858000"
        << std::endl;
  polar << " H                     0.33451000    -0.93139000    -0.44958000"
        << std::endl;
  polar << " H                    -1.06849000    -0.05369000     0.19212000"
        << std::endl;
  polar << "" << std::endl;
  polar << "   convergence    iter        energy       DeltaE   RMS-Dens  "
           "Diis-err    time"
        << std::endl;
  polar << " ---------------- ----- ----------------- --------- --------- "
           "---------  ------"
        << std::endl;
  polar << " d= 0,ls=0.0,diis     1    -40.2408052356 -5.37D+01  9.16D-09  "
           "2.20D-14     0.2"
        << std::endl;
  polar << " d= 0,ls=0.0,diis     2    -40.2408052356 -3.55D-14  4.75D-09  "
           "1.88D-14     0.3"
        << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "         Total DFT energy =      -40.240805235609" << std::endl;
  polar << "      One electron energy =      -79.605213751416" << std::endl;
  polar << "           Coulomb energy =       32.718690448038" << std::endl;
  polar << "    Exchange-Corr. energy =       -6.827125878236" << std::endl;
  polar << " Nuclear repulsion energy =       13.472843946005" << std::endl;
  polar << "" << std::endl;
  polar << " Numeric. integr. density =       10.000000095313" << std::endl;
  polar << "" << std::endl;
  polar << "     Total iterative time =      0.3s" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " Parallel integral file used       1 records with       0 large "
           "values"
        << std::endl;
  polar << "" << std::endl;
  polar << " *** CALLING NEW AORESP DRIVER FOR CLOSED AND OPEN SHELLS ***"
        << std::endl;
  polar << " Entering AOResponse driver routine" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "-------------------------------------------------------------------"
           "-------------"
        << std::endl;
  polar << "" << std::endl;
  polar << "          ****************" << std::endl;
  polar << "          *** RESPONSE ***" << std::endl;
  polar << "          ****************" << std::endl;
  polar << "" << std::endl;
  polar << " Response module for NWChem and dynamic CPKS solver" << std::endl;
  polar << " developed by J. Autschbach and coworkers, SUNY Buffalo"
        << std::endl;
  polar << " The methodology used in this program is described in "
        << std::endl;
  polar << " ChemPhysChem 12 (2011), 3224-3235 (main reference)" << std::endl;
  polar << " J. Chem. Phys. 123 (2005), 114103" << std::endl;
  polar << " J. Chem. Phys. 122 (2005), 224115" << std::endl;
  polar << " J. Chem. Phys. 122 (2005), 074105" << std::endl;
  polar << " Comp. Lett. 3 (2007), 131-150 (contact JA for a copy)"
        << std::endl;
  polar << " Please cite this work in publications based on results"
        << std::endl;
  polar << " obtained with this code. Thank you!" << std::endl;
  polar << "" << std::endl;
  polar << " For extension of response module to open shell" << std::endl;
  polar << " calculations please acknowledge:" << std::endl;
  polar << " F. Aquino, Northwestern University, Schatz Rsrch Group"
        << std::endl;
  polar << " The update to the methodology is described in" << std::endl;
  polar << " J. Phys. Chem. A 118 (2014) 517-525" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "          -----------------------------------------------"
        << std::endl;
  polar << "          Solving response equations for perturbing field"
        << std::endl;
  polar << "          -----------------------------------------------"
        << std::endl;
  polar << "" << std::endl;
  polar << " number of frequencies:     1" << std::endl;
  polar << " frequency in a.u.:  0.0000000E+00" << std::endl;
  polar << " Perturbing field: electric" << std::endl;
  polar << " Using El. Dipole Length Gauge" << std::endl;
  polar << "" << std::endl;
  polar << " Setting up CPKS with frequency omega =      0.00000000 a.u."
        << std::endl;
  polar << "" << std::endl;
  polar << " STATIC response" << std::endl;
  polar << "" << std::endl;
  polar << "                                NWChem CPHF Module" << std::endl;
  polar << "                                ------------------" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "  scftype          =     RHF " << std::endl;
  polar << "  nclosed          =        5" << std::endl;
  polar << "  nopen            =        0" << std::endl;
  polar << "  variables        =       60" << std::endl;
  polar << "  # of vectors     =        3" << std::endl;
  polar << "  tolerance        = 0.10D-03" << std::endl;
  polar << "  level shift      = 0.00D+00" << std::endl;
  polar << "  max iterations   =       50" << std::endl;
  polar << "  max subspace     =       30" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " #quartets = 4.186D+03 #integrals = 1.171D+04 #direct =  0.0% "
           "#cached =100.0%"
        << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " Integral file          = ./system.aoints.0" << std::endl;
  polar << " Record size in doubles =  65536        No. of integs per rec  =  "
           "43688"
        << std::endl;
  polar << " Max. records in memory =      2        Max. records in file   = "
           "826900"
        << std::endl;
  polar << " No. of bits per label  =      8        No. of bits per value  =   "
           "  64"
        << std::endl;
  polar << "" << std::endl;
  polar << " SCF residual:    1.0605013831017603E-008" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << "Iterative solution of linear equations" << std::endl;
  polar << "  No. of variables       60" << std::endl;
  polar << "  No. of equations        3" << std::endl;
  polar << "  Maximum subspace       30" << std::endl;
  polar << "        Iterations       50" << std::endl;
  polar << "       Convergence  1.0D-04" << std::endl;
  polar << "        Start time      0.5" << std::endl;
  polar << "" << std::endl;
  polar << "  IO offset    32.000000000000000     " << std::endl;
  polar << "  IO error message >End of File" << std::endl;
  polar << "  file_read_ga: failing reading from ./system.cphf_sol"
        << std::endl;
  polar << " Error in restart solution: ./system.cphf_sol                      "
           "                                                                   "
           "                                                                   "
           "                                                                   "
           "               "
        << std::endl;
  polar << "" << std::endl;
  polar << "   iter   nsub   residual    time" << std::endl;
  polar << "   ----  ------  --------  ---------" << std::endl;
  polar << "     1      3    6.75D-01       0.7" << std::endl;
  polar << "     2      6    7.75D-02       1.0" << std::endl;
  polar << "     3      9    1.66D-03       1.3" << std::endl;
  polar << "     4     12    1.77D-04       1.6" << std::endl;
  polar << "     5     15    9.21D-05       1.8" << std::endl;
  polar << "" << std::endl;
  polar << " Parallel integral file used       1 records with       0 large "
           "values"
        << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " Electric Dipole Response Matrix (nonzero elements):" << std::endl;
  polar << "" << std::endl;
  polar << "              1        2        3      " << std::endl;
  polar << "    1   11.4169  -0.0003   0.0001" << std::endl;
  polar << "    2   -0.0003  11.4172  -0.0002" << std::endl;
  polar << "    3    0.0001  -0.0002  11.4174" << std::endl;
  polar << "" << std::endl;
  polar << " ------------------------------------------" << std::endl;
  polar << " average:        11.41717 + I       0.00000" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " DFT Linear Response polarizability / au " << std::endl;
  polar << " Frequency  =       0.0000000 / au" << std::endl;
  polar << " Wavelength = -999999.0000000 / nm" << std::endl;
  polar << "              X              Y              Z" << std::endl;
  polar << " -----------------------------------------------" << std::endl;
  polar << " X      11.4168971     -0.0002704      0.0000932" << std::endl;
  polar << " Y      -0.0002704     11.4171575     -0.0002280" << std::endl;
  polar << " Z       0.0000932     -0.0002280     11.4174457" << std::endl;
  polar << " -----------------------------------------------" << std::endl;
  polar << " Eigenvalues =      11.4167244     11.4171385     11.4176375"
        << std::endl;
  polar << " Isotropic   =      11.4171668" << std::endl;
  polar << " Anisotropic =       0.0006466" << std::endl;
  polar << " -----------------------------------------------" << std::endl;
  polar << "" << std::endl;
  polar << " Magnetic Dipole Response Matrix (nonzero elements):" << std::endl;
  polar << " Optical rotation tensor Beta" << std::endl;
  polar << "" << std::endl;
  polar << " *** static G'=0. Use GPRIME input with finite omega" << std::endl;
  polar << " *** or better ORBETA without GPRIME for small or zero freq."
        << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " zero matrix" << std::endl;
  polar << "" << std::endl;
  polar << "" << std::endl;
  polar << " Total times  cpu:        2.2s     wall:        1.8s" << std::endl;

  polar.close();
  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> nwchem =
      std::unique_ptr<QMPackage>(QMPackages().Create("nwchem"));
  Logger log;
  nwchem->setLog(&log);
  nwchem->setRunDir(".");
  nwchem->setLogFileName("polar_nwchem.log");
  Eigen::Matrix3d polar_mat = nwchem->GetPolarizability();

  Eigen::Matrix3d polar_ref = Eigen::Matrix3d::Zero();
  polar_ref << 11.4168971, -0.0002704, 0.0000932, -0.0002704, 11.4171575,
      -0.0002280, 0.0000932, -0.0002280, 11.4174457;

  bool polar_check = polar_ref.isApprox(polar_mat, 1e-5);
  if (!polar_check) {
    std::cout << "res" << std::endl;
    std::cout << polar_mat << std::endl;
    std::cout << "ref " << std::endl;
    std::cout << polar_ref << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(ext_charges) {
  std::ofstream ext_charge("extcharges_nwchem.log");

  ext_charge << " argument  1 = system.nw" << std::endl;
  ext_charge << "                                         " << std::endl;
  ext_charge << "                                         " << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge
      << "              Northwest Computational Chemistry Package (NWChem) 6.8"
      << std::endl;
  ext_charge
      << "              ------------------------------------------------------"
      << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "           Job information" << std::endl;
  ext_charge << "           ---------------" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "    hostname        = nbwin1568" << std::endl;
  ext_charge << "    program         = nwchem" << std::endl;
  ext_charge << "    date            = Tue Jul 30 20:18:19 2019" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "    compiled        = Wed_Aug_15_19:14:19_2018" << std::endl;
  ext_charge
      << "    source          = /home/edo/debichem-team/nwchem/nwchem-6.8.1"
      << std::endl;
  ext_charge << "    nwchem branch   = 6.8.1" << std::endl;
  ext_charge << "    nwchem revision = v6.8-133-ge032219" << std::endl;
  ext_charge << "    ga revision     = 5.6.5" << std::endl;
  ext_charge << "    use scalapack   = T" << std::endl;
  ext_charge << "    input           = system.nw" << std::endl;
  ext_charge << "    prefix          = system." << std::endl;
  ext_charge << "    data base       = ./system.db" << std::endl;
  ext_charge << "    status          = restart" << std::endl;
  ext_charge << "    nproc           =        1" << std::endl;
  ext_charge << "    time left       =     -1s" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "           Directory information" << std::endl;
  ext_charge << "           ---------------------" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "  0 permanent = ." << std::endl;
  ext_charge << "  0 scratch   = ." << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "           Geometries in the database" << std::endl;
  ext_charge << "           --------------------------" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "       Name                              Natoms  Last Modified"
             << std::endl;
  ext_charge << "       --------------------------------  ------  "
                "------------------------"
             << std::endl;
  ext_charge << "    1  geometry                               5  Tue Jul 30 "
                "20:16:16 2019  "
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "  The geometry named \"geometry\" is the default for restart"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "           Basis sets in the database" << std::endl;
  ext_charge << "           --------------------------" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "       Name                              Natoms  Last Modified"
             << std::endl;
  ext_charge << "        --------------------------------  ------  "
                "------------------------"
             << std::endl;
  ext_charge << "    1  ao basis                               0  Tue Jul 30 "
                "20:16:16 2019  "
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "  The basis set named \"ao basis\" is the default AO basis "
                "for restart"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "                                NWChem Input Module"
             << std::endl;
  ext_charge << "                                -------------------"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Scaling coordinates for geometry \"geometry\" by  1.889725989"
             << std::endl;
  ext_charge << " (inverse scale =  0.529177249)" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "                             Geometry \"geometry\" -> \"\""
             << std::endl;
  ext_charge << "                             -------------------------"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Output coordinates in angstroms (scale by  1.889725989 to "
                "convert to a.u.)"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "  No.       Tag          Charge          X              Y     "
                "         Z"
             << std::endl;
  ext_charge << " ---- ---------------- ---------- -------------- "
                "-------------- --------------"
             << std::endl;
  ext_charge << "    1 H                    1.0000     0.52881000     "
                "0.16101000     0.93592000"
             << std::endl;
  ext_charge << "    2 C                    6.0000     0.00001000     "
                "0.00001000     0.00002000"
             << std::endl;
  ext_charge << "    3 H                    1.0000     0.20511000     "
                "0.82401000    -0.67858000"
             << std::endl;
  ext_charge << "    4 H                    1.0000     0.33451000    "
                "-0.93139000    -0.44958000"
             << std::endl;
  ext_charge << "    5 H                    1.0000    -1.06849000    "
                "-0.05369000     0.19212000"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "      Atomic Mass " << std::endl;
  ext_charge << "      ----------- " << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "      H                  1.007825" << std::endl;
  ext_charge << "      C                 12.000000" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Effective nuclear repulsion energy (a.u.)      13.4728439460"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "            Nuclear Dipole moment (a.u.) " << std::endl;
  ext_charge << "            ----------------------------" << std::endl;
  ext_charge << "        X                 Y               Z" << std::endl;
  ext_charge << " ---------------- ---------------- ----------------"
             << std::endl;
  ext_charge << "     0.0000000000    -0.0000000000     0.0000000000"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "            XYZ format geometry" << std::endl;
  ext_charge << "            -------------------" << std::endl;
  ext_charge << "     5" << std::endl;
  ext_charge << " geometry" << std::endl;
  ext_charge
      << " H                     0.52881000     0.16101000     0.93592000"
      << std::endl;
  ext_charge
      << " C                     0.00001000     0.00001000     0.00002000"
      << std::endl;
  ext_charge
      << " H                     0.20511000     0.82401000    -0.67858000"
      << std::endl;
  ext_charge
      << " H                     0.33451000    -0.93139000    -0.44958000"
      << std::endl;
  ext_charge
      << " H                    -1.06849000    -0.05369000     0.19212000"
      << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Summary of \"ao basis\" -> \"\" (cartesian)" << std::endl;
  ext_charge << " -------------------------------------------------------------"
                "-----------------"
             << std::endl;
  ext_charge << "       Tag                 Description            Shells   "
                "Functions and Types"
             << std::endl;
  ext_charge << " ---------------- ------------------------------  ------  "
                "---------------------"
             << std::endl;
  ext_charge
      << " *                           3-21G                    on all atoms "
      << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Successfully loaded bq information from bq.xyz" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "                   Bq Structure Information (Angstroms)"
             << std::endl;
  ext_charge << "                   ------------------------------------"
             << std::endl;
  ext_charge << " Name: default                                                "
                "                 "
             << std::endl;
  ext_charge << " Number of centers:                     2" << std::endl;
  ext_charge << "    1  Bq         5.00000000     0.00000000     0.00000000   "
                "charge        1.00000000"
             << std::endl;
  ext_charge << "    2  Bq        -5.00000000     0.00000000     0.00000000   "
                "charge       -1.00000000"
             << std::endl;
  ext_charge << " Total Bq charge:    0.0000000000000000     " << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "                                 NWChem DFT Module"
             << std::endl;
  ext_charge << "                                 -----------------"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;

  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "  Caching 1-el integrals " << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "            General Information" << std::endl;
  ext_charge << "            -------------------" << std::endl;
  ext_charge << "          SCF calculation type: DFT" << std::endl;
  ext_charge << "          Wavefunction type:  closed shell." << std::endl;
  ext_charge << "          No. of atoms     :     5" << std::endl;
  ext_charge << "          No. of electrons :    10" << std::endl;
  ext_charge << "           Alpha electrons :     5" << std::endl;
  ext_charge << "            Beta electrons :     5" << std::endl;
  ext_charge << "          Charge           :     0" << std::endl;
  ext_charge << "          Spin multiplicity:     1" << std::endl;
  ext_charge << "          Use of symmetry is: off; symmetry adaption is: off"
             << std::endl;
  ext_charge << "          Maximum number of iterations:  30" << std::endl;
  ext_charge << "          AO basis - number of functions:    17" << std::endl;
  ext_charge << "                     number of shells:    13" << std::endl;
  ext_charge << "          Convergence on energy requested:  1.00D-06"
             << std::endl;
  ext_charge << "          Convergence on density requested:  1.00D-05"
             << std::endl;
  ext_charge << "          Convergence on gradient requested:  5.00D-04"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "              XC Information" << std::endl;
  ext_charge << "              --------------" << std::endl;
  ext_charge << "                         PBE0 Method XC Functional"
             << std::endl;
  ext_charge
      << "                     Hartree-Fock (Exact) Exchange  0.250          "
      << std::endl;
  ext_charge
      << "          PerdewBurkeErnzerhof Exchange Functional  0.750          "
      << std::endl;
  ext_charge
      << "            Perdew 1991 LDA Correlation Functional  1.000 local    "
      << std::endl;
  ext_charge
      << "           PerdewBurkeErnz. Correlation Functional  1.000 non-local"
      << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;

  ext_charge << " Nuclear repulsion energy =   13.472843946004776     "
             << std::endl;
  ext_charge << " Bq nuclear interaction energy =  -4.3475130603090328E-003"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "      Superposition of Atomic Density Guess" << std::endl;
  ext_charge << "      -------------------------------------" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Sum of atomic energies:         -39.44855565" << std::endl;
  ext_charge << " Nuclear repulsion energy =   13.472843946004776     "
             << std::endl;
  ext_charge << " Bq nuclear interaction energy =  -4.3475130603090328E-003"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "      Non-variational initial energy" << std::endl;
  ext_charge << "      ------------------------------" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Total energy =     -40.025908" << std::endl;
  ext_charge << " 1-e energy   =     -78.205184" << std::endl;
  ext_charge << " 2-e energy   =      24.710780" << std::endl;
  ext_charge << " HOMO         =      -0.493492" << std::endl;
  ext_charge << " LUMO         =       0.190711" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Nuclear repulsion energy =   13.472843946004776     "
             << std::endl;
  ext_charge << " Bq nuclear interaction energy =  -4.3475130603090328E-003"
             << std::endl;
  ext_charge << "   Time after variat. SCF:      0.0" << std::endl;
  ext_charge << "   Time prior to 1st pass:      0.0" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " #quartets = 4.186D+03 #integrals = 1.171D+04 #direct =  0.0% "
                "#cached =100.0%"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Integral file          = ./system.aoints.0" << std::endl;
  ext_charge << " Record size in doubles =  65536        No. of integs per rec "
                " =  43688"
             << std::endl;
  ext_charge << " Max. records in memory =      2        Max. records in file  "
                " = 851951"
             << std::endl;
  ext_charge << " No. of bits per label  =      8        No. of bits per value "
                " =     64"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Nuclear repulsion energy =   13.472843946004776     "
             << std::endl;
  ext_charge << " Bq nuclear interaction energy =  -4.3475130603090328E-003"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "         Total DFT energy =      -40.244517738916"
             << std::endl;
  ext_charge << "      One electron energy =      -79.606936158045"
             << std::endl;
  ext_charge << "           Coulomb energy =       32.720789921507"
             << std::endl;
  ext_charge << "    Exchange-Corr. energy =       -6.826867935323"
             << std::endl;
  ext_charge << " Nuclear repulsion energy =       13.468496432944"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Numeric. integr. density =       10.000000179004"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "     Total iterative time =      0.5s" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " Nuclear repulsion energy =   13.472843946004776     "
             << std::endl;
  ext_charge << " Bq nuclear interaction energy =  -4.3475130603090328E-003"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << " center of mass" << std::endl;
  ext_charge << " --------------" << std::endl;
  ext_charge << " x =   0.00000702 y =   0.00000702 z =   0.00001403"
             << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "" << std::endl;
  ext_charge << "                                NWChem Input Module"
             << std::endl;
  ext_charge << "                                -------------------"
             << std::endl;
  ext_charge << " Total times  cpu:        1.4s     wall:        0.6s"
             << std::endl;
  ext_charge.close();

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> nwchem =
      std::unique_ptr<QMPackage>(QMPackages().Create("nwchem"));
  Logger log;
  Orbitals orb;
  nwchem->setLog(&log);
  nwchem->setRunDir(".");
  nwchem->setLogFileName("extcharges_nwchem.log");
  nwchem->ParseLogFile(orb);
  const QMMolecule& seg = orb.QMAtoms();
  double ang2bohr = votca::tools::conv::ang2bohr;
  QMMolecule ref("ref", 0);
  Eigen::Vector3d pos1 = {0.528800, 0.161000, 0.935900};
  QMAtom s1(0, "H", pos1 * ang2bohr);
  pos1 = {0.00001000, 0.00001000, 0.00002000};
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
  BOOST_CHECK_CLOSE(orb.getScaHFX(), 0.25, 1e-5);

  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (int i = 0; i < seg.size(); i++) {
    bool coord_check = ref[i].getPos().isApprox(seg[i].getPos(), 1e-4);
    BOOST_CHECK_EQUAL(coord_check, true);
    if (!coord_check) {
      std::cout << "result coord " << i << std::endl;
      std::cout << seg[i].getPos().transpose() * votca::tools::conv::bohr2ang
                << std::endl;
      std::cout << "ref coord " << i << std::endl;
      std::cout << ref[i].getPos().transpose() * votca::tools::conv::bohr2ang
                << std::endl;
    }
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
  double ref_tot = -40.244517738916;  // HF - Self energy ext charges
  BOOST_CHECK_CLOSE(orb.getDFTTotalEnergy(), ref_tot, 1e-5);
}

BOOST_AUTO_TEST_CASE(getcharges) {
  std::ofstream chelpg("charges_nwchem.log");
  chelpg << " argument  1 = system.nw" << std::endl;
  chelpg << "                                         " << std::endl;
  chelpg << "                                         " << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg
      << "              Northwest Computational Chemistry Package (NWChem) 6.8"
      << std::endl;
  chelpg
      << "              ------------------------------------------------------"
      << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                             Geometry \"geometry\" -> \"\""
         << std::endl;
  chelpg << "                             -------------------------"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Output coordinates in angstroms (scale by  1.889725989 to "
            "convert to a.u.)"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "  No.       Tag          Charge          X              Y         "
            "     Z"
         << std::endl;
  chelpg << " ---- ---------------- ---------- -------------- -------------- "
            "--------------"
         << std::endl;
  chelpg << "    1 H                    1.0000     0.52881000     0.16101000   "
            "  0.93592000"
         << std::endl;
  chelpg << "    2 C                    6.0000     0.00001000     0.00001000   "
            "  0.00002000"
         << std::endl;
  chelpg << "    3 H                    1.0000     0.20511000     0.82401000   "
            " -0.67858000"
         << std::endl;
  chelpg << "    4 H                    1.0000     0.33451000    -0.93139000   "
            " -0.44958000"
         << std::endl;
  chelpg << "    5 H                    1.0000    -1.06849000    -0.05369000   "
            "  0.19212000"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "      Atomic Mass " << std::endl;
  chelpg << "      ----------- " << std::endl;
  chelpg << "" << std::endl;
  chelpg << "      H                  1.007825" << std::endl;
  chelpg << "      C                 12.000000" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Effective nuclear repulsion energy (a.u.)      13.4728439460"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "            Nuclear Dipole moment (a.u.) " << std::endl;
  chelpg << "            ----------------------------" << std::endl;
  chelpg << "        X                 Y               Z" << std::endl;
  chelpg << " ---------------- ---------------- ----------------" << std::endl;
  chelpg << "     0.0000000000    -0.0000000000     0.0000000000" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "            XYZ format geometry" << std::endl;
  chelpg << "            -------------------" << std::endl;
  chelpg << "     5" << std::endl;
  chelpg << " geometry" << std::endl;
  chelpg << " H                     0.52881000     0.16101000     0.93592000"
         << std::endl;
  chelpg << " C                     0.00001000     0.00001000     0.00002000"
         << std::endl;
  chelpg << " H                     0.20511000     0.82401000    -0.67858000"
         << std::endl;
  chelpg << " H                     0.33451000    -0.93139000    -0.44958000"
         << std::endl;
  chelpg << " H                    -1.06849000    -0.05369000     0.19212000"
         << std::endl;
  chelpg << " Summary of \"ao basis\" -> \"\" (cartesian)" << std::endl;
  chelpg << " -----------------------------------------------------------------"
            "-------------"
         << std::endl;
  chelpg << "       Tag                 Description            Shells   "
            "Functions and Types"
         << std::endl;
  chelpg << " ---------------- ------------------------------  ------  "
            "---------------------"
         << std::endl;
  chelpg
      << " *                           3-21G                    on all atoms "
      << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                                 NWChem DFT Module" << std::endl;
  chelpg << "                                 -----------------" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Summary of \"ao basis\" -> \"ao basis\" (cartesian)" << std::endl;
  chelpg << " -----------------------------------------------------------------"
            "-------------"
         << std::endl;
  chelpg << "       Tag                 Description            Shells   "
            "Functions and Types"
         << std::endl;
  chelpg << " ---------------- ------------------------------  ------  "
            "---------------------"
         << std::endl;
  chelpg
      << " H                           3-21G                   2        2   2s"
      << std::endl;
  chelpg << " C                           3-21G                   5        9   "
            "3s2p"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Summary of \"ao basis\" -> \"ao basis\" (cartesian)" << std::endl;
  chelpg << " -----------------------------------------------------------------"
            "-------------"
         << std::endl;
  chelpg << "       Tag                 Description            Shells   "
            "Functions and Types"
         << std::endl;
  chelpg << " ---------------- ------------------------------  ------  "
            "---------------------"
         << std::endl;
  chelpg
      << " H                           3-21G                   2        2   2s"
      << std::endl;
  chelpg << " C                           3-21G                   5        9   "
            "3s2p"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "  Caching 1-el integrals " << std::endl;
  chelpg << "" << std::endl;
  chelpg << "            General Information" << std::endl;
  chelpg << "            -------------------" << std::endl;
  chelpg << "          SCF calculation type: DFT" << std::endl;
  chelpg << "          Wavefunction type:  closed shell." << std::endl;
  chelpg << "          No. of atoms     :     5" << std::endl;
  chelpg << "          No. of electrons :    10" << std::endl;
  chelpg << "           Alpha electrons :     5" << std::endl;
  chelpg << "            Beta electrons :     5" << std::endl;
  chelpg << "          Charge           :     0" << std::endl;
  chelpg << "          Spin multiplicity:     1" << std::endl;
  chelpg << "          Use of symmetry is: off; symmetry adaption is: off"
         << std::endl;
  chelpg << "          Maximum number of iterations:  30" << std::endl;
  chelpg << "          AO basis - number of functions:    17" << std::endl;
  chelpg << "                     number of shells:    13" << std::endl;
  chelpg << "          Convergence on energy requested:  1.00D-06" << std::endl;
  chelpg << "          Convergence on density requested:  1.00D-05"
         << std::endl;
  chelpg << "          Convergence on gradient requested:  5.00D-04"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "              XC Information" << std::endl;
  chelpg << "              --------------" << std::endl;
  chelpg << "                         PBE0 Method XC Functional" << std::endl;
  chelpg
      << "                     Hartree-Fock (Exact) Exchange  0.250          "
      << std::endl;
  chelpg
      << "          PerdewBurkeErnzerhof Exchange Functional  0.750          "
      << std::endl;
  chelpg
      << "            Perdew 1991 LDA Correlation Functional  1.000 local    "
      << std::endl;
  chelpg
      << "           PerdewBurkeErnz. Correlation Functional  1.000 non-local"
      << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                    Damping( 0%)  Levelshifting(0.5)       DIIS"
         << std::endl;
  chelpg
      << "                  --------------- ------------------- ---------------"
      << std::endl;
  chelpg << "          dE  on:    start            ASAP                start   "
         << std::endl;
  chelpg << "          dE off:    2 iters         30 iters            30 iters "
         << std::endl;
  chelpg << "" << std::endl;

  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Loading old vectors from job with title :" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "   Time after variat. SCF:      0.0" << std::endl;
  chelpg << "   Time prior to 1st pass:      0.0" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " #quartets = 4.186D+03 #integrals = 1.171D+04 #direct =  0.0% "
            "#cached =100.0%"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Integral file          = ./system.aoints.0" << std::endl;
  chelpg << " Record size in doubles =  65536        No. of integs per rec  =  "
            "43688"
         << std::endl;
  chelpg << " Max. records in memory =      2        Max. records in file   = "
            "851953"
         << std::endl;
  chelpg << " No. of bits per label  =      8        No. of bits per value  =  "
            "   64"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Grid_pts file          = ./system.gridpts.0" << std::endl;
  chelpg << " Record size in doubles =  12289        No. of grid_pts per rec  "
            "=   3070"
         << std::endl;
  chelpg << " Max. records in memory =     37        Max. recs in file   =   "
            "4543386"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "           Memory utilization after 1st SCF pass: " << std::endl;
  chelpg
      << "           Heap Space remaining (MW):       12.52            12520080"
      << std::endl;
  chelpg
      << "          Stack Space remaining (MW):       13.11            13106932"
      << std::endl;
  chelpg << "" << std::endl;
  chelpg << "   convergence    iter        energy       DeltaE   RMS-Dens  "
            "Diis-err    time"
         << std::endl;
  chelpg << " ---------------- ----- ----------------- --------- --------- "
            "---------  ------"
         << std::endl;
  chelpg << " d= 0,ls=0.0,diis     1    -40.2408052356 -5.37D+01  1.30D-07  "
            "4.19D-12     0.2"
         << std::endl;
  chelpg << " d= 0,ls=0.0,diis     2    -40.2408052356 -5.68D-14  6.43D-08  "
            "3.25D-12     0.3"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "         Total DFT energy =      -40.240805235609" << std::endl;
  chelpg << "      One electron energy =      -79.605210629040" << std::endl;
  chelpg << "           Coulomb energy =       32.718686983262" << std::endl;
  chelpg << "    Exchange-Corr. energy =       -6.827125535836" << std::endl;
  chelpg << " Nuclear repulsion energy =       13.472843946005" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Numeric. integr. density =       10.000000095312" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "     Total iterative time =      0.3s" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                       DFT Final Molecular Orbital Analysis"
         << std::endl;
  chelpg << "                       ------------------------------------"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " center of mass" << std::endl;
  chelpg << " --------------" << std::endl;
  chelpg << " x =   0.00000702 y =   0.00000702 z =   0.00001403" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " moments of inertia (a.u.)" << std::endl;
  chelpg << " ------------------" << std::endl;
  chelpg << "          11.339523676494           0.000128306894          "
            "-0.000262759041"
         << std::endl;
  chelpg << "           0.000128306894          11.339492365125           "
            "0.000161348049"
         << std::endl;
  chelpg << "          -0.000262759041           0.000161348049          "
            "11.338906849319"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "     Multipole analysis of the density" << std::endl;
  chelpg << "     ---------------------------------" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "     L   x y z        total         alpha         beta         "
            "nuclear"
         << std::endl;
  chelpg << "     -   - - -        -----         -----         ----         "
            "-------"
         << std::endl;
  chelpg << "     0   0 0 0     -0.000000     -5.000000     -5.000000     "
            "10.000000"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "     1   1 0 0     -0.000025     -0.000012     -0.000012      "
            "0.000000"
         << std::endl;
  chelpg << "     1   0 1 0     -0.000018     -0.000009     -0.000009     "
            "-0.000000"
         << std::endl;
  chelpg << "     1   0 0 1     -0.000024     -0.000012     -0.000012      "
            "0.000000"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "     2   2 0 0     -6.111637     -5.868528     -5.868528      "
            "5.625419"
         << std::endl;
  chelpg << "     2   1 1 0     -0.000031      0.000048      0.000048     "
            "-0.000127"
         << std::endl;
  chelpg << "     2   1 0 1      0.000122     -0.000070     -0.000070      "
            "0.000261"
         << std::endl;
  chelpg << "     2   0 2 0     -6.111696     -5.868573     -5.868573      "
            "5.625450"
         << std::endl;
  chelpg << "     2   0 1 1     -0.000045      0.000058      0.000058     "
            "-0.000160"
         << std::endl;
  chelpg << "     2   0 0 2     -6.111425     -5.868728     -5.868728      "
            "5.626031"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Parallel integral file used       1 records with       0 large "
            "values"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Task  times  cpu:        0.7s     wall:        0.3s" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                                NWChem Input Module" << std::endl;
  chelpg << "                                -------------------" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                     NWChem Electrostatic Potential Fit Module"
         << std::endl;
  chelpg << "                     -----------------------------------------"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Atom parameters" << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Number of atoms is                                    5"
         << std::endl;
  chelpg << " Number of basis functions is                         17"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Atomic radii" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "    1    0.100000" << std::endl;
  chelpg << "    6    0.147000" << std::endl;
  chelpg << " FASTESP  F" << std::endl;
  chelpg << "  using M.O. file = ./system.movecs" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "" << std::endl;
  chelpg << "    Atom              Coordinates                           Charge"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << "                                                  ESP   "
         << std::endl;
  chelpg << "                                                        "
         << std::endl;
  chelpg << " " << std::endl;
  chelpg << "    1 H     0.052881    0.016101    0.093592    0.116609"
         << std::endl;
  chelpg << "    2 C     0.000001    0.000001    0.000002   -0.473785"
         << std::endl;
  chelpg << "    3 H     0.020511    0.082401   -0.067858    0.116584"
         << std::endl;
  chelpg << "    4 H     0.033451   -0.093139   -0.044958    0.119898"
         << std::endl;
  chelpg << "    5 H    -0.106849   -0.005369    0.019212    0.120694"
         << std::endl;
  chelpg << "                                            ------------"
         << std::endl;
  chelpg << "                                                0.000000"
         << std::endl;
  chelpg << "" << std::endl;
  chelpg << " Task  times  cpu:        0.1s     wall:        0.0s" << std::endl;
  chelpg << " Summary of allocated global arrays" << std::endl;
  chelpg << "-----------------------------------" << std::endl;
  chelpg << " Total times  cpu:        0.8s     wall:        0.4s" << std::endl;

  chelpg.close();
  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> nwchem =
      std::unique_ptr<QMPackage>(QMPackages().Create("nwchem"));
  Logger log;
  nwchem->setLog(&log);
  nwchem->setRunDir(".");
  nwchem->setLogFileName("charges_nwchem.log");
  StaticSegment seg = nwchem->GetCharges();

  double ang2bohr = votca::tools::conv::ang2bohr;
  StaticSegment ref("ref", 0);
  Eigen::Vector3d pos1 = {0.052881, 0.016101, 0.093592};
  StaticSite s1(0, "H", pos1 * ang2bohr);
  s1.setCharge(0.116609);
  pos1 = {0.000001, 0.000001, 0.000002};
  StaticSite s2(1, "C", pos1 * ang2bohr);
  s2.setCharge(-0.473785);
  pos1 = {0.020511, 0.082401, -0.067858};
  StaticSite s3(2, "H", pos1 * ang2bohr);
  s3.setCharge(0.116584);
  pos1 = {0.033451, -0.093139, -0.044958};
  StaticSite s4(3, "H", pos1 * ang2bohr);
  s4.setCharge(0.119898);
  pos1 = {-0.106849, -0.005369, 0.019212};
  StaticSite s5(4, "H", pos1 * ang2bohr);
  s5.setCharge(0.120694);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);

  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (int i = 0; i < seg.size(); i++) {
    BOOST_CHECK_EQUAL(ref[i].Q().isApprox(seg[i].Q(), 1e-5), true);
    bool check_pos = ref[i].getPos().isApprox(seg[i].getPos(), 1e-5);
    BOOST_CHECK_EQUAL(check_pos, true);
    if (!check_pos) {
      std::cout << "res" << i << std::endl;
      std::cout << seg[i].getPos().transpose() << std::endl;
      std::cout << "ref " << i << std::endl;
      std::cout << ref[i].getPos().transpose() << std::endl;
    }
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
}
BOOST_AUTO_TEST_CASE(opt_test) {
  std::ofstream opt("opt_nwchem.log");
  opt << " argument  1 = system.nw" << std::endl;
  opt << "                                         " << std::endl;
  opt << "                                         " << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "              Northwest Computational Chemistry Package (NWChem) 6.8"
      << std::endl;
  opt << "              ------------------------------------------------------"
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "           Job information" << std::endl;
  opt << "           ---------------" << std::endl;
  opt << "" << std::endl;
  opt << "    hostname        = nbwin1568" << std::endl;
  opt << "    program         = nwchem" << std::endl;
  opt << "    date            = Tue Jul 30 19:54:24 2019" << std::endl;
  opt << "" << std::endl;
  opt << "    compiled        = Wed_Aug_15_19:14:19_2018" << std::endl;
  opt << "    source          = /home/edo/debichem-team/nwchem/nwchem-6.8.1"
      << std::endl;
  opt << "    nwchem branch   = 6.8.1" << std::endl;
  opt << "    nwchem revision = v6.8-133-ge032219" << std::endl;
  opt << "    ga revision     = 5.6.5" << std::endl;
  opt << "    use scalapack   = T" << std::endl;
  opt << "    input           = system.nw" << std::endl;
  opt << "    prefix          = system." << std::endl;
  opt << "    data base       = ./system.db" << std::endl;
  opt << "    status          = restart" << std::endl;
  opt << "    nproc           =        1" << std::endl;
  opt << "    time left       =     -1s" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "    Theory    = dft                             " << std::endl;
  opt << "    Operation = optimize                      " << std::endl;
  opt << "    Status    = ok                            " << std::endl;
  opt << "    Qmmm      = F" << std::endl;
  opt << "    Ignore    = F" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << " Scaling coordinates for geometry \"geometry\" by  1.889725989"
      << std::endl;
  opt << " (inverse scale =  0.529177249)" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "                             Geometry \"geometry\" -> \"\""
      << std::endl;
  opt << "                             -------------------------" << std::endl;
  opt << "" << std::endl;
  opt << " Output coordinates in angstroms (scale by  1.889725989 to convert "
         "to a.u.)"
      << std::endl;
  opt << "" << std::endl;
  opt << "  No.       Tag          Charge          X              Y            "
         "  Z"
      << std::endl;
  opt << " ---- ---------------- ---------- -------------- -------------- "
         "--------------"
      << std::endl;
  opt << "    1 H                    1.0000     0.55881000     0.16101000     "
         "0.93592000"
      << std::endl;
  opt << "    2 C                    6.0000     0.03001000     0.00001000     "
         "0.00002000"
      << std::endl;
  opt << "    3 H                    1.0000     0.23511000     0.82401000    "
         "-0.67858000"
      << std::endl;
  opt << "    4 H                    1.0000     0.36451000    -0.93139000    "
         "-0.44958000"
      << std::endl;
  opt << "    5 H                    1.0000    -1.33849000    -0.05369000     "
         "0.19212000"
      << std::endl;
  opt << "" << std::endl;
  opt << "      Atomic Mass " << std::endl;
  opt << "      ----------- " << std::endl;
  opt << "" << std::endl;
  opt << "      H                  1.007825" << std::endl;
  opt << "      C                 12.000000" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << " Effective nuclear repulsion energy (a.u.)      12.7378283140"
      << std::endl;
  opt << "" << std::endl;
  opt << "            Nuclear Dipole moment (a.u.) " << std::endl;
  opt << "            ----------------------------" << std::endl;
  opt << "        X                 Y               Z" << std::endl;
  opt << " ---------------- ---------------- ----------------" << std::endl;
  opt << "     0.0000000000    -0.0000000000     0.0000000000" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "            XYZ format geometry" << std::endl;
  opt << "            -------------------" << std::endl;
  opt << "     5" << std::endl;
  opt << " geometry" << std::endl;
  opt << " H                     0.55881000     0.16101000     0.93592000"
      << std::endl;
  opt << " C                     0.03001000     0.00001000     0.00002000"
      << std::endl;
  opt << " H                     0.23511000     0.82401000    -0.67858000"
      << std::endl;
  opt << " H                     0.36451000    -0.93139000    -0.44958000"
      << std::endl;
  opt << " H                    -1.33849000    -0.05369000     0.19212000"
      << std::endl;
  opt << "" << std::endl;
  opt << "  Caching 1-el integrals " << std::endl;
  opt << "" << std::endl;
  opt << "            General Information" << std::endl;
  opt << "            -------------------" << std::endl;
  opt << "          SCF calculation type: DFT" << std::endl;
  opt << "          Wavefunction type:  closed shell." << std::endl;
  opt << "          No. of atoms     :     5" << std::endl;
  opt << "          No. of electrons :    10" << std::endl;
  opt << "           Alpha electrons :     5" << std::endl;
  opt << "            Beta electrons :     5" << std::endl;
  opt << "          Charge           :     0" << std::endl;
  opt << "          Spin multiplicity:     1" << std::endl;
  opt << "          Use of symmetry is: off; symmetry adaption is: off"
      << std::endl;
  opt << "          Maximum number of iterations:  60" << std::endl;
  opt << "          This is a Direct SCF calculation." << std::endl;
  opt << "          AO basis - number of functions:    17" << std::endl;
  opt << "                     number of shells:    13" << std::endl;
  opt << "          Convergence on energy requested:  1.00D-06" << std::endl;
  opt << "          Convergence on density requested:  1.00D-05" << std::endl;
  opt << "          Convergence on gradient requested:  5.00D-04" << std::endl;
  opt << "" << std::endl;
  opt << "              XC Information" << std::endl;
  opt << "              --------------" << std::endl;
  opt << "                         PBE0 Method XC Functional" << std::endl;
  opt << "                     Hartree-Fock (Exact) Exchange  0.250          "
      << std::endl;
  opt << "          PerdewBurkeErnzerhof Exchange Functional  0.750          "
      << std::endl;
  opt << "            Perdew 1991 LDA Correlation Functional  1.000 local    "
      << std::endl;
  opt << "           PerdewBurkeErnz. Correlation Functional  1.000 non-local"
      << std::endl;
  opt << "" << std::endl;
  opt << "             Grid Information" << std::endl;
  opt << "             ----------------" << std::endl;
  opt << "          Grid used for XC integration:  medium    " << std::endl;
  opt << "          Radial quadrature: Mura-Knowles        " << std::endl;
  opt << "          Angular quadrature: Lebedev. " << std::endl;
  opt << "          Tag              B.-S. Rad. Rad. Pts. Rad. Cut. Ang. Pts."
      << std::endl;
  opt << "          ---              ---------- --------- --------- ---------"
      << std::endl;
  opt << "          H                   0.35       45           6.0       434"
      << std::endl;
  opt << "          C                   0.70       49           5.0       434"
      << std::endl;
  opt << "          Grid pruning is: on " << std::endl;
  opt << "          Number of quadrature shells:   229" << std::endl;
  opt << "          Spatial weights used:  Erf1" << std::endl;
  opt << "" << std::endl;
  opt << "          Convergence Information" << std::endl;
  opt << "          -----------------------" << std::endl;
  opt << "          Convergence aids based upon iterative change in "
      << std::endl;
  opt << "          total energy or number of iterations. " << std::endl;
  opt << "          Levelshifting, if invoked, occurs when the " << std::endl;
  opt << "          HOMO/LUMO gap drops below (HL_TOL):  1.00D-02" << std::endl;
  opt << "          DIIS, if invoked, will attempt to extrapolate "
      << std::endl;
  opt << "          using up to (NFOCK): 10 stored Fock matrices." << std::endl;
  opt << "" << std::endl;
  opt << "                    Damping( 0%)  Levelshifting(0.5)       DIIS"
      << std::endl;
  opt << "                  --------------- ------------------- ---------------"
      << std::endl;
  opt << "          dE  on:    start            ASAP                start   "
      << std::endl;
  opt << "          dE off:    2 iters         60 iters            60 iters "
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "      Screening Tolerance Information" << std::endl;
  opt << "      -------------------------------" << std::endl;
  opt << "          Density screening/tol_rho:  1.00D-10" << std::endl;
  opt << "          AO Gaussian exp screening on grid/accAOfunc:  14"
      << std::endl;
  opt << "          CD Gaussian exp screening on grid/accCDfunc:  20"
      << std::endl;
  opt << "          XC Gaussian exp screening on grid/accXCfunc:  20"
      << std::endl;
  opt << "          Schwarz screening/accCoul:  1.00D-08" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << " Loading old vectors from job with title :" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "   Time after variat. SCF:      4.8" << std::endl;
  opt << "   Time prior to 1st pass:      4.8" << std::endl;
  opt << "" << std::endl;
  opt << " Grid_pts file          = ./system.gridpts.0" << std::endl;
  opt << " Record size in doubles =  12289        No. of grid_pts per rec  =   "
         "3070"
      << std::endl;
  opt << " Max. records in memory =     37        Max. recs in file   =   "
         "4543376"
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "           Memory utilization after 1st SCF pass: " << std::endl;
  opt << "           Heap Space remaining (MW):       12.65            12651164"
      << std::endl;
  opt << "          Stack Space remaining (MW):       13.11            13106940"
      << std::endl;
  opt << "" << std::endl;
  opt << "   convergence    iter        energy       DeltaE   RMS-Dens  "
         "Diis-err    time"
      << std::endl;
  opt << " ---------------- ----- ----------------- --------- --------- "
         "---------  ------"
      << std::endl;
  opt << " d= 0,ls=0.0,diis     1    -40.2409134539 -5.36D+01  7.16D-05  "
         "1.86D-06     5.0"
      << std::endl;
  opt << " d= 0,ls=0.0,diis     2    -40.2409137270 -2.73D-07  2.15D-05  "
         "1.19D-07     5.1"
      << std::endl;
  opt << " d= 0,ls=0.0,diis     3    -40.2409137321 -5.11D-09  9.56D-06  "
         "8.59D-08     5.1"
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "         Total DFT energy =      -40.240913732145" << std::endl;
  opt << "      One electron energy =      -79.457186509192" << std::endl;
  opt << "           Coulomb energy =       32.639969221881" << std::endl;
  opt << "    Exchange-Corr. energy =       -6.817708659018" << std::endl;
  opt << " Nuclear repulsion energy =       13.394012214183" << std::endl;
  opt << "" << std::endl;
  opt << " Numeric. integr. density =        9.999999713353" << std::endl;
  opt << "" << std::endl;
  opt << "     Total iterative time =      0.3s" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "  Caching 1-el integrals " << std::endl;
  opt << "" << std::endl;
  opt << "            General Information" << std::endl;
  opt << "            -------------------" << std::endl;
  opt << "          SCF calculation type: DFT" << std::endl;
  opt << "          Wavefunction type:  closed shell." << std::endl;
  opt << "          No. of atoms     :     5" << std::endl;
  opt << "          No. of electrons :    10" << std::endl;
  opt << "           Alpha electrons :     5" << std::endl;
  opt << "            Beta electrons :     5" << std::endl;
  opt << "          Charge           :     0" << std::endl;
  opt << "          Spin multiplicity:     1" << std::endl;
  opt << "          Use of symmetry is: off; symmetry adaption is: off"
      << std::endl;
  opt << "          Maximum number of iterations:  60" << std::endl;
  opt << "          This is a Direct SCF calculation." << std::endl;
  opt << "          AO basis - number of functions:    17" << std::endl;
  opt << "                     number of shells:    13" << std::endl;
  opt << "          Convergence on energy requested:  1.00D-06" << std::endl;
  opt << "          Convergence on density requested:  1.00D-05" << std::endl;
  opt << "          Convergence on gradient requested:  5.00D-04" << std::endl;
  opt << "" << std::endl;
  opt << "              XC Information" << std::endl;
  opt << "              --------------" << std::endl;
  opt << "                         PBE0 Method XC Functional" << std::endl;
  opt << "                     Hartree-Fock (Exact) Exchange  0.250          "
      << std::endl;
  opt << "          PerdewBurkeErnzerhof Exchange Functional  0.750          "
      << std::endl;
  opt << "            Perdew 1991 LDA Correlation Functional  1.000 local    "
      << std::endl;
  opt << "           PerdewBurkeErnz. Correlation Functional  1.000 non-local"
      << std::endl;
  opt << "" << std::endl;

  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << " Effective nuclear repulsion energy (a.u.)      13.3931649813"
      << std::endl;
  opt << "" << std::endl;

  opt << "" << std::endl;
  opt << "  The DFT is already converged " << std::endl;
  opt << "" << std::endl;
  opt << "         Total DFT energy =    -40.240914039306" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "            General Information" << std::endl;
  opt << "            -------------------" << std::endl;
  opt << "          SCF calculation type: DFT" << std::endl;
  opt << "          Wavefunction type:  closed shell." << std::endl;
  opt << "          No. of atoms     :     5" << std::endl;
  opt << "          No. of electrons :    10" << std::endl;
  opt << "           Alpha electrons :     5" << std::endl;
  opt << "            Beta electrons :     5" << std::endl;
  opt << "          Charge           :     0" << std::endl;
  opt << "          Spin multiplicity:     1" << std::endl;
  opt << "          Use of symmetry is: off; symmetry adaption is: off"
      << std::endl;
  opt << "          Maximum number of iterations:  60" << std::endl;
  opt << "          This is a Direct SCF calculation." << std::endl;
  opt << "          AO basis - number of functions:    17" << std::endl;
  opt << "                     number of shells:    13" << std::endl;
  opt << "          Convergence on energy requested:  1.00D-06" << std::endl;
  opt << "          Convergence on density requested:  1.00D-05" << std::endl;
  opt << "          Convergence on gradient requested:  5.00D-04" << std::endl;
  opt << "" << std::endl;
  opt << "              XC Information" << std::endl;
  opt << "              --------------" << std::endl;
  opt << "                         PBE0 Method XC Functional" << std::endl;
  opt << "                     Hartree-Fock (Exact) Exchange  0.250          "
      << std::endl;
  opt << "          PerdewBurkeErnzerhof Exchange Functional  0.750          "
      << std::endl;
  opt << "            Perdew 1991 LDA Correlation Functional  1.000 local    "
      << std::endl;
  opt << "           PerdewBurkeErnz. Correlation Functional  1.000 non-local"
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "                            NWChem DFT Gradient Module" << std::endl;
  opt << "                            --------------------------" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "  charge          =   0.00" << std::endl;
  opt << "  wavefunction    = closed shell" << std::endl;
  opt << "" << std::endl;

  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "      ----------------------" << std::endl;
  opt << "      Optimization converged" << std::endl;
  opt << "      ----------------------" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "  Step       Energy      Delta E   Gmax     Grms     Xrms     Xmax   "
         "Walltime"
      << std::endl;
  opt << "  ---- ---------------- -------- -------- -------- -------- -------- "
         "--------"
      << std::endl;
  opt << "@    5     -40.24091404 -3.1D-07  0.00003  0.00002  0.00043  0.00075 "
         "     5.7"
      << std::endl;
  opt << "                                     ok       ok       ok       ok  "
      << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << "                         Geometry \"geometry\" -> \"geometry\""
      << std::endl;
  opt << "                         ---------------------------------"
      << std::endl;
  opt << "" << std::endl;
  opt << " Output coordinates in angstroms (scale by  1.889725989 to convert "
         "to a.u.)"
      << std::endl;
  opt << "" << std::endl;
  opt << "  No.       Tag          Charge          X              Y            "
         "  Z"
      << std::endl;
  opt << " ---- ---------------- ---------- -------------- -------------- "
         "--------------"
      << std::endl;
  opt << "    1 H                    1.0000     0.48487682     0.15901714     "
         "0.95134757"
      << std::endl;
  opt << "    2 C                    6.0000    -0.03003067    -0.00002910    "
         "-0.00001387"
      << std::endl;
  opt << "    3 H                    1.0000     0.19360402     0.82779797    "
         "-0.67842703"
      << std::endl;
  opt << "    4 H                    1.0000     0.31013087    -0.93869906    "
         "-0.44600443"
      << std::endl;
  opt << "    5 H                    1.0000    -1.10863104    -0.04813694     "
         "0.17299775"
      << std::endl;
  opt << "" << std::endl;
  opt << "      Atomic Mass " << std::endl;
  opt << "      ----------- " << std::endl;
  opt << "" << std::endl;
  opt << "" << std::endl;
  opt << " Total times  cpu:       12.0s     wall:        5.7s" << std::endl;

  opt.close();

  QMPackageFactory::RegisterAll();
  std::unique_ptr<QMPackage> nwchem =
      std::unique_ptr<QMPackage>(QMPackages().Create("nwchem"));
  Logger log;
  nwchem->setLog(&log);
  nwchem->setRunDir(".");
  Orbitals orb;
  nwchem->setLogFileName("opt_nwchem.log");
  nwchem->ParseLogFile(orb);
  double ref_tot = -40.24091404;  // HF - Self energy ext charges

  BOOST_CHECK_CLOSE(ref_tot, orb.getDFTTotalEnergy(), 1e-5);
  const QMMolecule& seg = orb.QMAtoms();
  double ang2bohr = votca::tools::conv::ang2bohr;
  QMMolecule ref("ref", 0);
  Eigen::Vector3d pos1 = {0.48487682, 0.15901714, 0.95134757};
  QMAtom s1(0, "H", pos1 * ang2bohr);
  pos1 = {-0.03003067, -0.00002910, -0.00001387};
  QMAtom s2(1, "C", pos1 * ang2bohr);
  pos1 = {0.19360402, 0.82779797, -0.67842703};
  QMAtom s3(2, "H", pos1 * ang2bohr);
  pos1 = {0.31013087, -0.93869906, -0.44600443};
  QMAtom s4(3, "H", pos1 * ang2bohr);
  pos1 = {-1.10863104, -0.04813694, 0.17299775};
  QMAtom s5(4, "H", pos1 * ang2bohr);
  ref.push_back(s1);
  ref.push_back(s2);
  ref.push_back(s3);
  ref.push_back(s4);
  ref.push_back(s5);
  BOOST_CHECK_EQUAL(seg.size(), ref.size());
  for (int i = 0; i < seg.size(); i++) {
    bool coord_check = ref[i].getPos().isApprox(seg[i].getPos(), 1e-5);
    BOOST_CHECK_EQUAL(coord_check, true);
    if (!coord_check) {
      std::cout << "result coord " << i << std::endl;
      std::cout << seg[i].getPos().transpose() * votca::tools::conv::bohr2ang
                << std::endl;
      std::cout << "ref coord " << i << std::endl;
      std::cout << ref[i].getPos().transpose() * votca::tools::conv::bohr2ang
                << std::endl;
    }
    BOOST_CHECK_EQUAL(ref[i].getElement(), seg[i].getElement());
  }
}

BOOST_AUTO_TEST_SUITE_END()
