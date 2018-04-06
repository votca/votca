/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE test_hdf5
#include <boost/test/unit_test.hpp>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>

BOOST_AUTO_TEST_SUITE(test_hdf5)
using namespace votca::xtp;
BOOST_AUTO_TEST_CASE(checkpoint_file_test) {
    votca::xtp::CheckpointFile cpf("xtp_testing.hdf5");

    // Write orbitals
    votca::xtp::Orbitals orbWrite;

    int basisSetSize = 17;
    int occupiedLevels = 4;
    int unoccupiedLevels = 13;
    int numElectrons = 12;


    Eigen::VectorXd moeTest = Eigen::VectorXd::Zero(17);
    Eigen::MatrixXd mocTest = Eigen::MatrixXd::Zero(17, 17);

    QMAtom atoms[100];
    std::vector<QMAtom*> atomsTest;

    for (size_t p = 0; p < 100; ++p)
        atomsTest.push_back(atoms+p);

    double qmEnergy = -2.1025e-3;

    std::string qmPackage = "NOPE";
    double selfEnergy = 3.14159e23;

    std::string dftBasis = "AWESOME basis*,, 2/.8";
    std::string auxBasis = "cos(theta) = pretty okay basis";

    int rpaMin = '?';
    int rpaMax = 1e3;

    int qpMin = -91091;
    int qpMax = 75918275;

    unsigned int bseVmin = -6019386;
    unsigned int bseVmax = 1092581;

    unsigned int bseCmin = 2718L;
    unsigned int bseCmax = 42;


    votca::xtp::MatrixXfd eh_dTest = votca::xtp::MatrixXfd::Zero(32, 290);
    votca::xtp::MatrixXfd eh_xTest = votca::xtp::MatrixXfd::Zero(3, 22);
    votca::xtp::VectorXfd BSE_singlet_energiesTest = votca::xtp::VectorXfd::Zero(25);
    Eigen::MatrixXd vxcTest = Eigen::MatrixXd::Zero(200,200);

    std::string someECP = "aye aye Cap'n";


    orbWrite.setBasisSetSize(basisSetSize);
    orbWrite.setNumberOfLevels(occupiedLevels, unoccupiedLevels);
    orbWrite.setNumberOfElectrons(numElectrons);
    orbWrite.MOEnergies() = moeTest;
    orbWrite.MOCoefficients() = mocTest;
    orbWrite.QMAtoms() = atomsTest;
    orbWrite.setECP(someECP);
    orbWrite.setBasisSetSize(17);
    orbWrite.setNumberOfLevels(4, 13);

    orbWrite.eh_d() = eh_dTest;
    orbWrite.eh_x() = eh_xTest;

    orbWrite.BSESingletEnergies() = BSE_singlet_energiesTest;

    orbWrite.AOVxc() = vxcTest;


    orbWrite.WriteToCpt(cpf, "Test Orbital");

    // Read Orbitals
    votca::xtp::Orbitals orbRead;


    BOOST_AUTO_TEST_SUITE_END()}
