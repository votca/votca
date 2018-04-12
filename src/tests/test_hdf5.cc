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

    unsigned int bseVmin = -6019386;
    unsigned int bseVmax = 1092581;

    unsigned int bseCmin = 2718L;
    unsigned int bseCmax = 42;

    double scaHfx = 3.14159;

    std::string bseType = "A+";


    Eigen::MatrixXd vxcTest = Eigen::MatrixXd::Zero(200,200);
    std::string someECP = "aye aye Cap'n";

    Eigen::MatrixXd QPpertEnergiesTest = Eigen::MatrixXd::Zero(31, 42);
    Eigen::MatrixXd QPdiagEnergiesTest = Eigen::VectorXd::Zero(21);
    Eigen::MatrixXd QPdiagCoefficientsTest = Eigen::MatrixXd::Identity(31, 42);


    MatrixXfd eh_dTest = MatrixXfd::Zero(32, 290);
    MatrixXfd eh_xTest = MatrixXfd::Zero(3, 22);
    VectorXfd BSESingletEnergiesTest = VectorXfd::Zero(25);
    MatrixXfd BSESingletCoefficientsTest = MatrixXfd::Zero(25, 38);
    MatrixXfd BSESingletCoefficientsARTest = MatrixXfd::Zero(42, 42);

    VectorXfd BSETripletEnergiesTest = VectorXfd::Zero(33);
    MatrixXfd BSETripletCoefficientsTest = MatrixXfd::Zero(33,31);

    std::vector <votca::tools::vec> transitionDipolesTest;
    for (size_t i =0; i < 1000; ++i){
        transitionDipolesTest.push_back(votca::tools::vec(1,2,3));
    }

    // Write orbitals
    Orbitals orbWrite;

    orbWrite.setBasisSetSize(basisSetSize);
    orbWrite.setNumberOfLevels(occupiedLevels, unoccupiedLevels);
    orbWrite.setNumberOfElectrons(numElectrons);
    orbWrite.MOEnergies() = moeTest;
    orbWrite.MOCoefficients() = mocTest;
    //orbWrite.QMAtoms() = atomsTest;
    orbWrite.setQMEnergy(qmEnergy);
    orbWrite.setQMpackage(qmPackage);
    orbWrite.setSelfEnergy(selfEnergy);
    orbWrite.setDFTbasis(dftBasis);
    orbWrite.setAuxbasis(auxBasis);
    orbWrite.setRPAindices(rpaMin, rpaMax);
    // no need to write qpmin, qpmax
    orbWrite.setBSEindices(bseVmin, bseVmax, bseCmin, bseCmax, 3);
    orbWrite.setScaHFX(scaHfx);
    orbWrite.setBSEtype(bseType);
    orbWrite.setECP(someECP);
    orbWrite.QPpertEnergies() = QPpertEnergiesTest;
    orbWrite.QPdiagEnergies() = QPdiagEnergiesTest;
    orbWrite.QPdiagCoefficients() = QPdiagCoefficientsTest;
    orbWrite.setBasisSetSize(17);
    orbWrite.setNumberOfLevels(4, 13);
    orbWrite.eh_d() = eh_dTest;
    orbWrite.eh_x() = eh_xTest;
    orbWrite.BSESingletEnergies() = BSESingletEnergiesTest;
    orbWrite.BSESingletCoefficients() = BSESingletCoefficientsTest;
    orbWrite.BSESingletCoefficientsAR() = BSESingletCoefficientsARTest;
    orbWrite.TransitionDipoles() = transitionDipolesTest;
    orbWrite.BSETripletEnergies() = BSETripletEnergiesTest;
    orbWrite.BSETripletCoefficients() = BSETripletCoefficientsTest;

    orbWrite.WriteToCpt(cpf, "Test Orbital");


    // Read Orbitals
    Orbitals orbRead;

    orbRead.ReadFromCpt(cpf, "Test Orbital");

    std::cout << orbRead.getQMpackage() << std::endl;
    std::cout << orbRead.getBasisSetSize() << std::endl;
    for (auto const& v:orbRead.TransitionDipoles()){
        std::cout << "(" << v.getX() << ", "
                  << v.getY() << ", " << v.getZ() << ")"
                  << std::endl;
    }

    BOOST_AUTO_TEST_SUITE_END()}
