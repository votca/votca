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
#include <cassert>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/checkpointwriter.h>
#include <votca/xtp/checkpointreader.h>
#include <votca/xtp/checkpoint.h>

BOOST_AUTO_TEST_SUITE(test_hdf5)
using namespace votca::xtp;
BOOST_AUTO_TEST_CASE(checkpoint_file_test) {

    int basisSetSize = 17;
    int occupiedLevels = 4;
    int unoccupiedLevels = 13;
    int numElectrons = 12;


    Eigen::VectorXd moeTest = Eigen::VectorXd::Random(17);
    Eigen::MatrixXd mocTest = Eigen::MatrixXd::Random(17,17);

    std::vector<QMAtom> atomsTest(1000);

    double qmEnergy = -2.1025e-3;

    std::string qmPackage = "NOPE";
    double selfEnergy = 3.14159e23;

    std::string dftBasis = "AWESOME basis*,, 2/.8";
    std::string auxBasis = "cos(theta) = pretty okay basis";

    int rpaMin = '?';
    int rpaMax = 1e3;

    unsigned int bseVmin = -6019386;
    unsigned int bseCmax = 42;

    double scaHfx = 3.14159;

    bool useTDA = true;


    Eigen::MatrixXd vxcTest = Eigen::MatrixXd::Random(200,200);
    std::string someECP = "aye aye Cap'n";

    Eigen::MatrixXd QPpertEnergiesTest = Eigen::MatrixXd::Random(31, 42);
    Eigen::MatrixXd QPdiagEnergiesTest = Eigen::VectorXd::Random(21);
    Eigen::MatrixXd QPdiagCoefficientsTest = Eigen::MatrixXd::Identity(42, 42);


    MatrixXfd eh_dTest = MatrixXfd::Random(32, 290);
    MatrixXfd eh_xTest = MatrixXfd::Random(3, 22);
    VectorXfd BSESingletEnergiesTest = VectorXfd::Random(25);
    MatrixXfd BSESingletCoefficientsTest = MatrixXfd::Random(25, 38);
    MatrixXfd BSESingletCoefficientsARTest = MatrixXfd::Random(42, 42);

    VectorXfd BSETripletEnergiesTest = VectorXfd::Random(33);
    MatrixXfd BSETripletCoefficientsTest = MatrixXfd::Random(33,31);


    std::vector <votca::tools::vec> transitionDipolesTest;
    for (size_t i =0; i < 1000; ++i){
        transitionDipolesTest.push_back(votca::tools::vec(1,2,3));
    }

    {
        // Write orbitals
        Orbitals orbWrite;

        orbWrite.setBasisSetSize(basisSetSize);
        orbWrite.setNumberOfLevels(occupiedLevels, unoccupiedLevels);
        orbWrite.setNumberOfElectrons(numElectrons);
        orbWrite.MOEnergies() = moeTest;
        orbWrite.MOCoefficients() = mocTest;

        for (auto const& qma:atomsTest){
            orbWrite.AddAtom(qma);
        }

        orbWrite.setQMEnergy(qmEnergy);
        orbWrite.setQMpackage(qmPackage);
        orbWrite.setSelfEnergy(selfEnergy);
        orbWrite.setDFTbasis(dftBasis);
        orbWrite.setAuxbasis(auxBasis);
        orbWrite.setRPAindices(rpaMin, rpaMax);
        // no need to write qpmin, qpmax
        orbWrite.setBSEindices(bseVmin, bseCmax, 3);
        orbWrite.setScaHFX(scaHfx);
        orbWrite.setTDAApprox(useTDA);
        orbWrite.setECP(someECP);
        orbWrite.QPpertEnergies() = QPpertEnergiesTest;
        orbWrite.QPdiagEnergies() = QPdiagEnergiesTest;
        orbWrite.QPdiagCoefficients() = QPdiagCoefficientsTest;
        orbWrite.eh_t() = eh_dTest;
        orbWrite.eh_s() = eh_xTest;
        orbWrite.BSESingletEnergies() = BSESingletEnergiesTest;
        orbWrite.BSESingletCoefficients() = BSESingletCoefficientsTest;
        orbWrite.BSESingletCoefficientsAR() = BSESingletCoefficientsARTest;
        orbWrite.TransitionDipoles() = transitionDipolesTest;
        orbWrite.BSETripletEnergies() = BSETripletEnergiesTest;
        orbWrite.BSETripletCoefficients() = BSETripletCoefficientsTest;

        orbWrite.WriteToCpt("xtp_testing.hdf5");

    }
    // Read Orbitals
    Orbitals orbRead;
    orbRead.ReadFromCpt("xtp_testing.hdf5");

    double tol = 1e-6;

    // Test the read values
    BOOST_CHECK_EQUAL(orbRead.getBasisSetSize() , basisSetSize);
    BOOST_CHECK_EQUAL(orbRead.getNumberOfLevels() , occupiedLevels + unoccupiedLevels);
    BOOST_CHECK_EQUAL(orbRead.getNumberOfElectrons() , numElectrons);
    BOOST_CHECK(orbRead.MOEnergies().isApprox(moeTest, tol));

    BOOST_CHECK(orbRead.MOCoefficients().isApprox(mocTest, tol));
    BOOST_CHECK_CLOSE(orbRead.getQMEnergy(), qmEnergy, tol);
    BOOST_CHECK_EQUAL(orbRead.getQMpackage(), qmPackage);
    BOOST_CHECK_CLOSE(orbRead.getSelfEnergy(), selfEnergy, tol);
    BOOST_CHECK_EQUAL(orbRead.getDFTbasis(), dftBasis);
    BOOST_CHECK_EQUAL(orbRead.getAuxbasis(), auxBasis);
    BOOST_CHECK_EQUAL(orbRead.getRPAmin(), rpaMin);
    BOOST_CHECK_EQUAL(orbRead.getRPAmax(), rpaMax);

    BOOST_CHECK_EQUAL(orbRead.getBSEvmin(), bseVmin);
    BOOST_CHECK_EQUAL(orbRead.getBSEcmax(), bseCmax);

    BOOST_CHECK_CLOSE(orbRead.getScaHFX(), scaHfx, tol);
    BOOST_CHECK_EQUAL(orbRead.getTDAApprox(), useTDA);
    BOOST_CHECK_EQUAL(orbRead.getECP(), someECP);
    BOOST_CHECK(orbRead.QPpertEnergies().isApprox(QPpertEnergiesTest, tol));
    BOOST_CHECK(orbRead.QPdiagEnergies().isApprox(QPdiagEnergiesTest, tol));
    BOOST_CHECK(orbRead.QPdiagCoefficients().isApprox(QPdiagCoefficientsTest));
    BOOST_CHECK(orbRead.eh_t().isApprox(eh_dTest, tol));
    BOOST_CHECK(orbRead.eh_s().isApprox(eh_xTest, tol));
    BOOST_CHECK(orbRead.BSESingletEnergies().isApprox(BSESingletEnergiesTest, tol));
    BOOST_CHECK(orbRead.BSESingletCoefficients().isApprox(BSESingletCoefficientsTest, tol));
    BOOST_CHECK(orbRead.BSESingletCoefficientsAR().isApprox(BSESingletCoefficientsARTest, tol));
    BOOST_CHECK(orbRead.BSETripletEnergies().isApprox(BSETripletEnergiesTest, tol));
    BOOST_CHECK(orbRead.BSETripletCoefficients().isApprox(BSETripletCoefficientsTest, tol));

    BOOST_REQUIRE_EQUAL(orbRead.TransitionDipoles().size(), transitionDipolesTest.size());

    for (size_t c = 0; c<transitionDipolesTest.size(); ++c){
        BOOST_CHECK(
            orbRead.TransitionDipoles()[c].isClose(transitionDipolesTest[c], tol));

    }

    BOOST_REQUIRE_EQUAL(orbRead.QMAtoms().size(), atomsTest.size());

    for (size_t i = 0; i<atomsTest.size(); ++i){
        auto atomRead = *(orbRead.QMAtoms()[i]);
        auto atomTest = atomsTest[i];
        BOOST_CHECK_EQUAL(atomRead.getAtomID(), atomTest.getAtomID());
        BOOST_CHECK(atomRead.getPos().isClose(atomTest.getPos(), tol));
        BOOST_CHECK_EQUAL(atomRead.getNuccharge(), atomTest.getNuccharge());
        BOOST_CHECK_EQUAL(atomRead.getPartialcharge(), atomTest.getPartialcharge());
        // no way to get qmatom index
    }
}

BOOST_AUTO_TEST_CASE(open_file_error){
    BOOST_REQUIRE_THROW(CheckpointFile cpf("/bin/mr/root/man.pls",
                                           CheckpointAccessLevel::READ),
                        std::runtime_error);
}

BOOST_AUTO_TEST_CASE(checkpoint_open_non_existing_loc) {
    CheckpointFile cpf ("testin_yo.ab", CheckpointAccessLevel::MODIFY);
    BOOST_REQUIRE_THROW(CheckpointReader r = cpf.getReader("/some/bulshit"),
                        std::runtime_error);

}

BOOST_AUTO_TEST_CASE(read_non_exisiting_matrix){

    CheckpointFile cpf("xtp_testing.hdf5", CheckpointAccessLevel::READ);
    CheckpointReader r = cpf.getReader("/QMdata");

    Eigen::MatrixXd someMatrix;

    BOOST_REQUIRE_THROW(r(someMatrix, "someMatrix012'5915.jb"),
                        std::runtime_error);
}

BOOST_AUTO_TEST_CASE(read_non_existing_scalar){
    CheckpointFile cpf("xtp_testing.hdf5", CheckpointAccessLevel::READ);
    CheckpointReader r = cpf.getReader("/QMdata");

    float someThing = 0;
    BOOST_REQUIRE_THROW(r(someThing, "someThing"), std::runtime_error);

}

BOOST_AUTO_TEST_SUITE_END()
